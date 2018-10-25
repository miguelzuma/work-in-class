#!/usr/bin/python

import numpy as np
import os
from classy import Class
from classy import CosmoSevereError
import sys
import scipy.integrate as integrate
import wicmath as wicm
from wicmath import fit_pade
from scipy.interpolate import interp1d
import cosmo_extra

class Binning():
    def __init__(self, fname, outdir='./'):
        self._cosmo = Class()
        self._fname = fname
        self._outdir = outdir
        self._set_default_values()

    def _set_full_filenames(self, filesuffixes):
        """
        Return a list with the fullnames of the others, increasing their number
        in case they exist

        Additionally, set the full file names, of self._fparamsname and
        self.fshootname.
        """
        fullfilenames = []
        for suffix in filesuffixes + ['params', 'shooting']:
            fullfilenames.append(os.path.join(self._outdir, self._fname + '-' + suffix))

        i = 0

        bools = [True] * len(fullfilenames)

        while 1:
            for n, f in enumerate(fullfilenames):
                bools[n] = os.path.exists(f + '-%s.txt' % i)

            if True not in bools:
                break

            i += 1

        self._fparamsname = fullfilenames[-2] + '-%s.txt' % i
        self._fshootname = fullfilenames[-1] + '-%s.txt' % i

        return [f + '-%s.txt' % i for f in fullfilenames[:-2]]

    def _set_default_values(self):
        """
        Set default values of parameters lists.
        """
        self.params_smg = []
        self.params_2_smg = []
        self.h = []
        self.Omega_cdm = []

        self.gravity_model = []

        self._params = {"Omega_Lambda": 0,
                        "Omega_fld": 0,
                        "Omega_smg": -1,
                        'output': 'mPk',   #
                        'z_max_pk': 1000}  # Added for relative errors in f.

        self._computed = False
        self._path = []
        self._binType = ''

    def set_Pade(self, n_num, m_den, xvar='a', xReverse=False, accuracy=1e-3, increase=False, maxfev=0):
        """
        Set what Pade polynomial orders, temporal variable and its ordering use.
        """
        self.reset()
        self._PadeOrder = [n_num, m_den]
        self._Pade_xvar = xvar
        self._Pade_xReverse = xReverse
        self._Pade_maxfev = maxfev
        self._Pade_increase = increase
        self._Pade_accuracy = accuracy
        self._binType = 'Pade'

    def set_fit(self, fit_function, n_coeffs, variable_to_fit, fit_function_label='', z_max_pk=1000, bounds=(-np.inf, np.inf), p0=[], xvar='ln(1+z)'):
        """
        Set the fitting_function and the number of coefficients.

        variable_to_fit must be one of 'F'or 'w'.

        fit_function_label will be written in the header of fit files.
        """
        self.reset()
        self._fit_function = fit_function
        self._n_coeffs = n_coeffs
        self._list_variables_to_fit = ['F', 'w', 'logRho', 'logX', 'X']
        if variable_to_fit in self._list_variables_to_fit:
            self._variable_to_fit = variable_to_fit
        else:
            raise ValueError('variable_to_fit must be one of {}'.format(self._list_variables_to_fit))

        self._fit_function_label = fit_function_label
        self._binType = 'fit'
        self._fFitname = self._set_full_filenames(['fit-' + variable_to_fit])[0]
        self._params.update({'z_max_pk': z_max_pk})
        self._fit_bounds = bounds
        self._p0 = p0
        list_fit_xvar = ['ln(1+z)', 'lna', 'a', '(1-a)']
        if xvar in list_fit_xvar:
            self._fit_xvar = xvar
        else:
            raise ValueError('xvar must be one of {}'.format(list_fit_xvar))

    def set_bins(self, zbins, abins):
        """
        Set what bins to use and reset to avoid confusions.
        """
        self.reset()
        self._zbins = zbins
        self._abins = abins
        self._binType = 'bins'
        self._fwzname, self._fwaname = self._set_full_filenames(['wz-bins', 'wa-bins'] )

    def _read_from_file(self, path):
        """
        Return params for class from files used in quintessence Marsh.
        """
        with open(path) as f:
            f.readline()
            header = f.readline().split()[3:]  # Remove '#', 'w0', and 'wa'

        columns = np.loadtxt(path, unpack=True)[2:]  # Remove columns w0, wa

        for index_h, head in enumerate(header):
            if head[-1] == 'h':
                break

        self.params_smg = zip(*columns[:index_h])
        self.params_2_smg = [list(row[~np.isnan(row)]) for row in
                             np.array(zip(*columns[index_h+2:]))]
        self.h = columns[index_h]
        self.Omega_cdm = columns[index_h+1]

        self.gravity_model = os.path.basename(path).split('-')[0]

    def _params_from_row(self, row):
        """
        Set parameters.
        """
        params = self._params
        params.update({
            'parameters_smg': str(self.params_smg[row]).strip('[()]'),
            'h': self.h[row],
            'Omega_cdm': self.Omega_cdm[row],
            'gravity_model': self.gravity_model
        })

        if len(self.params_2_smg):
            params.update({
                'parameters_2_smg': str(self.params_2_smg[row]).strip('[()]')
            })

        return params

    def compute_bins(self, params):
        """
        Compute the w_i bins for the model with params.
        """
        wzbins = np.empty(len(self._zbins))
        wabins = np.empty(len(self._abins))
        self._params = params
        self._cosmo.set(params)
        try:
            self._cosmo.compute()
            for n, z in enumerate(self._zbins):
                wzbins[n] = self._cosmo.w_smg(z)
            for n, a in enumerate(self._abins):
                wabins[n] = self._cosmo.w_smg(1./a-1.)
            shoot = self._cosmo.get_current_derived_parameters(['tuning_parameter'])['tuning_parameter']
        except Exception as e:
            self._cosmo.struct_cleanup()
            self._cosmo.empty()
            raise e

        self._cosmo.struct_cleanup()
        self._cosmo.empty()

        return wzbins, wabins, shoot

    def _compute_common_init(self, params):
        """
        Common first steps for compute methods
        """
        self._params.update(params)
        self._cosmo.set(self._params)

        try:
            self._cosmo.compute()
            b = self._cosmo.get_background()
            shoot = self._cosmo.get_current_derived_parameters(['tuning_parameter'])['tuning_parameter']
        except Exception as e:
            self._cosmo.struct_cleanup()
            self._cosmo.empty()
            raise e

        return b, shoot

    def _fit(self, X, Y):
        """
        Fits self._fit_function to X, Y, with self._n_coeffs.
        """
        return wicm.fit(self._fit_function, X, Y, self._n_coeffs, bounds=self._fit_bounds, p0=self._p0)

    def _get_fit_xvar_and_log1Plusz(self, z):
        """
        Return the X array to fit and log(1+z)
        """
        if self._fit_xvar == 'ln(1+z)':
            X = np.log(z + 1)
            lna = X
        elif self._fit_xvar == 'lna':
            X = - np.log(z + 1)
            lna = - X
        elif self._fit_xvar == 'a':
            X = 1. / (1 + z)
            lna = np.log(z + 1)
        elif self._fit_xvar == '(1-a)':
            X = (z / (1 + z))
            lna = np.log(z + 1)

        return X, lna

    def compute_fit_coefficients_for_F(self, params):
        """
        Returns the coefficients of the polynomial fit of f(a) = \int w and the
        maximum and minimum residual in absolute value.
        """
        b, shoot = self._compute_common_init(params)

        # Compute the exact \int dlna a
        ###############################
        z, w = b['z'], b['w_smg']

        Fint = []
        lna = -np.log(1+z)[::-1]
        for i in lna:
            Fint.append(integrate.trapz(w[::-1][lna>=i], lna[lna>=i]))
        Fint = np.array(Fint)

        #####

        zlim = self._params['z_max_pk']
        # Note that lna is log(1+z). I used this name because is more convenient
        X, lna = self._get_fit_xvar_and_log1Plusz(z[z <= zlim])

        #####################
        zTMP = z[z<=zlim]
        Y1 = Fint[::-1][z<zlim]  # Ordered as in CLASS

        #####################

        # Fit to fit_function
        #####################
        popt1, yfit1 = self._fit(X, Y1)

        # Obtain max. rel. dev. for DA and f.
        #####################

        rhoDE_fit = b['(.)rho_smg'][-1]*np.exp(-3 * yfit1) *(1+zTMP)**3   ###### CHANGE WITH CHANGE OF FITTED THING

        Xw_fit, w_fit = wicm.diff(lna, yfit1)
        w_fit = -interp1d(Xw_fit, w_fit, bounds_error=False, fill_value='extrapolate')(lna)

        DA_reldev, f_reldev = self._compute_maximum_relative_error_DA_f(rhoDE_fit, w_fit)

        # Free structures
        ###############
        self._cosmo.struct_cleanup()
        self._cosmo.empty()

        return np.concatenate([popt1, [DA_reldev, f_reldev]]), shoot

    def compute_fit_coefficients_for_logX(self, params):
        """
        Returns the coefficients of the polynomial fit of log(rho/rho_0) = -3
        \int dlna (w+1) and the maximum and minimum residual in absolute value.
        """
        b, shoot = self._compute_common_init(params)

        # Compute the exact -3 \int dlna (w + 1)
        ###############################
        z = b['z']

        logX = np.log(b['(.)rho_smg']/b['(.)rho_smg'][-1])

        #####

        zlim = self._params['z_max_pk']
        # Note that lna is log(1+z). I used this name because is more convenient
        X, lna = self._get_fit_xvar_and_log1Plusz(z[z <= zlim])

        #####################
        Y1 = logX[z <= zlim]

        #####################

        # Fit to fit_function
        #####################
        popt1, yfit1 = self._fit(X, Y1)

        # Obtain max. rel. dev. for DA and f.
        #####################

        rhoDE_fit = b['(.)rho_smg'][-1] * np.exp(yfit1)   ###### CHANGE WITH CHANGE OF FITTED THING

        Xw_fit, ThreewPlus1 = wicm.diff(lna, yfit1)
        w_fit = ThreewPlus1 / 3. - 1  # The minus sign is taken into account by the CLASS ordering.
        w_fit = interp1d(Xw_fit, w_fit, bounds_error=False, fill_value='extrapolate')(lna)

        DA_reldev, f_reldev = self._compute_maximum_relative_error_DA_f(rhoDE_fit, w_fit)

        # Free structures
        ###############
        self._cosmo.struct_cleanup()
        self._cosmo.empty()

        return np.concatenate([popt1, [DA_reldev, f_reldev]]), shoot

    def compute_fit_coefficients_for_X(self, params):
        """
        Returns the coefficients of the polynomial fit of rho/rho_0 = exp[-3
        \int dlna (w+1)] and the maximum and minimum residual in absolute value.
        """
        b, shoot = self._compute_common_init(params)

        # Compute the exact -3 \int dlna (w + 1)
        ###############################
        z = b['z']

        Y = b['(.)rho_smg']/b['(.)rho_smg'][-1]

        #####

        zlim = self._params['z_max_pk']
        # Note that lna is log(1+z). I used this name because is more convenient
        X, lna = self._get_fit_xvar_and_log1Plusz(z[z <= zlim])

        #####################
        Y1 = Y[z <= zlim]

        #####################

        # Fit to fit_function
        #####################
        popt1, yfit1 = self._fit(X, Y1)

        # Obtain max. rel. dev. for DA and f.
        #####################

        rhoDE_fit = b['(.)rho_smg'][-1] * yfit1   ###### CHANGE WITH CHANGE OF FITTED THING

        Xw_fit, ThreewPlus1 = wicm.diff(lna, yfit1) / (yfit1[:-1] + np.diff(yfit1)/2.)
        w_fit = ThreewPlus1 / 3. - 1  # The minus sign is taken into account by the CLASS ordering.
        w_fit = interp1d(Xw_fit, w_fit, bounds_error=False, fill_value='extrapolate')(lna)

        DA_reldev, f_reldev = self._compute_maximum_relative_error_DA_f(rhoDE_fit, w_fit)

        # Free structures
        ###############
        self._cosmo.struct_cleanup()
        self._cosmo.empty()

        return np.concatenate([popt1, [DA_reldev, f_reldev]]), shoot

    def compute_fit_coefficients_for_logRho(self, params):
        """
        Returns the coefficients of the fit of ln(rho_de) and the maximum and
        minimum residual in absolute value.
        """
        b, shoot = self._compute_common_init(params)

        # Compute the exact log(rho)
        ###############################
        z = b['z']

        logX = np.log(b['(.)rho_smg'])

        #####

        zlim = self._params['z_max_pk']
        # Note that lna is log(1+z). I used this name because is more convenient
        X, lna = self._get_fit_xvar_and_log1Plusz(z[z <= zlim])

        #####################
        Y1 = logX[z <= zlim]

        #####################

        # Fit to fit_function
        #####################

        popt1, yfit1 = self._fit(X, Y1)

        # Obtain max. rel. dev. for DA and f.
        #####################

        rhoDE_fit = np.exp(yfit1)   ###### CHANGE WITH CHANGE OF FITTED THING

        Xw_fit, ThreewPlus1 = wicm.diff(lna, yfit1 - yfit1[-1])
        w_fit = ThreewPlus1 / 3. - 1  # The minus sign is taken into account by the CLASS ordering.
        w_fit = interp1d(Xw_fit, w_fit, bounds_error=False, fill_value='extrapolate')(lna)

        DA_reldev, f_reldev = self._compute_maximum_relative_error_DA_f(rhoDE_fit, w_fit)

        # Free structures
        ###############
        self._cosmo.struct_cleanup()
        self._cosmo.empty()

        return np.concatenate([popt1, [DA_reldev, f_reldev]]), shoot

    def compute_fit_coefficients_for_w(self, params):
        """
        Returns the coefficients of the polynomial fit of f(a) = \int w and the
        maximum and minimum residual in absolute value.
        """
        b, shoot = self._compute_common_init(params)

        # Compute the exact \int dlna a
        ###############################
        z, w = b['z'], b['w_smg']

        zlim = self._params['z_max_pk']
        # Note that lna is log(1+z). I used this name because is more convenient
        X, lna = self._get_fit_xvar_and_log1Plusz(z[z <= zlim])

        #####################
        zTMP = z[z <= zlim]
        Y1 = w[z <= zlim]

        # Fit to fit_function
        #####################
        popt1, yfit1 = self._fit(X, Y1)

        # Obtain max. rel. dev. for DA and f.
        #####################
        Fint = []
        lna = -np.log(1+zTMP)[::-1]
        for i in lna:
            Fint.append(integrate.trapz(yfit1[::-1][lna >= i], lna[lna >= i]))
        Fint = np.array(Fint[::-1])

        rhoDE_fit = b['(.)rho_smg'][-1] * np.exp(-3 * Fint) * (1+zTMP)**3   # CHANGE WITH CHANGE OF FITTED THING

        # TODO: needed?
        # Xw_fit, w_fit = X, yfit1
        # w_fit = interp1d(Xw_fit, w_fit, bounds_error=False, fill_value='extrapolate')(X)
        w_fit = yfit1

        DA_reldev, f_reldev = self._compute_maximum_relative_error_DA_f(rhoDE_fit, w_fit)

        # Free structures
        ###############
        self._cosmo.struct_cleanup()
        self._cosmo.empty()

        return np.concatenate([popt1, [DA_reldev, f_reldev]]), shoot

    def _compute_maximum_relative_error_DA_f(self, rhoDE_fit, w_fit):
        """
        Return the relative error for the diameter angular distance and the
        growth factor, f.

        rhoDE_fit = array
        wfit = interp1d(w)
        """

        b = self._cosmo.get_background()

        # Compute the exact growth rate
        #####################
        z_max_pk = self._params['z_max_pk']
        zlim = z_max_pk
        z, w = b['z'], b['w_smg']
        zTMP = z[z <= zlim]
        rhoM = (b['(.)rho_b'] + b['(.)rho_cdm'])
        rhoR = (b['(.)rho_g'] + b['(.)rho_ur'])
        DA = b['ang.diam.dist.']

        OmegaDEwF_exact = interp1d(z[z <= z_max_pk], (b['(.)rho_smg']/b['(.)rho_crit']*w)[z <= z_max_pk])
        OmegaMF = interp1d(z[z <= z_max_pk], (rhoM/b['(.)rho_crit'])[z <= z_max_pk])

        time_boundaries = [z[z <= z_max_pk][0], z[z <= z_max_pk][-1]]

        # Use LSODA integrator as some solutions were wrong with RK45 and OK
        # with this.
        f = integrate.solve_ivp(cosmo_extra.fprime(OmegaDEwF_exact, OmegaMF), time_boundaries, [cosmo_extra.growthrate_at_z(self._cosmo, z_max_pk)],
                                method='LSODA', dense_output=True)


        # Compute D_A for fitted model
        ################
        H_fit = np.sqrt(rhoM[z <= zlim] + rhoR[z <= zlim] + rhoDE_fit)

        DA_fit = cosmo_extra.angular_distance(z[z <= zlim], H_fit[zTMP <= zlim])


        # Compute the growth rate for fitted model
        ###############

        OmegaMF_fit = interp1d(zTMP, 1-rhoDE_fit/H_fit**2-rhoR[z <= zlim]/H_fit**2)   ####### THIS FITS OBSERVABLES CORRECTLY
        # OmegaMF_fit = interp1d(zTMP, rhoM[z<=zlim]/H_fit**2)      ####### THIS FITS OBSERVABLES CORRECTLY
        OmegaDEwF_fit = interp1d(zTMP, rhoDE_fit/H_fit**2 * w_fit)

        f_fit = integrate.solve_ivp(cosmo_extra.fprime(OmegaDEwF_fit, OmegaMF_fit), [zTMP[0], zTMP[-1]], [cosmo_extra.growthrate_at_z(self._cosmo, zTMP[0])],
                                    method='LSODA', dense_output=True)

        # Obtain rel. deviations.
        ################

        # Remove close to 0 points as rel.dev diverges. z = 0.05 is the lowest
        # redshift observed and is done in BOSS survey. arXiv: 1308.4164
        # DA_reldev = max(np.abs(DA_fit[zTMP>=0.04]/DA[ (z>=0.04) & (z<=zlim)] - 1))
        DA_reldev = max(np.abs(DA_fit/DA[z <= zlim] - 1))
        f_reldev = max(np.abs(f_fit.sol(zTMP)[0]/f.sol(zTMP)[0] - 1))

        return DA_reldev, f_reldev

    def compute_Pade_coefficients(self, params):
        """
        Returns the Pade coefficients for w computed from params and the maximum
        and minimum residual in absolute value.
        """
        self._params = params
        self._cosmo.set(params)

        try:
            self._cosmo.compute()
            b = self._cosmo.get_background()
            shoot = self._cosmo.get_current_derived_parameters(['tuning_parameter'])['tuning_parameter']
        except Exception as e:
            self._cosmo.struct_cleanup()
            self._cosmo.empty()
            raise e

        self._cosmo.struct_cleanup()
        self._cosmo.empty()

        xDict = {'z': b['z'],
                 'z+1': b['z']+1,
                 'a': 1./(b['z']+1),
                 'log(a)': -np.log(b['z']+1),
                 'log(z+1)': np.log(b['z']+1)}

        X = xDict[self._Pade_xvar]
        w = b['w_smg']

        if self._Pade_xReverse:
            X = X[::-1]
            w = w[::-1]

        PadeOrder = np.array(self._PadeOrder)

        if not self._Pade_increase:
            reduceOrder = [[0, 0], [1, 0], [0, 1], [2, 0], [2, 1], [3, 1]]
            orderList = PadeOrder - reduceOrder

        else:
            orderList = [[1, 1], [2, 0], [3, 0], [2, 1], [2, 2], [3, 1], [4, 0],
                         [2, 3], [3, 2], [4, 1], [5, 0], [3, 3], [4, 2], [5, 1],
                         [3, 4], [4, 3], [5, 2], [3, 5], [4, 4], [5, 3], [4, 5],
                         [5, 4], [5, 5]]

        r = np.array([np.inf])
        for order in orderList:
            # Increase order of Pade up to [5/5].
            try:
                padeCoefficientsTMP, padeFitTMP = fit_pade(X, w, *order,
                                                           maxfev=self._Pade_maxfev)
                rTMP = np.abs(padeFitTMP/w - 1.)
                if self._Pade_increase and (np.max(rTMP) > self._Pade_accuracy):
                    if np.max(rTMP) < np.max(r):
                        padeCoefficients = padeCoefficientsTMP
                        r = rTMP
                    continue
                else:
                    padeCoefficients = padeCoefficientsTMP
                    r = rTMP
                    break
            except Exception as e:
                if (order == orderList[-1]) and (len(r) == 1):
                    raise e

                continue

        zeros = (PadeOrder - order)

        numCoefficients = np.append(padeCoefficients[:order[0] + 1], [0.]*zeros[0])
        denCoefficients = np.append(padeCoefficients[order[0] + 1:], [0.]*zeros[1])
        padeCoefficients = np.concatenate([numCoefficients, denCoefficients])

        return np.concatenate([padeCoefficients, [np.min(r), np.max(r)]]), shoot

    def compute_bins_from_params(self, params_func, number_of_rows):
        """
        Compute the w_i bins for the models given by the function
        params_func iterated #iterations.
        """
        self._create_output_files()

        wzbins = []
        wabins = []
        params = []
        shoot = []

        for row in range(number_of_rows):
            sys.stdout.write("{}/{}\n".format(row+1, number_of_rows))
            params_tmp = params_func().copy()

            try:
                wzbins_tmp, wabins_tmp, shoot_tmp = self.compute_bins(params_tmp)
                wzbins.append(wzbins_tmp)
                wabins.append(wabins_tmp)
                params.append(params_tmp)
                shoot.append(shoot_tmp)
                # Easily generalizable. It could be inputted a list with the
                # desired derived parameters and store the whole dictionary.
            except Exception as e:
                sys.stderr.write(str(self._params) + '\n')
                sys.stderr.write(str(e))
                sys.stderr.write('\n')
                continue

            if len(wzbins) == 5:
                self._save_computed(params, shoot, [wzbins, wabins])

                params = []
                wzbins = []
                wabins = []
                shoot = []

        self._save_computed(params, shoot, [wzbins, wabins])

    def compute_Pade_from_params(self, params_func, number_of_rows):
        """
        Compute the w_i bins for the models given by the function
        params_func iterated #iterations.
        """
        self._create_output_files()

        wbins = []
        params = []
        shoot = []

        for row in range(number_of_rows):
            sys.stdout.write("{}/{}\n".format(row+1, number_of_rows))
            params_tmp = params_func().copy()

            try:
                wbins_tmp, shoot_tmp = self.compute_Pade_coefficients(params_tmp)
                wbins.append(wbins_tmp)
                params.append(params_tmp)
                shoot.append(shoot_tmp)
                # Easily generalizable. It could be inputted a list with the
                # desired derived parameters and store the whole dictionary.
            except Exception as e:
                sys.stderr.write(str(self._params) + '\n')
                sys.stderr.write(str(e))
                sys.stderr.write('\n')
                continue

            if len(wbins) == 5:
                self._save_computed(params, shoot, wbins)

                params = []
                wbins = []
                shoot = []

        self._save_computed(params, shoot, wbins)

    def compute_fit_from_params(self, params_func, number_of_rows):
        """
        Compute the fit for the models given by the function
        params_func iterated #iterations.

        The variable to fit is chosen in self.set_fit
        """
        # TODO: If this grows, consider creating a separate method
        if self._variable_to_fit == 'F':
            fit_variable_function = self.compute_fit_coefficients_for_F
        elif self._variable_to_fit == 'w':
            fit_variable_function = self.compute_fit_coefficients_for_w
        elif self._variable_to_fit == 'logRho':
            fit_variable_function = self.compute_fit_coefficients_for_logRho
        elif self._variable_to_fit == 'logX':
            fit_variable_function = self.compute_fit_coefficients_for_logX
        elif self._variable_to_fit == 'X':
            fit_variable_function = self.compute_fit_coefficients_for_X

        self._create_output_files()


        coeffs = []
        params = []
        shoot = []

        for row in range(number_of_rows):
            sys.stdout.write("{}/{}\n".format(row+1, number_of_rows))
            # params_tmp = params_func().copy()

            try:
                coeffs_tmp, shoot_tmp = fit_variable_function(params_func())
                coeffs.append(coeffs_tmp)
                params.append(self._params.copy())
                shoot.append(shoot_tmp)
                # Easily generalizable. It could be inputted a list with the
                # desired derived parameters and store the whole dictionary.
            except Exception as e:
                sys.stderr.write(str(self._params) + '\n')
                sys.stderr.write(str(e))
                sys.stderr.write('\n')
                continue

            if len(coeffs) == 5:
                self._save_computed(params, shoot, coeffs)

                params = []
                coeffs = []
                shoot = []

        self._save_computed(params, shoot, coeffs)

    def compute_bins_from_file(self, path):
        """
        Compute the w_i bins for the models given in path.
        """
        if self._computed is True:
            print("Bins already computed. Use reset if you want to compute it again")
            return

        self._path = path

        self._read_from_file(path)

        def params_gen(length):
            row = 0
            while row < length:
                yield self._params_from_row(row)
                row += 1

        params = params_gen(len(self.params_smg))

        self.compute_bins_from_params(params.next,
                                      len(self.params_smg))

    def _create_output_files(self):
        """
        Initialize the output files.
        """
        # TODO: Add check if files exist
        with open(self._fparamsname, 'a') as f:
            f.write('# ' + "Dictionary of params to use with cosmo.set()" + '\n')

        with open(self._fshootname, 'a') as f:
            f.write('# ' + "Shooting variable value" + '\n')

        if self._binType == 'bins':
            with open(self._fwzname, 'a') as f:
                f.write('# ' + "Bins on redshift" + '\n')
                f.write('# ' + str(self._zbins).strip('[]').replace('\n', '') + '\n')

            with open(self._fwaname, 'a') as f:
                f.write('# ' + "Bins on scale factor" + '\n')
                f.write('# ' + str(self._abins).strip('[]').replace('\n', '') + '\n')
        elif self._binType == 'Pade':
            with open(self._fPadename, 'a') as f:
                f.write('# ' + "Pade fit for temporal variable {} \n".format(self._Pade_xvar))
                coeff_header_num = ['num_{}'.format(n) for n in range(self._PadeOrder[0] + 1)]
                coeff_header_den = ['den_{}'.format(n + 1) for n in range(self._PadeOrder[1])]
                res_header = ['min(residual)', 'max(residual)']
                f.write('# ' + ' '.join(coeff_header_num + coeff_header_den + res_header) + '\n')
        elif self._binType == 'fit':
            with open(self._fFitname, 'a') as f:
                f.write('# ' + "{} fit for temporal variable {} of {}\n".format(self._fit_function_label, self._fit_xvar, self._variable_to_fit))
                coeff_header_num = ['c_{}'.format(n) for n in range(self._n_coeffs)]
                res_header = ['max(rel.dev. D_A)', 'max(rel.dev. f)']
                f.write('# ' + ' '.join(coeff_header_num + res_header) + '\n')

    def _save_computed(self, params, shoot, wbins):
        """
        Save stored iterations in file.
        """
        with open(self._fparamsname, 'a') as f:
            for i in params:
                f.write(str(i)+'\n')

        with open(self._fshootname, 'a') as f:
            np.savetxt(f, shoot)

        if self._binType == 'bins':
            wzbins, wabins = wbins
            with open(self._fwzname, 'a') as f:
                np.savetxt(f, wzbins)

            with open(self._fwaname, 'a') as f:
                np.savetxt(f, wabins)
        elif self._binType == 'Pade':
            with open(self._fPadename, 'a') as f:
                np.savetxt(f, wbins)
        elif self._binType == 'fit':
            with open(self._fFitname, 'a') as f:
                np.savetxt(f, wbins)

    def reset(self):
        """
        Reset class
        """
        self._cosmo.struct_cleanup()
        self._cosmo.empty()
        self._set_default_values()
