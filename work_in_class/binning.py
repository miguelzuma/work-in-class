#!/usr/bin/python

import numpy as np
import os
from classy import Class
import sys
import scipy.integrate as integrate
import wicmath as wicm
from wicmath import fit_pade
from wicmath import fit, Taylor
from scipy.interpolate import interp1d

class Binning():
    def __init__(self, fname, outdir='./'):
        self._cosmo = Class()
        # self._zbins = zbins
        # self._abins = abins
        self._set_file_names(fname, outdir)
        self._set_default_values()

    def _set_file_names(self, fname, outdir):
        """
        Set the file names and add a number if exists
        """

        fwzname = os.path.join(outdir, fname+'-wz-bins')
        fwaname = os.path.join(outdir, fname+'-wa-bins')
        fparamsname = os.path.join(outdir, fname+'-params')
        fshootname = os.path.join(outdir, fname+'-shooting')
        fPadename = os.path.join(outdir, fname+'-Pade')
        fFitname = os.path.join(outdir, fname+'-fit-F')

        i = 0

        while (os.path.exists(fwzname + '-%s.txt' % i) or
              os.path.exists(fwaname + '-%s.txt' % i) or
              os.path.exists(fparamsname + '-%s.txt' % i) or
              os.path.exists(fshootname + '-%s.txt' % i) or
              os.path.exists(fPadename + '-%s.txt' % i) or
              os.path.exists(fFitname + '-%s.txt' % i)):
            i += 1

        self._fwzname = fwzname + '-%s.txt' % i
        self._fwaname = fwaname + '-%s.txt' % i
        self._fparamsname = fparamsname + '-%s.txt' % i
        self._fshootname= fshootname + '-%s.txt' % i
        self._fPadename= fPadename + '-%s.txt' % i
        self._fFitname= fFitname + '-%s.txt' % i

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
                        "Omega_smg": -1}

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

    def set_n_coeffs_fit(self, n_coeffs, fix_origin=False):
        """
        Set the number of coefficients for Taylor fit of F(a) ~ \int dlna w
        """
        self.reset()
        self._n_coeffs = n_coeffs
        self._fix_origin = fix_origin
        self._binType = 'fit'
        self._params.update({'output': 'mPk', 'z_max_pk': 1000})

    def set_fit(self, fit_function, n_coeffs, variable_to_fit, fit_function_label=''):
        """
        Set the fitting_function and the number of coefficients.

        fit_function_label will be written in the header of fit files.
        """
        self.reset()
        self._n_coeffs = n_coeffs
        self._fix_origin = False
        self._variable_to_fit = variable_to_fit
        self._binType = 'fit'

    def set_bins(self, zbins, abins):
        """
        Set what bins to use and reset to avoid confusions.
        """
        self.reset()
        self._zbins = zbins
        self._abins = abins
        self._binType = 'bins'

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


    def compute_f_coefficients(self, params):
        """
        Returns the coefficients of the polynomial fit of f(a) = \int w and the
        maximum and minimum residual in absolute value.
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

        #####
        # NOTE: Emilio's functions from classy.pyx. self should not be used like this, but let's do an exception.
        #####

        def growthrate_at_k_and_z(cosmo, k, z):
            dz = 0.005
            f=-0.5*(1+z)*(-3.*np.log(cosmo.pk(k,z))+4.*np.log(cosmo.pk(k,z+dz))-np.log(cosmo.pk(k,z+2.*dz)))/(2.*dz)
            if (f==0.):
                raise CosmoSevereError(
                    "The step in redshift is too small to compute growth rate")
            return f

        def growthrate_at_z(cosmo, z):
            k_fid = 0.01
            return growthrate_at_k_and_z(cosmo, k_fid, z)
        #########

        def fprime(OmegaDEw, OmegaM):
            """
            Return value for df/dz.

            w, OmegaDE and OmegaM must be functions of z
            """
            def wrap(z, f):
                try:
                    #output = 1.5 * OmegaM(z) - 0.5 * (1 - 3*w(z)*OmegaDE(z)) * f - f**2
                    output = ((0.5 * (1 - 3*OmegaDEw(z))) * f + f**2 - 1.5 * OmegaM(z)) / (1+z)
                except Exception as e:
                    raise e
                return output
            return wrap


        # Compute the exact growth rate
        #####################
        z_max_pk = self._params['z_max_pk']
        z, w = b['z'], b['w_smg']
        rhoDE = b['(.)rho_smg']
        rhoM = (b['(.)rho_b'] + b['(.)rho_cdm'])
        rhoR = (b['(.)rho_g'] + b['(.)rho_ur'])
        DA = b['ang.diam.dist.']

        OmegaDEwF_exact = interp1d(z[z<=z_max_pk], (b['(.)rho_smg']/b['(.)rho_crit']*w)[z<=z_max_pk])
        OmegaMF = interp1d(z[z<=z_max_pk], (rhoM/b['(.)rho_crit'])[z<=z_max_pk])

        time_boundaries = [z[z<=z_max_pk][0], z[z<=z_max_pk][-1]]

        # Use LSODA integrator as some solutions were wrong with RK45 and OK
        # with this.
        f = integrate.solve_ivp(fprime(OmegaDEwF_exact, OmegaMF), time_boundaries, [growthrate_at_z(self._cosmo, z_max_pk)],
                                method='LSODA', dense_output=True)

        # Compute the exact \int dlna a
        ###############################
        Fint = []
        lna = -np.log(1+z)[::-1]
        for i in lna:
            Fint.append(integrate.trapz(w[::-1][lna>=i], lna[lna>=i]))
        Fint = np.array(Fint)

        #####

        zlim = 1000
        X = np.log(z + 1)[z<=zlim]

        #####################
        zTMP = z[z<=zlim]
        Y1 = Fint[::-1][z<zlim]  # Ordered as in CLASS

        #####################

        # Fit to fit_function
        #####################
        if not self._fix_origin:
            popt1, yfit1 = fit(Taylor, X, Y1, self._n_coeffs)
        else:
            def Taylor_fix_origin(x, c):
                return Taylor(x, np.concatenate([[0.], c]))
            popt1, yfit1 = fit(Taylor_fix_origin, X, Y1, self._n_coeffs)

        # Compute D_A for fitted model
        ################
        rhoDE_fit = b['(.)rho_smg'][-1]*np.exp(-3 * yfit1) *(1+zTMP)**3   ###### CHANGE WITH CHANGE OF FITTED THING

        H_fit = np.sqrt(rhoM[z<=zlim] + rhoR[z<=zlim] + rhoDE_fit)

        DA_fit = []
        for i in z[z<zlim]:
            #DA_fit.append(1/(1+i)*integrate.trapz(1/H_fit[zTMP<=i][::-1], zTMP[zTMP<=i][::-1]))
            DA_fit.append(1/(1+i)*integrate.simps(1/H_fit[zTMP<=i][::-1], zTMP[zTMP<=i][::-1], even='last'))
        DA_fit = np.array(DA_fit)

        # Compute the growth rate for fitted model
        ###############
        Xw_fit, w_fit = wicm.diff(X, yfit1)
        w_fit = -interp1d(Xw_fit, w_fit, bounds_error=False, fill_value='extrapolate')(X)

        OmegaMF_fit = interp1d(zTMP, 1-rhoDE_fit/H_fit**2-rhoR[z<=zlim]/H_fit**2)   ####### THIS FITS OBSERVABLES CORRECTLY
        #OmegaMF_fit = interp1d(zTMP, rhoM[z<=zlim]/H_fit**2)      ####### THIS FITS OBSERVABLES CORRECTLY
        OmegaDEwF_fit = interp1d(zTMP, rhoDE_fit/H_fit**2 * w_fit)

        f_fit = integrate.solve_ivp(fprime(OmegaDEwF_fit, OmegaMF_fit), [zTMP[0], zTMP[-1]], [growthrate_at_z(self._cosmo, zTMP[0])],
                                    method='LSODA', dense_output=True)

        # Free structures
        ###############
        self._cosmo.struct_cleanup()
        self._cosmo.empty()

        # Obtain rel. deviations.
        ################

        # Remove close to 0 points as rel.dev diverges. z = 0.05 is the lowest
        # redshift observed and is done in BOSS survey. arXiv: 1308.4164
        # DA_reldev = max(np.abs(DA_fit[zTMP>=0.04]/DA[ (z>=0.04) & (z<=zlim)] - 1))
        DA_reldev = max(np.abs(DA_fit/DA[z<=zlim] - 1))
        f_reldev = max(np.abs(f_fit.sol(zTMP)[0]/f.sol(zTMP)[0] - 1))

        return np.concatenate([popt1, [DA_reldev, f_reldev]]), shoot

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

    def compute_f_from_params(self, params_func, number_of_rows):
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
            # params_tmp = params_func().copy()

            try:
                wbins_tmp, shoot_tmp = self.compute_f_coefficients(params_func())
                wbins.append(wbins_tmp)
                params.append(self._params.copy())
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
        if self._variable_to_fit == 'f':
            self._params.update({'output': 'mPk', 'z_max_pk': 1000})
            fit_variable_function = self.compute_f_coefficients
        elif self._variable_to_fit == 'w':
            self.compute_w_coefficients,

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
            print "Bins already computed. Use reset if you want to compute it again"
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
                f.write('# ' + "Taylor fit for temporal variable ln(a) of \int dlna w\n")
                if self._fix_origin:
                    coeff_header_num = ['c_{}'.format(n) for n in range(1, self._n_coeffs + 1)]
                else:
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
