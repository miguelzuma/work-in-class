#!/usr/bin/python

import numpy as np
import os
from classy import Class
import sys
from wicmath import fit_pade

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

        i = 0

        while (os.path.exists(fwzname + '-%s.txt' % i) or
              os.path.exists(fwaname + '-%s.txt' % i) or
              os.path.exists(fparamsname + '-%s.txt' % i) or
              os.path.exists(fshootname + '-%s.txt' % i) or
              os.path.exists(fPadename + '-%s.txt' % i)):
            i += 1

        self._fwzname = fwzname + '-%s.txt' % i
        self._fwaname = fwaname + '-%s.txt' % i
        self._fparamsname = fparamsname + '-%s.txt' % i
        self._fshootname= fshootname + '-%s.txt' % i
        self._fPadename= fPadename + '-%s.txt' % i

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

    def set_Pade(self, n_num, m_den, xvar='a', xReverse=False, accuracy=1e-3, increase=False, maxfev=0):
        """
        Set what Pade polynomial orders, temporal variable and its ordering use.
        """
        self._PadeOrder = [n_num, m_den]
        self._Pade_xvar = xvar
        self._Pade_xReverse = xReverse
        self._Pade_maxfev = maxfev
        self._Pade_increase = increase
        self._Pade_accuracy = accuracy
        self.reset()

    def set_bins(self, zbins, abins):
        """
        Set what bins to use and reset to avoid confusions.
        """
        self._zbins = zbins
        self._abins = abins
        self.reset()

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
        self._create_output_files(True)

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
                self._save_computed(params, shoot, wbins, True)

                params = []
                wbins = []
                shoot = []

        self._save_computed(params, shoot, wbins, True)

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

    def _create_output_files(self, Pade=False):
        """
        Initialize the output files.
        """
        # TODO: Add check if files exist
        with open(self._fparamsname, 'a') as f:
            f.write('# ' + "Dictionary of params to use with cosmo.set()" + '\n')

        with open(self._fshootname, 'a') as f:
            f.write('# ' + "Shooting variable value" + '\n')

        if not Pade:
            with open(self._fwzname, 'a') as f:
                f.write('# ' + "Bins on redshift" + '\n')
                f.write('# ' + str(self._zbins).strip('[]').replace('\n', '') + '\n')

            with open(self._fwaname, 'a') as f:
                f.write('# ' + "Bins on scale factor" + '\n')
                f.write('# ' + str(self._abins).strip('[]').replace('\n', '') + '\n')
        else:
            with open(self._fPadename, 'a') as f:
                f.write('# ' + "Pade fit for temporal variable {} \n".format(self._Pade_xvar))
                coeff_header_num = ['num_{}'.format(n) for n in range(self._PadeOrder[0] + 1)]
                coeff_header_den = ['den_{}'.format(n + 1) for n in range(self._PadeOrder[1])]
                res_header = ['min(residual)', 'max(residual)']
                f.write('# ' + ' '.join(coeff_header_num + coeff_header_den + res_header) + '\n')

    def _save_computed(self, params, shoot, wbins, Pade=False):
        """
        Save stored iterations in file.
        """
        with open(self._fparamsname, 'a') as f:
            for i in params:
                f.write(str(i)+'\n')

        with open(self._fshootname, 'a') as f:
            np.savetxt(f, shoot)

        if not Pade:
            wzbins, wabins = wbins
            with open(self._fwzname, 'a') as f:
                np.savetxt(f, wzbins)

            with open(self._fwaname, 'a') as f:
                np.savetxt(f, wabins)
        else:
            with open(self._fPadename, 'a') as f:
                np.savetxt(f, wbins)

    def reset(self):
        """
        Reset class
        """
        self._cosmo.struct_cleanup()
        self._cosmo.empty()
        self._set_default_values()
