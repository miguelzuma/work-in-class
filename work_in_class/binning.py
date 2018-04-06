#!/usr/bin/python

import numpy as np
import os
from classy import Class
import sys
from scipy.optimize import curve_fit

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

    def set_PadeOrder(self, n_num, m_den):
        """
        Set what bins to use and reset to avoid confusions.
        """
        self._PadeOrder = [n_num, m_den]
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

        """
        Fit padde wrapper functions by Emilio Bellini.
        """

        def wrapper(x, n_num, n_den, *args):
            coeff_num = list(args[0:n_num+1])
            coeff_den = list(args[n_num+1:n_num+n_den+1])
            return pade_approx(x, coeff_num, coeff_den)

        def pade_approx(x, coeff_num, coeff_den):  # the actual fit function
            num = 0.
            den = 1.
            for i, coeff in enumerate(coeff_num):
                num = num + coeff*(x**i)
            for i, coeff in enumerate(coeff_den):
                den = den + coeff*(x**(i+1))
            return num/den

        def fit_pade(xdata, ydata, n_num, n_den):
            p0 = np.ones(n_num+n_den+1)
            popt, _ = curve_fit(lambda x, *p0: wrapper(xdata, n_num, n_den, *p0), xdata, ydata, p0=p0)
            yfit = wrapper(xdata, n_num, n_den, *popt)
            return popt, yfit

        """
        End of fit Pade wrapper.
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

        abins = 1./(b['z']+1)
        w = b['w_smg']

        padeCoefficients, padeFit = fit_pade(abins, w, *self._PadeOrder)

        r = np.abs(padeFit/w - 1.)

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
                f.write('# ' + "Pade coefficients up to {}th order (num) | {}th order (den) | min(residual) | max(residual)".format(*self._PadeOrder) + '\n')


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