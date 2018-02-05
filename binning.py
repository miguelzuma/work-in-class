#!/usr/bin/python

import numpy as np
import os
from classy import Class
import sys


class Binning():
    def __init__(self, zbins, abins, fname, outdir='./'):
        self._cosmo = Class()
        self._zbins = zbins
        self._abins = abins
        self._fwzname = os.path.join(outdir, fname+'-wz-bins.txt')
        self._fwaname = os.path.join(outdir, fname+'-wa-bins.txt')
        self._fparamsname = os.path.join(outdir, fname+'-params.txt')
        self._set_default_values()

    def _set_empty_wbins(self):
        """
        Initialize wzbins and wabins arrays to ones
        """
        self.wzbins = np.ones((5, len(self._zbins)))  # Each 5 rows it is saved on disk
        self.wabins = np.ones((5, len(self._abins)))

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
        except Exception as e:
            self._cosmo.struct_cleanup()
            self._cosmo.empty()
            sys.stderr.write(str(self._params) + '\n')
            sys.stderr.write(str(e))
            sys.stderr.write('\n')
            raise Exception

        self._cosmo.struct_cleanup()
        self._cosmo.empty()

        return wzbins, wabins

    def compute_bins_from_params(self, params_func, number_of_rows):
        """
        Compute the w_i bins for the models given by the function
        params_func iterated #iterations.
        """
        self._create_output_files()

        wzbins = []
        wabins = []
        params = []

        for row in range(number_of_rows):
            sys.stdout.write("{}/{}\n".format(row+1, number_of_rows))
            params_tmp = params_func()

            try:
                wzbins_tmp, wabins_tmp = self.compute_bins(params_tmp)
                wzbins.append(wzbins_tmp)
                wabins.append(wabins_tmp)
                params.append(params_tmp)
            except Exception:
                continue

            if len(wzbins) == 5:
                self._save_computed(params, wzbins, wabins)

                params = []
                wzbins = []
                wabins = []

        self._save_computed(params, wzbins, wabins)

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

        with open(self._fwzname, 'a') as f:
            f.write('# ' + "Bins on redshift" + '\n')
            f.write('# ' + str(self._zbins).strip('[]').replace('\n', '') + '\n')

        with open(self._fwaname, 'a') as f:
            f.write('# ' + "Bins on scale factor" + '\n')
            f.write('# ' + str(self._abins).strip('[]').replace('\n', '') + '\n')

    def _save_computed(self, params, wzbins, wabins):
        """
        Save stored iterations in file.
        """
        if self._computed is True:
            print "Bins not yet computed. Compute them first"
            return

        with open(self._fparamsname, 'a') as f:
            for i in params:
                f.write(str(i)+'\n')

        with open(self._fwzname, 'a') as f:
            np.savetxt(f, wzbins)

        with open(self._fwaname, 'a') as f:
            np.savetxt(f, wabins)

    def reset(self):
        """
        Reset class
        """
        self._cosmo.struct_cleanup()
        self._cosmo.empty()
        self._set_default_values()
