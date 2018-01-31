#!/usr/bin/python

import numpy as np
import os
from classy import Class


class Binning():
    def __init__(self, zbins, abins, fname):
        self._cosmo = Class()
        self._zbins = zbins
        self._abins = abins
        self._fwzname = fname+'-wz-bins.txt'
        self._fwaname = fname+'-wa-bins.txt'
        self._fparamsname = fname+'-params.txt'
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

    def _set_params(self, row):
        """
        Set parameters.
        """
        self._params.update({
            'parameters_smg': str(self.params_smg[row]).strip('[()]'),
            'h': self.h[row],
            'Omega_cdm': self.Omega_cdm[row],
            'gravity_model': self.gravity_model
        })

        if len(self.params_2_smg):
            self._params.update({
                'parameters_2_smg': str(self.params_2_smg[row]).strip('[()]')
            })

        self._cosmo.set(self._params)

    def compute_bins(self, path):
        """
        Compute the w_i bins for the models given in path.
        """
        if self._computed is True:
            print "Bins already computed. Use reset if you want to compute it again"
            return

        self._path = path

        self._read_from_file(path)

        self._create_output_files()

        numbers_of_rows = len(self.params_smg)

        wzbins = np.empty((5, len(self._zbins)))
        wabins = np.empty((5, len(self._abins)))
        params = []
        failed_indexes = []

        for row in range(numbers_of_rows):
            print "{}/{}".format(row+1, numbers_of_rows)
            self._set_params(row)
            index_row = row % 5
            try:
                self._cosmo.compute()
                for n, z in enumerate(self._zbins):
                    wzbins[index_row][n] = self._cosmo.w_smg(z)
                for n, a in enumerate(self._abins):
                    wabins[index_row][n] = self._cosmo.w_smg(1./a-1.)
                params.append(self._params)
            except Exception as e:
                failed_indexes.append(index_row)
                self._cosmo.struct_cleanup()
                print self._params
                print e
                print "\n"
                continue

            self._cosmo.empty()
            self._cosmo.struct_cleanup()

            if row and (not index_row):
                self._save_computed(params,
                                    np.delete(wzbins, failed_indexes, axis=0),
                                    np.delete(wabins, failed_indexes, axis=0))

                params = []
                failed_indexes = []
                wzbins = np.empty((5, len(self._zbins)))
                wabins = np.empty((5, len(self._abins)))

        wzbins[:index_row]  # Remove the unset rows
        wabins[:index_row]
        self._save_computed()
        self._set_empty_wbins()

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
        self._cosmo.empty()
        self._cosmo.struct_cleanup()
        self._set_default_values()
