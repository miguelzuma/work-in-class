#!/usr/bin/python

import numpy as np
import os
from classy import Class


class Binning():
    def __init__(self, zbins=[], abins=[]):
        self._cosmo = Class()
        self._zbins = zbins
        self._abins = abins
        self._set_default_values()

    def _set_default_values(self):
        """
        Set default values of parameters lists.
        """
        self.params_smg = []
        self.params_smg_2 = []
        self.h = []
        self.Omega_cdm = []

        self.gravity_theory = []

        self._params = {"Omega_Lambda": 0,
                        "Omega_fld": 0,
                        "Omega_smg": -1}

        self.wzbins = np.ones(len(self._zbins))
        self.wabins = np.ones(len(self._abins))

        self._computed = False
        self._path = []

    def set_bins(self, zbins, abins):
        """
        Set what bins to use
        """
        self._zbins = zbins
        self._abins = abins

    def _read_from_file(self, path):
        """
        Return params for class from files used in quintessence Marsh.
        """
        with open(path) as f:
            f.readline()
            header = f.readline()[3:].split()  # Remove '#', 'w0', and 'wa'

        columns = np.loadtxt(path, unpack=True)[2:]  # Remove columns w0, wa

        for index_h, head in enumerate(header):
            if head[-1] == 'h':
                break

        self.params_smg = zip(*columns[:index_h])
        self.params_smg_2 = zip(*columns[index_h+2:])
        self.h = columns[index_h]
        self.Omega_cdm = columns[index_h+1]

        self.gravity_theory = os.path.basename(path).split('-')[0]

    def _set_params(self, row):
        """
        Set parameters.
        """
        self._params.update({
            'parameters_smg': str(self.params_smg[row]).strip('()'),
            'parameter_smg_2': str(self.params_smg_2[row]).strip('()'),
            'h': self.h[row],
            'Omega_cdm': self.Omega_cdm[row],
            'gravity_theory': self.gravity_theory
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

        for row in len(self.params_smg):
            self._set_params(row)
            self._cosmo.compute()
            for z in self._zbins:
                self.wzbins = self._cosmo.w(z)
            for a in self._zbins:
                self.wabins = self._cosmo.w(1./a-1.)

        self._cosmo.empty()
        self._cosmo.struct_cleanup()

    def savetxt(self, fname, extension):
        """
        Save the w bins computed to two different files
        'fname'-wz-bins.'extension' and 'fname'-wa-bins.'extension'.
        """
        if self.computed is True:
            print "Bins not yet computed. Compute them first"
            return

        headerZ = '\n'.join(["w(z_i) bins correspondig to models in file " +
                             self._path, str(self.wzbins)])

        np.savetxt(fname+'wz-bins'+extension, self.wzbins, header=headerZ)

        headerA = '\n'.join(["w(a_i) bins correspondig to models in file " +
                             self._path, str(self.wabins)])

        np.savetxt(fname+'wa-bins'+extension, self.wabins, header=headerA)

    def reset(self):
        """
        Reset class
        """
        self._cosmo.empty()
        self._cosmo.struct_cleanup()
        self._set_default_values()
