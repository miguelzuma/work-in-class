#!/usr/bin/python
# coding: utf-8

from abc import ABCMeta, abstractmethod
from classy import Class
from classy import CosmoComputationError
from classy import CosmoSevereError
from classy import CosmoGuessingError
import numpy as np
from datetime import datetime


class Theory(object):

    __metaclass__ = ABCMeta

    def __init__(self):
        self.data = []
        self.data_compu_error = []
        self.data_guess_error = []
        self.name = self.__class__.__name__
        self.params = {
            'Omega_Lambda': 0,
            'Omega_fld': 0,
            'Omega_smg': -1}
        self.parameters_smg = []
        self.parameters_2_smg = []

    def compute_data(self, points):
        cosmo = Class()
        while len(self.data) < points:
            print "######## Point: {} of {}".format(len(self.data) + 1, points)

            # The distribution limits are taken from Marsh et al.
            Omega_cdm = np.random.uniform(low=0.15, high=0.35)
            h = np.random.uniform(low=0.6, high=0.8)
            #

            self.compute_parameters()
            self.params['parameters_smg'] = str(self.parameters_smg).strip('[]')
            if self.parameters_2_smg:
                self.params['parameters_2_smg'] = str(self.parameters_2_smg).strip('[]')
            self.params['h'] = h
            self.params['Omega_cdm'] = Omega_cdm
            print self.parameters_smg + [h, Omega_cdm]

            cosmo.set(self.params)
            try:
                cosmo.compute()

            except CosmoSevereError:
                print 'CosmoSevere!!!!'
                break

            except CosmoComputationError:
                print 'CosmoCompu!!!'
                self.data_compu_error.append(self.parameters_smg + [h, Omega_cdm])
                cosmo.struct_cleanup()
                cosmo.empty()
                continue

            except CosmoGuessingError:
                print 'CosmoGuessing!!!'
                self.data_guess_error.append(self.parameters_smg + [h, Omega_cdm])
                cosmo.struct_cleanup()
                cosmo.empty()
                continue

            self.data.append([cosmo.w0_smg(), cosmo.wa_smg()] + self.parameters_smg + [h, Omega_cdm])
            cosmo.struct_cleanup()
            cosmo.empty()

    def savetxt(self):
        time = datetime.now().strftime("%Y%m%d%H%M")
        header_compu_error = "Init values for which computation fails.\n"
        header_guess_error = "Init values for which guessing fails.\n"

        np.savetxt('./output/' + self.name + '-w0_wa-' + time + '.dat',
                   self.data, header=self.header_common + self.header)

        if len(self.data_compu_error):
            np.savetxt('./output/' + self.name + '-w0_wa-' + time +
                       '-errors-computation.dat', self.data_compu_error,
                       header=header_compu_error + self.header_common +
                       self.header_error)

        if len(self.data_guess_error):
            np.savetxt('./output/' + self.name + '-w0_wa-' + time +
                       '-errors-guessing.dat', self.data_guess_error,
                       header=header_guess_error + self.header_common +
                       self.header_error)

    @abstractmethod
    def compute_parameters(self):
        pass


class quintessence_monomial(Theory):
    def __init__(self):
        Theory.__init__(self)
        self.params['gravity_model'] = 'quintessence_monomial'
        self.header_common = "N \in U(1,7), log_10(V0) \in (-1,1), log_10(phi_prime_ini) \in U[0,2) and log_10(phi_ini) \in U[0,1.7), h \in U(0.6, 0.8)\n"
        self.header = "1:w0    2:wa    3:N    4:V0     5:phi_prime_ini    6:phi_ini  7:h    8:Omega_cdm"
        self.header_error = "1:N    2:V0    3:phi_prime_ini    4:phi_ini  5:h   6:Omega_cdm"

    def compute_parameters(self):
        N = np.random.uniform(low=1, high=7)
        V0 = 10 ** np.random.uniform(low=-1, high=1)  # It will be overwritten by tunning_index unless
        phi_prime_ini = 10 ** np.random.uniform(low=0, high=2)
        phi_ini = 10 ** np.random.uniform(low=0, high=1.7)  # 1.7 max phi_ini=50

        self.parameters_smg = [V0, N, phi_prime_ini, phi_ini]

        return self.parameters_smg


def main():
    parser = argparse.ArgumentParser(
        description="Explore the parameter space of different quintessence\
        models and obtain wa-w0 data as in Marsh et al.")

    parser.add_argument("theory", help="Quintessence model to explore")
    parser.add_argument("-p", "--points", type=int, help="Number of points to compute")

    args = parser.parse_args()
    d = {'quintessence_monomial': quintessence_monomial}

    th = d[args.theory]()
    th.compute_data(args.points)
    th.savetxt()

if __name__ == "__main__":
    import argparse
    main()
