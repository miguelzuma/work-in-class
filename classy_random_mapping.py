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
            # 'input_verbose': 10}
        self.parameters_smg = []
        self.parameters_2_smg = []
        self.parameters_2_smg_short = []
        self.cosmo = Class()

    def compute_model(self, parameters, debug=False):

        self.parameters_smg = parameters[:len(self.parameters_smg)]
        self.params['parameters_smg'] = str(self.parameters_smg).strip('[]')
        if self.parameters_2_smg:
            self.parameters_2_smg = parameters[len(self.parameters_smg):-2]
            self.params['parameters_2_smg'] = str(self.parameters_2_smg_short).strip('[]')
        self.params['h'] = parameters[-2]
        self.params['Omega_cdm'] = parameters[-1]

        if debug:
            del self.params['Omega_smg']
            self.params['Omega_smg_debug'] = -1

        self.cosmo.set(self.params)

        if debug:
            self.params['Omega_smg'] = -1

        self.cosmo.compute()


    def model_clean(self):
        self.cosmo.struct_cleanup()
        self.cosmo.empty()

    def compute_data(self, points):
        while len(self.data) < points:
            print "######## Point: {} of {}".format(len(self.data) + 1, points)

            # The distribution limits are taken from Marsh et al.
            Omega_cdm = np.random.uniform(low=0.15, high=0.35)
            h = np.random.uniform(low=0.6, high=0.8)
            #

            self.compute_parameters()
            self.params['parameters_smg'] = str(self.parameters_smg).strip('[]')
            if self.parameters_2_smg:
                self.params['parameters_2_smg'] = str(self.parameters_2_smg_short).strip('[]')
            self.params['h'] = h
            self.params['Omega_cdm'] = Omega_cdm
            print self.parameters_smg + [h, Omega_cdm]
            print self.parameters_2_smg

            self.cosmo.set(self.params)
            try:
                self.cosmo.compute()

            except CosmoSevereError, e:
                print "CosmoSevere!!!"
                print e
                break

            except CosmoComputationError, e:
                print "CosmoCompu!!!"
                print e
                self.data_compu_error.append(self.parameters_smg + [h, Omega_cdm]
                                             + self.parameters_2_smg)
                self.cosmo.struct_cleanup()
                self.cosmo.empty()
                continue

            except CosmoGuessingError, e:
                print "CosmoGuessing!!!"
                print e
                if "root must be bracketed in zriddr." in str(e):
                    self.data_guess_error.append(self.parameters_smg + [h, Omega_cdm]
                                                + self.parameters_2_smg)
                self.cosmo.struct_cleanup()
                self.cosmo.empty()
                continue

            self.data.append([self.cosmo.w0_smg(), self.cosmo.wa_smg()] +
                             self.parameters_smg + [h, Omega_cdm] +
                             self.parameters_2_smg)

            self.cosmo.struct_cleanup()
            self.cosmo.empty()

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
        #self.header_common = "h \in U(0.6, 0.8), Omega_cdm \in U(0.15, 0.35), N \in Uz(1,7), log_10(V0) \in (-1,1), phi_prime_ini = 1.e-100 and phi_ini \in U[0,4),\n"
        self.header_common = "h \in U(0.6, 0.8), Omega_cdm \in U(0.15, 0.35), N \in U(1,7), log_10(V0) \in (-1,1), phi_prime_ini = 1.e-100 and phi_ini \in U[0,10),\n"
        self.header = "1:w0    2:wa    3:N    4:V0     5:phi_prime_ini    6:phi_ini  7:h    8:Omega_cdm"
        self.header_error = "1:N    2:V0    3:phi_prime_ini    4:phi_ini  5:h   6:Omega_cdm"

    def compute_parameters(self):
        # N = np.random.uniform(low=1, high=7)
        # phi_prime_ini = 10 ** np.random.uniform(low=0, high=2)
        #N = np.random.randint(low=1, high=7)
        N = np.random.uniform(low=1, high=7) #Not as in Marsh (randint)
        V0 = 10 ** np.random.uniform(low=-1, high=1)  # It will be overwritten by tunning_index unless
        phi_prime_ini = 1.e-100  # In Marsh et al. phi'=0, but hi_class find nan with it.
        phi_ini = np.random.uniform(low=0, high=10) #Not as in Marsh (high=4)

        self.parameters_smg = [N, V0, phi_prime_ini, phi_ini]

        return [self.parameters_smg, self.parameters_2_smg]


class quintessence_axion(Theory):
    def __init__(self):
        Theory.__init__(self)
        self.params['gravity_model'] = 'quintessence_axion'
        self.header_common = "h \in U(0.6, 0,8), Omega_cdm \in U(0.15, 0.35), n_max \in Uz(10, 20), log_10(E_F, E_NP) \in U(-3,-1), log_10(V0) \in (-1,1), phi_prime_ini = 1.e-100 and phi_ini \in U[-pi/E_F, pi/E_F), zeta's \in N(0,1) (normal dist).\n "
        self.header = "1:w0    2:wa      3:phi_prime_ini    4:phi_ini   5:V0     6:E_F   7:E_NP  8:n_max  9:h   10:Omega_cdm    11:zeta_2    12: zeta_3 ... zeta_(n_max)\n"
        self.header_error = "1:phi_prime_ini    2:phi_ini   3:V0     4:E_F   5:E_NP     6:n_max  7:h  8:Omega_cdm   9:zeta_2    10: zeta_3 ... zeta_(n_max)\n"

    def compute_parameters(self):
        # phi_prime_ini = 10 ** np.random.uniform(low=0, high=2)
        V0 = 10 ** np.random.uniform(low=-1, high=1)  # It will be overwritten by tunning_index unless
        E_F = 10 ** np.random.uniform(low=-3, high=-1)
        E_NP = 10 ** np.random.uniform(low=-3, high=-1)
        n_max = np.random.randint(low=10, high=21)  # The highest number is 20.
        phi_prime_ini = 1.e-100  # In Marsh et al. phi'=0, but hi_class find nan with it.
        phi_ini = np.random.uniform(low=-np.pi/E_F, high=np.pi/E_F)

        self.parameters_smg = [phi_prime_ini, phi_ini, V0, E_F, E_NP, n_max]
        self.parameters_2_smg = [np.nan] * 19  # np.savetxt need all arrays of same size.
        self.parameters_2_smg_short = list(np.random.standard_normal(n_max - 1))
        for index, item in enumerate(self.parameters_2_smg_short):
            self.parameters_2_smg[index] = item

        return [self.parameters_smg, self.parameters_2_smg]


class quintessence_eft(Theory):
    def __init__(self):
        Theory.__init__(self)
        self.params['gravity_model'] = 'quintessence_eft'
        self.header_common = "h \in U(0.6, 0,8), Omega_cdm \in U(0.15, 0.35),n_min & n_Q \in Uz(5,10), log_10(E_F) \in U(-3,-1), log_10(V0) \in (-1,1), phi_prime_ini = 1.e-100 and phi_ini \in U[-1/E_F, 1/E_F), zeta's \in N(0,1) (normal dist).\n"
        self.header = "1:w0    2:wa      3:phi_prime_ini    4:phi_ini   5:V0     6:E_F   7:n_min     8:n_Q   9:h     10: Omega_cdm      11:zeta_2    12: zeta_4 ... zeta_(n_Q + 2)\n"
        self.header_error = "1:phi_prime_ini    2:phi_ini   3:V0     4:E_F   5:n_min     6:n_Q   7:h   8:Omega_cdm      9:zeta_2    10: zeta_4 ... zeta_(n_Q + 2)\n"

    def compute_parameters(self):
        V0 = 10 ** np.random.uniform(low=-1, high=1)  # It will be overwritten by tunning_index unless
        E_F = 10 ** np.random.uniform(low=-3, high=-1)
        n_min = np.random.randint(low=5, high=11)  # The highest number is 10.
        n_Q = np.random.randint(low=5, high=11)  # The highest number is 10.
        # phi_prime_ini = 10 ** np.random.uniform(low=0, high=2)
        phi_prime_ini = 1.e-100  # In Marsh et al. phi'=0, but hi_class find nan with it.
        phi_ini = np.random.uniform(low=-1/E_F, high=1/E_F)

        self.parameters_smg = [phi_prime_ini, phi_ini, V0, E_F, n_min, n_Q]
        self.parameters_2_smg = [np.nan] * 12  # np.savetxt need all arrays of same size.
        self.parameters_2_smg_short = list(np.random.standard_normal(2 + n_Q))
        for index, item in enumerate(self.parameters_2_smg_short):
            self.parameters_2_smg[index] = item

        return [self.parameters_smg, self.parameters_2_smg]


class quintessence_modulus(Theory):
    def __init__(self):
        Theory.__init__(self)
        self.params['gravity_model'] = 'quintessence_modulus'
        self.header_common = "h \in U(0.6, 0,8), Omega_cdm \in U(0.15, 0.35), n_max \in Uz(10,20), p_D \in Uz(1,5), alpha \in U(0, 1), log_10(E_D) \in U(-3,-1), log_10(V0) \in (-1,1), phi_prime_ini = 1.e-100 and phi_ini \in U[-1,1), zeta's \in N(0,1) (normal dist).\n"
        self.header = "1:w0    2:wa      3:phi_prime_ini    4:phi_ini   5:V0     6:E_D   7:p_D   8:n_max     9:alpha     10:h    11:Omega_cdm   12:zeta_0    13: zeta_1 ... zeta_(n_max)\n"
        self.header_error = "1:phi_prime_ini    2:phi_ini   3:V0     4:E_D   5:p_D   6:n_max     7:alpha   8:h    9:Omega_cdm   10:zeta_0    11: zeta_1 ... zeta_(n_max)\n"

    def compute_parameters(self):
        # phi_prime_ini = 10 ** np.random.uniform(low=0, high=2)
        phi_prime_ini = 1.e-100  # In Marsh et al. phi'=0, but hi_class find nan with it.
        phi_ini = np.random.uniform(low=-1, high=1)
        V0 = 10 ** np.random.uniform(low=-1, high=1)  # It will be overwritten by tunning_index unless
        E_D = 10 ** np.random.uniform(low=-3, high=-1)
        p_D = np.random.randint(low=1, high=6)  # The highest number is 5.
        n_max = np.random.randint(low=10, high=21)  # The highest number is 20.
        alpha = np.random.uniform(low=0, high=1)

        self.parameters_smg = [phi_prime_ini, phi_ini, V0, E_D, p_D, n_max, alpha]
        self.parameters_2_smg = [np.nan] * 21  # np.savetxt need all arrays of same size.
        self.parameters_2_smg_short = list(np.random.standard_normal(1 + n_max))
        for index, item in enumerate(self.parameters_2_smg_short):
            self.parameters_2_smg[index] = item

        return [self.parameters_smg, self.parameters_2_smg]


class quintessence_binomial(Theory):
    def __init__(self):
        Theory.__init__(self)
        self.params['gravity_model'] = 'quintessence_binomial'
        self.header_common = "h \in U(0.6, 0,8), Omega_cdm \in U(0.15, 0.35), N_i \in Uz(1,7), log_10(V_i) \in (-1,1), phi_prime_ini = 1.e-100 and phi_ini \in U[0,4).\n"
        self.header = "1:w0    2:wa      3:phi_prime_ini    4:phi_ini    5:N1    6:V1    7:N2    8:V2    9:h    10:Omega_cdm\n"
        self.header_error = "1:phi_prime_ini    2:phi_ini     3:N1    4:V1    5:N2    6:V2    7:h     8:Omega_cdm\n"

    def compute_parameters(self):
        # phi_prime_ini = 10 ** np.random.uniform(low=0, high=2)
        phi_prime_ini = 1.e-100 # In Marsh et al. phi'=0, but hi_class find nan with it.
        phi_ini = np.random.uniform(low=0, high=4)
        N1, N2 = np.random.uniform(low=1, high=7, size=2)
        V1, V2 = 10 ** np.random.uniform(low=-1, high=1, size=2)  # It will be overwritten by tunning_index unless
        # Sampleamos log_10(phi_prime_ini).

        self.parameters_smg = [phi_prime_ini, phi_ini, N1, V1, N2, V2]

        return [self.parameters_smg, self.parameters_2_smg]


class alpha_attractor_canonical(Theory):
    def __init__(self):
        Theory.__init__(self)
        self.params['gravity_model'] = 'alpha_attractor_canonical'
        self.header_common = "h \in U(0.6, 0,8), Omega_cdm \in U(0.15, 0.35), phi_prime_ini = 1.e-100, phi_ini \in U[0,10), log_10(alpha) \in U(-2, 2), log_10(c)\in U(-2,2), p \in U(0,10), n \in U(0,10)' .\n"
        self.header = "1:w0    2:wa      3:phi_prime_ini    4:phi_ini    5:alpha    6:c    7:p    8:n    9:h    10:Omega_cdm\n"
        self.header_error = "1:phi_prime_ini    2:phi_ini     3:alpha    4:c    5:p    6:n    7:h     8:Omega_cdm\n"

    def compute_parameters(self):
        phi_prime_ini = 1.e-100
        phi_ini = np.random.uniform(low=0, high=10)
        alpha = 10 ** np.random.uniform(low=-3, high=3)
        c = 1 #10 ** np.random.uniform(low=-2, high=2)  # It will be overwritten by tunning_index unless
        p = np.random.uniform(0, 10)
        n = np.random.uniform(0, 10)

        self.parameters_smg = [phi_prime_ini, phi_ini, alpha, c, p, n]

        return [self.parameters_smg, self.parameters_2_smg]


def main():
    parser = argparse.ArgumentParser(
        description="Explore the parameter space of different quintessence\
        models and obtain wa-w0 data as in Marsh et al.")

    parser.add_argument("theory", help="Quintessence model to explore")
    parser.add_argument("-p", "--points", type=int, help="Number of points to compute")

    args = parser.parse_args()
    d = {'quintessence_monomial': quintessence_monomial,
         'quintessence_binomial': quintessence_binomial,
         'quintessence_axion': quintessence_axion,
         'quintessence_eft': quintessence_eft,
         'quintessence_modulus': quintessence_modulus,
         'alpha_attractor_canonical': alpha_attractor_canonical}

    th = d[args.theory]()
    th.compute_data(args.points)
    th.savetxt()

if __name__ == "__main__":
    import argparse
    main()
