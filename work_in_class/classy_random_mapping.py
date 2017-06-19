#!/usr/bin/python
# coding: utf-8

from abc import ABCMeta  # , abstractmethod
from classy import Class
from classy import CosmoComputationError
from classy import CosmoSevereError
# from classy import CosmoGuessingError
import numpy as np
from datetime import datetime


class Theory(object):

    __metaclass__ = ABCMeta

    def __init__(self):
        self.data = []
        self.data_compu_error = []
        self.data_guess_error = []
        self.name = self.__class__.__name__
        self.common_params_names = ["w0", "wa", "h", "Omega_cdm"]
        # The distribution limits are taken from Marsh et al.
        self.params_dists = {"h": [[0.6, 0.8], "U"],
                             "Omega_cdm": [[0.15, 0.35], "U"]}
        self.params = {
            'Omega_Lambda': 0,
            'Omega_fld': 0,
            'Omega_smg': -1}
            # 'input_verbose': 10}
        self.parameters_smg = []
        self.parameters_2_smg = []
        self.parameters_2_smg_short = []
        self.cosmo = Class()

    def __header(self, kind):
        if kind == "noerror":
            params_names = self.common_params_names + self.model_params_names
        else:
            params_names = self.common_params_names[2:] + self.model_params_names

        header = ""
        for number, name in enumerate(params_names):
            header += "{}:{} ".format(number, name)

        return header.strip() + "\n"  # strip to remove trailin whitespace removed

    def header_noerror(self):
        self.header_noerror = self.__header("noerror")

    def header_error(self):
        self.header_error = self.__header("error")

    def header_intervals(self):
        header = ""
        for param, dist_list in self.params_dists.iteritems():
            if "log" in dist_list[1]:
                header += "log_10({}) \in {}{} ".format(param, dist_list[1][3:], dist_list[0])
            else:
                header += "{} \in {}{} ".format(param, dist_list[1], dist_list[0])

        self.header_intervals = header + "\n"

    def update_headers(self):
        self.header_noerror()
        self.header_error()
        self.header_intervals()

    def update_model_params_names(self, model_params_names):
        self.model_params_names = model_params_names

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

            h, Omega_cdm = self.compute_cosmological_params()
            # TODO: Improve managament of cosmological parameters.

            self.compute_parameters()
            self.params['parameters_smg'] = str(self.parameters_smg).strip('[]')
            if self.parameters_2_smg:
                self.params['parameters_2_smg'] = str(self.parameters_2_smg_short).strip('[]')
            self.params['h'] = h
            self.params['Omega_cdm'] = Omega_cdm
            print [h, Omega_cdm] + self.parameters_smg
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

            #except CosmoGuessingError, e:
            #    print "CosmoGuessing!!!"
            #    print e
            #    if "root must be bracketed in zriddr." in str(e):
            #        self.data_guess_error.append(self.parameters_smg + [h, Omega_cdm]
            #                                    + self.parameters_2_smg)
            #    self.cosmo.struct_cleanup()
            #    self.cosmo.empty()
            #    continue

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
                   self.data, header=self.header_intervals + self.header_noerror)

        if len(self.data_compu_error):
            np.savetxt('./output/' + self.name + '-w0_wa-' + time +
                       '-errors-computation.dat', self.data_compu_error,
                       header=header_compu_error + self.header_intervals +
                       self.header_error)

        if len(self.data_guess_error):
            np.savetxt('./output/' + self.name + '-w0_wa-' + time +
                       '-errors-guessing.dat', self.data_guess_error,
                       header=header_guess_error + self.header_intervals +
                       self.header_error)

    def __compute_params(self, kind):
        params_values = []

        if kind == "common":
            params_names = self.common_params_names
        elif kind == "smg":
            params_names = self.model_params_names

        for param_name in params_names:
            try:
                param_dist_type = self.params_dists[param_name][1]
                param_dist_lims = self.params_dists[param_name][0]
            except:
                continue

            if 'U' in param_dist_type:
                param_value = np.random.uniform(low=param_dist_lims[0], high=param_dist_lims[1])
            elif 'Uz' in param_dist_type:
                param_value = np.random.randint(low=param_dist_lims[0], high=param_dist_lims[1])
            elif 'N' in param_dist_type:
                param_value = np.random.standard_normal()

            if "log" in param_dist_type:
                param_value = 10**param_value

            params_values.append(param_value)

        return params_values

    def compute_parameters_smg(self):
            self.parameters_smg = self.__compute_params("smg")

    def compute_cosmological_params(self):
            return self.__compute_params("common")

    def compute_parameters(self):
        self.compute_parameters_smg()


class quintessence_monomial(Theory):
    def __init__(self):
        Theory.__init__(self)
        self.params['gravity_model'] = 'quintessence_monomial'
        self.model_params_names = ["N", "V0", "phi_prime_ini", "phi_ini"]
        self.params_dists.update({"N": [[1, 7], "U"],
                                  "V0": [[-1, 1], "logU"],
                                  "phi_prime_ini": [[0, 0], "U"],
                                  "phi_ini": [[0, 10], "U"]})

        self.update_headers()

#class quintessence_axion(Theory):
#    def __init__(self):
#        Theory.__init__(self)
#        self.params['gravity_model'] = 'quintessence_axion'
#        self.header_intervals = "h \in U(0.6, 0,8), Omega_cdm \in U(0.15, 0.35), n_max \in Uz(10, 20), log_10(E_F, E_NP) \in U(-3,-1), log_10(V0) \in (-1,1), phi_prime_ini = 1.e-100 and phi_ini \in U[-pi/E_F, pi/E_F), zeta's \in N(0,1) (normal dist).\n "
#        self.header_noerror = "1:w0    2:wa      3:phi_prime_ini    4:phi_ini   5:V0     6:E_F   7:E_NP  8:n_max  9:h   10:Omega_cdm    11:zeta_2    12: zeta_3 ... zeta_(n_max)\n"
#        self.header_error = "1:phi_prime_ini    2:phi_ini   3:V0     4:E_F   5:E_NP     6:n_max  7:h  8:Omega_cdm   9:zeta_2    10: zeta_3 ... zeta_(n_max)\n"
#
#    def compute_parameters(self):
#        # phi_prime_ini = 10 ** np.random.uniform(low=0, high=2)
#        V0 = 10 ** np.random.uniform(low=-1, high=1)  # It will be overwritten by tunning_index unless
#        E_F = 10 ** np.random.uniform(low=-3, high=-1)
#        E_NP = 10 ** np.random.uniform(low=-3, high=-1)
#        n_max = np.random.randint(low=10, high=21)  # The highest number is 20.
#        phi_prime_ini = 1.e-100  # In Marsh et al. phi'=0, but hi_class find nan with it.
#        phi_ini = np.random.uniform(low=-np.pi/E_F, high=np.pi/E_F)
#
#        self.parameters_smg = [phi_prime_ini, phi_ini, V0, E_F, E_NP, n_max]
#        self.parameters_2_smg = [np.nan] * 19  # np.savetxt need all arrays of same size.
#        self.parameters_2_smg_short = list(np.random.standard_normal(n_max - 1))
#        for index, item in enumerate(self.parameters_2_smg_short):
#            self.parameters_2_smg[index] = item
#
#        return [self.parameters_smg, self.parameters_2_smg]
#
#
#class quintessence_eft(Theory):
#    def __init__(self):
#        Theory.__init__(self)
#        self.params['gravity_model'] = 'quintessence_eft'
#        self.header_intervals = "h \in U(0.6, 0,8), Omega_cdm \in U(0.15, 0.35),n_min & n_Q \in Uz(5,10), log_10(E_F) \in U(-3,-1), log_10(V0) \in (-1,1), phi_prime_ini = 1.e-100 and phi_ini \in U[-1/E_F, 1/E_F), zeta's \in N(0,1) (normal dist).\n"
#        self.header_noerror = "1:w0    2:wa      3:phi_prime_ini    4:phi_ini   5:V0     6:E_F   7:n_min     8:n_Q   9:h     10: Omega_cdm      11:zeta_2    12: zeta_4 ... zeta_(n_Q + 2)\n"
#        self.header_error = "1:phi_prime_ini    2:phi_ini   3:V0     4:E_F   5:n_min     6:n_Q   7:h   8:Omega_cdm      9:zeta_2    10: zeta_4 ... zeta_(n_Q + 2)\n"
#
#    def compute_parameters(self):
#        V0 = 10 ** np.random.uniform(low=-1, high=1)  # It will be overwritten by tunning_index unless
#        E_F = 10 ** np.random.uniform(low=-3, high=-1)
#        n_min = np.random.randint(low=5, high=11)  # The highest number is 10.
#        n_Q = np.random.randint(low=5, high=11)  # The highest number is 10.
#        # phi_prime_ini = 10 ** np.random.uniform(low=0, high=2)
#        phi_prime_ini = 1.e-100  # In Marsh et al. phi'=0, but hi_class find nan with it.
#        phi_ini = np.random.uniform(low=-1/E_F, high=1/E_F)
#
#        self.parameters_smg = [phi_prime_ini, phi_ini, V0, E_F, n_min, n_Q]
#        self.parameters_2_smg = [np.nan] * 12  # np.savetxt need all arrays of same size.
#        self.parameters_2_smg_short = list(np.random.standard_normal(2 + n_Q))
#        for index, item in enumerate(self.parameters_2_smg_short):
#            self.parameters_2_smg[index] = item
#
#        return [self.parameters_smg, self.parameters_2_smg]
#
#
#class quintessence_modulus(Theory):
#    def __init__(self):
#        Theory.__init__(self)
#        self.params['gravity_model'] = 'quintessence_modulus'
#        self.header_intervals = "h \in U(0.6, 0,8), Omega_cdm \in U(0.15, 0.35), n_max \in Uz(10,20), p_D \in Uz(1,5), alpha \in U(0, 1), log_10(E_D) \in U(-3,-1), log_10(V0) \in (-1,1), phi_prime_ini = 1.e-100 and phi_ini \in U[-1,1), zeta's \in N(0,1) (normal dist).\n"
#        self.header_noerror = "1:w0    2:wa      3:phi_prime_ini    4:phi_ini   5:V0     6:E_D   7:p_D   8:n_max     9:alpha     10:h    11:Omega_cdm   12:zeta_0    13: zeta_1 ... zeta_(n_max)\n"
#        self.header_error = "1:phi_prime_ini    2:phi_ini   3:V0     4:E_D   5:p_D   6:n_max     7:alpha   8:h    9:Omega_cdm   10:zeta_0    11: zeta_1 ... zeta_(n_max)\n"
#
#    def compute_parameters(self):
#        # phi_prime_ini = 10 ** np.random.uniform(low=0, high=2)
#        phi_prime_ini = 1.e-100  # In Marsh et al. phi'=0, but hi_class find nan with it.
#        phi_ini = np.random.uniform(low=-1, high=1)
#        V0 = 10 ** np.random.uniform(low=-1, high=1)  # It will be overwritten by tunning_index unless
#        E_D = 10 ** np.random.uniform(low=-3, high=-1)
#        p_D = np.random.randint(low=1, high=6)  # The highest number is 5.
#        n_max = np.random.randint(low=10, high=21)  # The highest number is 20.
#        alpha = np.random.uniform(low=0, high=1)
#
#        self.parameters_smg = [phi_prime_ini, phi_ini, V0, E_D, p_D, n_max, alpha]
#        self.parameters_2_smg = [np.nan] * 21  # np.savetxt need all arrays of same size.
#        self.parameters_2_smg_short = list(np.random.standard_normal(1 + n_max))
#        for index, item in enumerate(self.parameters_2_smg_short):
#            self.parameters_2_smg[index] = item
#
#        return [self.parameters_smg, self.parameters_2_smg]
#


class quintessence_binomial(Theory):
    def __init__(self):
        Theory.__init__(self)
        self.params['gravity_model'] = 'quintessence_binomial'
        self.header_intervals = "h \in U(0.6, 0,8), Omega_cdm \in U(0.15, 0.35), N_i \in Uz(1,7), log_10(V_i) \in (-1,1), phi_prime_ini = 1.e-100 and phi_ini \in U[0,4).\n"
        self.header_noerror = "1:w0    2:wa      3:phi_prime_ini    4:phi_ini    5:N1    6:V1    7:N2    8:V2    9:h    10:Omega_cdm\n"
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
        self.model_params_names = ["phi_prime_ini", "f_ini", "alpha", "c", "p", "n"]
        self.params_dists.update({"phi_prime_ini": [[0, 0], "U"],
                                      "f_ini": [[0, 10], "U"],
                                      "alpha": [[0, 10], "U"],
                                      "c": [[1, 1], "U"],
                                      "p": [[1, 10], "U"],
                                      "n": [[1, 10], "U"]})
        self.update_headers()


class galileon(Theory):
    def __init__(self, n):
        Theory.__init__(self)
        self.params['gravity_model'] = 'galileon'
        self.model_params_names = ["xi", "c_1", "c_2", "c_3", "c_4", "c_5", "phi_ini"]
        self.smg_params_dists.update({"xi": [[-10, 10], "U"],
                                      "c_1": [[0, 0], "U"],
                                      "c_2": [[-1, -1], "U"],
                                      "c_3": [[1e-3, 1e-3], "U"],
                                      "c_4": [[1e-3, 1e-3], "U"],
                                      "c_5": [[1e-3, 1e-3], "U"],
                                      "phi_ini": [[0, 0], "U"]
                                      })

        self.smg_params_dists["c_{}".format(n)][0] = [-5.5, 5]


class prado_nelson(Theory):
    def __init__(self):
        Theory.__init__(self)
        self.params['gravity_model'] = 'prado_nelson'
        self.model_params_names = ["phi_prime_ini", "phi_ini", "alpha", "beta",
                                   "m", "n", "c", "lambda"]
        self.params_dists.update({"phi_prime_ini": [[-20, -15], "logU"],
                                      "phi_ini": [[0, 0], "U"],
                                      "alpha": [[-10, 1], "logU"],
                                      "beta": [[-1, 1], "U"],
                                      "m": [[0, 5], "U"],
                                      "n": [[0, 5], "U"],
                                      "c": [[-5, 5], "U"],
                                      "lambda": [[0.7, 0.7], "U"]})

        self.update_headers()


def main():
    parser = argparse.ArgumentParser(
        description="Explore the parameter space of different quintessence\
        models and obtain wa-w0 data as in Marsh et al.")

    parser.add_argument("theory", help="Quintessence model to explore")
    parser.add_argument("-p", "--points", type=int, help="Number of points to compute")

    args = parser.parse_args()
    d = {'quintessence_monomial': quintessence_monomial,
         'quintessence_binomial': quintessence_binomial,
         # 'quintessence_axion': quintessence_axion,
         # 'quintessence_eft': quintessence_eft,
         # 'quintessence_modulus': quintessence_modulus,
         'alpha_attractor_canonical': alpha_attractor_canonical,
         'gal3': galileon,
         'gal4': galileon,
         'gal5': galileon,
         'prado_nelson': prado_nelson}

    if 'gal' in args.theory:
        th = d[args.theory](int(args.theory[-1]))
    else:
        th = d[args.theory]()
    th.compute_data(args.points)
    th.savetxt()

if __name__ == "__main__":
    import argparse
    main()
