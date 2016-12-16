#!/usr/bin/python

import numpy as np
import model_independent_variables as miv
import sys


def __data_gen(data_list):
    for line in data_list:
        yield line


def __try_loadtxt(data, usecols, kind):
    """Try loading data. If variable labels were wrong, print the possible
    options and exit"""
    try:
        x = np.loadtxt(data, usecols=usecols, unpack=True)
        return x

    except TypeError:  # Raised when var_col_dic returns [] (wrong key)
        sys.exit("Check you are using the correct file: needed " + kind)


def __V(theory, IC, phi_smg):
    c = 2.99792458e8  # Speed of light
    M_H = 1.e5 / c  # Hubble rate

    def no_theory():
        sys.exit("Theory not yet implemented.")

    def quintessence_monomial():
        IC['keys'] = ['V0', 'N']
        V0 = IC['V0']
        N = IC['N']
        V = M_H ** 2 * np.multiply(V0, np.power(phi_smg, N))

        return V

    def quintessence_binomial():
        IC['keys'] = ['V1', 'N1', 'V2', 'N2']
        V1 = IC['V1']
        N1 = IC['N1']
        V2 = IC['V2']
        N2 = IC['N2']
        V = M_H ** 2 * np.add(np.multiply(V1, np.power(phi_smg, N1)),
                                np.multiply(V2, np.power(phi_smg, N2)))

        return V

    def quintessence_eft():
        IC['keys'] = ['V0', 'E_F', 'n_min', 'n_Q', 'zeta_2', 'zeta_4']
        V0 = IC['V0']
        E_F = IC['E_F']
        n_min = IC['n_min']
        n_Q = IC['n_Q']
        zeta_2 = IC['zeta_2']
        zeta_4 = IC['zeta_4']
        n_max = n_min + n_Q - 1

        for n in np.arange(n_min, n_max + 1):
            zeta_n = 'zeta_{}'.format(n)
            IC['keys'].append(zeta_n)

        f = zeta_2 * (E_F * phi_smg) ** 2 + zeta_4 * (E_F * phi_smg) ** 4
        fsum = 0
        for n in np.arange(n_min, n_max + 1):
            zeta_n = 'zeta_{}'.format(n)
            fsum += IC[zeta_n] * (E_F * phi_smg) ** n

        V = M_H ** 2 * V0 * (f + fsum)

        return V

    def quintessence_axion():
        IC['keys'] = ['V0', 'E_F', 'E_NP', 'n_max']
        V0 = IC['V0']
        E_F = IC['E_F']
        E_NP = IC['E_NP']
        n_max = IC['n_max']

        for n in np.arange(2, n_max + 1):
            zeta_n = 'zeta_{}'.format(n)
            IC['keys'].append(zeta_n)

        f = 1 + np.cos(E_F * phi_smg)
        fsum = 0
        for n in np.arange(2, n_max + 1):
            zeta_n = 'zeta_{}'.format(n)
            fsum += IC[zeta_n] * E_NP ** (n - 1) * np.cos(n * E_F * phi_smg)

        V = M_H ** 2 * V0 * (f + fsum)

        return V

    def quintessence_modulus():
        IC['keys'] = ['V0', 'E_D', 'p_D', 'n_max', 'alpha']
        V0 = IC['V0']
        E_D = IC['E_D']
        p_D = IC['p_D']
        n_max = IC['n_max']
        alpha = IC['alpha']

        for n in np.arange(0, n_max + 1):
            zeta_n = 'zeta_{}'.format(n)
            IC['keys'].append(zeta_n)

        fsum = 0
        for n in np.arange(0, n_max + 1):
            zeta_n = 'zeta_{}'.format(n)
            fsum += IC[zeta_n] * E_D ** n * np.exp(alpha * (p_D - n) * phi_smg)

        V = M_H ** 2 * V0 * fsum

        return V

    def quintessence_exponential():
        IC['keys'] = ['lambda']
        lbd = IC['lambda']
        V = np.exp(-np.multiply(lbd, phi_smg))

        return V


    theories = {'quintessence_monomial': quintessence_monomial,
                'quintessence_binomial': quintessence_binomial,
                'quintessence_eft': quintessence_eft,
                'quintessence_axion': quintessence_axion,
                'quintessence_modulus': quintessence_modulus,
                'quintessence_exponential': quintessence_exponential}

    try:
        V = theories.get(theory, no_theory)()

    except KeyError:
        sys.exit("IC must be a dictionary with keys: {}".format(IC['keys']))

    return V


def w(theory):
    def w_in(data, var_col_dic, IC):
        """Fill dictionary Y with the equation of state for
        quintessence_monomial. Note data must be a generator or list"""

        usecols = (var_col_dic['z'], var_col_dic["phi_prime_smg"],
                   var_col_dic['phi_smg'])

        z, phi_prime_smg, phi_smg = __try_loadtxt(data, usecols, "background")

        V = __V(theory, IC, phi_smg)

        phi_dot = np.multiply(phi_prime_smg, np.add(z, 1))
        phi_dot2 = np.power(phi_dot, 2)

        w = np.divide(np.subtract(phi_dot2, np.multiply(2, V)),
                      np.add(phi_dot2, np.multiply(2, V)))

        legend = '(phi_dot^2 - 2V)/(phi_dot^2 + 2V)'

        return w, legend
    return w_in


def alphaK(choice, theory):

    # def alphaK1(data, var_col_dic, IC):
    #     """Fill dictionary Y with the evolution of alpha_K for
    #     quintessence_monomial. Note data must be a generator or list."""

    #     data_list = list(data)

    #     usecols = (var_col_dic['rho_b'], var_col_dic['rho_cdm'],
    #             var_col_dic['rho_crit'])

    #     data1 = __data_gen(data_list)

    #     rho_b, rho_cdm, rho_crit = __try_loadtxt(data1, usecols, "background")

    #     Omega_m = np.divide(np.add(rho_b, rho_cdm), rho_crit)

    #     data2 = __data_gen(data_list)

    #     wx, wx_legend = w(theory)(data2, var_col_dic, IC)

    #     alphaK = np.multiply(np.subtract(1, Omega_m), np.add(1, wx))

    #     legend = '(1-Omega_m)(1+w)'

    #     return alphaK, legend

    # def alphaK2(data, var_col_dic, IC):
    #     """Fill dictionary Y with the evolution of alpha_K for
    #     quintessence_monomial. Note data must be a generator or list."""

    #     data_list = list(data)

    #     usecols = (var_col_dic['rho_b'], var_col_dic['rho_cdm'],
    #             var_col_dic['rho_crit'])

    #     data1 = __data_gen(data_list)

    #     rho_b, rho_cdm, rho_crit = __try_loadtxt(data1, usecols, "background")

    #     Omega_m = np.divide(np.add(rho_b, rho_cdm), rho_crit)

    #     data2 = __data_gen(data_list)

    #     x_data, wx, legend = miv.w('z', data2, var_col_dic)

    #     alphaK = np.multiply(np.subtract(1, Omega_m), np.add(1, wx))

    #     legend = '(1-Omega_m)(1+p/rho)'

    #     return alphaK, legend

    def alphaK1(data, var_col_dic, IC):
        """Fill dictionary Y with the evolution of alpha_K for
        quintessence_monomial. Note data must be a generator or list."""

        data_list = list(data)

        usecols = (var_col_dic['rho_smg'], var_col_dic['rho_crit'])

        data1 = __data_gen(data_list)

        rho_smg, rho_crit = __try_loadtxt(data1, usecols, "background")

        Omega_smg = np.divide(rho_smg, rho_crit)

        data2 = __data_gen(data_list)

        wx, wx_legend = w(theory)(data2, var_col_dic, IC)

        alphaK = np.multiply(np.multiply(3, Omega_smg), np.add(1, wx))

        legend = 'Omega_smg(1+w)'

        return alphaK, legend

    def alphaK2(data, var_col_dic, IC):
        """Fill dictionary Y with the evolution of alpha_K for
        quintessence_monomial. Note data must be a generator or list."""

        data_list = list(data)

        usecols = (var_col_dic['rho_smg'], var_col_dic['rho_crit'])

        data1 = __data_gen(data_list)

        rho_smg, rho_crit = __try_loadtxt(data1, usecols, "background")

        Omega_smg = np.divide(rho_smg, rho_crit)

        data2 = __data_gen(data_list)

        x_data, wx, legend = miv.w('z', data2, var_col_dic)

        alphaK = np.multiply(np.multiply(3, Omega_smg), np.add(1, wx))

        legend = 'Omega_smg(1+p/rho)'

        return alphaK, legend

    alpha = {
        1: alphaK1,
        2: alphaK2
        # 3: alphaK3,
        # 4: alphaK4
    }

    return alpha[choice]
