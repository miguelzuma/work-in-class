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


def w(theory):
    def __V(theory, IC, phi_smg):
        try:
            if theory == 'quintessence_monomial':
                IC_keys = ['V0', 'N']
                # This is the value V0* printed in screen when class is run.
                V0 = IC['V0']
                N = IC['N']
                V = np.multiply(V0, np.power(phi_smg, N))

            elif theory == 'quintessence_exponential':
                IC_keys = ['lambda']
                lbd = IC['lambda']
                V = np.exp(-np.multiply(lbd, phi_smg))

            else:
                sys.exit("Theory not yet implemented.")

        except KeyError:
            sys.exit("IC must be a dictionary with keys: " + IC_keys)

        return V

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


def w0_wa(data, var_col_dic, IC):
    """Return w0, wa values during the evolution of phi, where w(a) \simeq w0 +
    (a0-a) w'(a0). Note data must be a generator or list."""

    usecols = (var_col_dic['z'], var_col_dic["phi_prime_smg"],
               var_col_dic['phi_smg'], var_col_dic['rho_smg'],
               var_col_dic['H'])

    z, phi_prime_smg, phi_smg, rho_smg, H = __try_loadtxt(data, usecols,
                                                       "background")

    ######
    # The following variables are exclusive of monomial quintessence
    ######
    IC_keys = ['V0', 'N']
    try:
        V0 = IC['V0']  # This is the value V0* printed in screen when class is
                       # run.
        N = IC['N']
    except KeyError:
        sys.exit("IC must be a dictionary with keys: " + IC_keys)

    P = np.power(phi_smg, N)
    P_phi = np.multiply(N, np.power(phi_smg, N - 1))
    V = np.multiply(V0, P)  # Quintessence_monomial
    ######

    phi_dot = np.multiply(phi_prime_smg, np.add(z, 1))
    phi_dot2 = np.power(phi_dot, 2)

    w0 = np.divide(np.subtract(phi_dot2, np.multiply(2, V)),
                   np.add(phi_dot2, np.multiply(2, V)))

    pre_factor = np.divide(np.multiply(np.add(1, z), np.divide(2 * V0, H)),
                           np.power(rho_smg, 2))

    summand1 = np.multiply(np.multiply(np.multiply(3, P), phi_dot2), H)

    summand2 = np.multiply(np.multiply(P_phi, phi_dot), rho_smg)

    wa = np.multiply(pre_factor, np.add(summand1, summand2))

    return w0, wa


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
