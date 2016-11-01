#!/usr/bin/python

import numpy as np


def w(data, var_col_dic):
    """Fill dictionary Y with the equation of state for
    quintessence_monomial"""

    usecols = (var_col_dic['z'], var_col_dic["phi_prime_smg"],
               var_col_dic['phi_smg'])

    z, phi_prime_smg, phi_smg = np.loadtxt(data, usecols=usecols,
                                           unpack=True)

    ######
    # The following variables are exclusive of monomial quintessence
    ######
    V0 = 1.74646e-07  # This is the value V0* printed in screen when class is
                      # run.
    N = 2
    V = np.multiply(V0, np.power(phi_smg, N))  # Quintessence_monomial
    ######

    phi_dot = np.multiply(phi_prime_smg, np.add(z, 1))
    phi_dot2 = np.power(phi_dot, 2)

    w = np.divide(np.subtract(phi_dot2, np.multiply(2, V)),
                  np.add(phi_dot2, np.multiply(2, V)))

    legend = '(phi_dot^2 - 2V)/(phi_dot^2 + 2V)'

    return w, legend


def w0_wa(data, var_col_dic):
    """Return w0, wa values during the evolution of phi, where w(a) \simeq w0 +
    (a0-a) w'(a0)."""

    usecols = (var_col_dic['z'], var_col_dic["phi_prime_smg"],
               var_col_dic['phi_smg'], var_col_dic['rho_smg'],
               var_col_dic['H'])

    z, phi_prime_smg, phi_smg, rho_smg, H = np.loadtxt(data, usecols=usecols,
                                                       unpack=True)

    ######
    # The following variables are exclusive of monomial quintessence
    ######
    V0 = 1.74646e-07  # This is the value V0* printed in screen when class is
                      # run.
    N = 2
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
