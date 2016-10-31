#!/usr/bin/python

import numpy as np


def w(var_col_dic, data):
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

    y_cdata = np.divide(np.subtract(phi_dot2, np.multiply(2, V)),
                        np.add(phi_dot2, np.multiply(2, V)))

    y_clegend = '(phi_dot^2 - 2V)/(phi_dot^2 + 2V)'

    return y_cdata, y_clegend
