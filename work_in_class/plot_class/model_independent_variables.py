#!/usr/bin/python

import numpy as np
import sys


def __try_loadtxt(data, usecols, kind):
    """Try loading data. If variable labels were wrong, print the possible
    options and exit"""
    try:
        x = np.loadtxt(data, usecols=usecols, unpack=True)
        return x

    except TypeError:  # Raised when var_col_dic returns [] (wrong key)
        sys.exit("Check you are using the correct file: needed " + kind)


def columns(x, y):
    def columns_data(X, data, var_col_dic):
        X = x  # X kept for compatibility. Needed three arguments

        usecols = (var_col_dic[X], var_col_dic[y])

        x_data, y_data = __try_loadtxt(data, usecols, 'any.')

        legend = ''

        return x_data, y_data, legend

    return columns_data


def alphaK(kineticity_safe):
    def alphaK_corrected(x, data, var_col_dic):

        usecols = (var_col_dic[x], var_col_dic['kineticity_smg'])

        x_data, alpha_K = __try_loadtxt(data, usecols, 'background')

        alpha_K = np.subtract(alpha_K, kineticity_safe)

        legend = 'alpha_k'

        return x_data, alpha_K, legend

    return alphaK_corrected


def w(x, data, var_col_dic):

    usecols = (var_col_dic[x], var_col_dic['p_smg'], var_col_dic['rho_smg'])

    x_data, p_smg, rho_smg = __try_loadtxt(data, usecols, 'background')
    w = np.divide(p_smg, rho_smg)

    legend = 'p/rho'

    return x_data, w, legend
