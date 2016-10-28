#!/usr/bin/python

import numpy as np
import class_headers_to_dict as chd
import matplotlib.pyplot as plt
import variable_title_dics as vtd
from wicmath import relative_deviation
from list_class_headers import check_two_header_lists
import sys


def __filter_lower_x(filename, min_x):
    """In case min_x is provided, filter all values lower than it"""

    with open(filename, 'r') as f:
        for i, line in enumerate(f):
            if not line.startswith('#'):
                if float(line.split()[0]) >= min_x:
                    yield line


def __plot_alone(X, Y):
    f, ax_plt = plt.subplots(1, 1)

    ax_plt.plot(X['data'], Y['data'])
    ax_plt.set_xscale(X['scale'])
    ax_plt.set_yscale(Y['scale'])
    ax_plt.set_xlabel(X['label'])
    ax_plt.set_ylabel(Y['label'])

    return 1


def __plot_compare(X, Y):
    f, (ax_plt, ax_dev) = plt.subplots(2, 1, sharex='col')

    ax_plt.set_xscale(X['scale'])
    ax_plt.set_yscale(Y['scale'])
    ax_plt.set_ylabel(Y['label'])
    ax_plt.plot(X['data'], Y['data'], label=Y['legend'])
    ax_plt.plot(X['cdata'], Y['cdata'], '--', label=Y['clegend'])

    ax_plt.legend(loc='center', bbox_to_anchor=(0.5, -0.06),
          ncol=3, fancybox=True, shadow=True)

    ax_dev.set_xscale(X['scale'])
    ax_dev.set_xlabel(X['label'])
    ax_dev.set_ylabel('rel. dev. [%]')
    x_rel, y_rel = relative_deviation(X['data'], Y['data'], X['cdata'],
                                      Y['cdata'])
    ax_dev.plot(x_rel, y_rel)
    # x_rel, y_rel = relative_deviation(X['data'], Y['data'], X['data'],
    #                                   Y['data'])
    # ax_dev.plot(x_rel, y_rel, '--')

    plt.tight_layout()

    return 1


def __try_loadtxt(filename, min_x, usecols):
    """Try loading data. If variable labels were wrong, print the possible
    options and exit"""
    try:
        x, y = np.loadtxt(__filter_lower_x(filename, min_x), usecols=usecols,
                          unpack=True)

        return x, y

    except TypeError:
        title_var_dic = chd.class_file_title_var_dic(filename)
        title_var_str = ""
        for title, var in title_var_dic.iteritems():
            title_var_str += ''.join("'" + title + "'" + ': ' + "'" + var +
                                     "'" + '\n')

        sys.exit("You probably used wrong variable labels. See bellow the" +
                 " avaibles ones for the present file: " +
                 filename.split('/')[-1] + '\n' + title_var_str)


def plot_utility(filename, x='', y='', min_x=0, x_scale='log', y_scale='log',
                 x_label='', y_label='', y_legend='', compare_with='',
                 compare_with_legend='', output=''):
    """Given a CLASS output file plot its selected variables x, y. If variable
    labels are not set or misspelled, print the possible ones"""
    X = {'var':  x, 'min': min_x, 'scale': x_scale}
    Y = {'var':  y, 'scale': y_scale}

    var_col_dic = chd.class_headers_to_dict(filename)

    usecols = (var_col_dic[x], var_col_dic[y])  # If x/y is not a valid key, it
                                                # returns an empty list. (see
                                                # chd module)

    X['data'], Y['data'] = __try_loadtxt(filename, min_x, usecols)
    X['label'] = x_label or vtd.var_title_dic[x]
    Y['label'] = y_label or vtd.var_title_dic[y]

    Y['legend'] = y_legend or filename.split('/')[-1]  # Label for lengend

    if not compare_with:
        __plot_alone(X, Y)
    else:
        if not check_two_header_lists(filename, compare_with):
            sys.exit("File headers don't coincide")

        X['cdata'], Y['cdata'] = __try_loadtxt(compare_with, min_x, usecols)

        Y['clegend'] = compare_with_legend or compare_with.split('/')[-1]
        __plot_compare(X, Y)

    if output:
        plt.savefig(output)

    else:
        plt.show()

    plt.close()
