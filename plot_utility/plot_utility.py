#!/usr/bin/python

import class_headers_to_dict as chd
import matplotlib.pyplot as plt
import variable_title_dics as vtd
import wicmath
# from list_class_headers import check_two_header_lists
import quintessence
import numpy as np
import matplotlib.cm as cm
import sys


def __colors(length):
    color = iter(cm.rainbow(np.linspace(0, 1, length)))

    return color


def __plot_close(output):
    if output:
        plt.savefig(output)

    else:
        plt.show()

    plt.close()


def __filter_lower_x(filename, x_col, min_x):
    """Filter values lower than min_x in x_col"""

    with open(filename, 'r') as f:
        for i, line in enumerate(f):
            if not line.startswith('#'):
                if float(line.split()[x_col]) >= min_x:
                    yield line


def __filter_data(filename, col, min_x):
    """In case min_x is provided, filter all values lower than it"""

    if not min:
        return filename

    else:
        return __filter_lower_x(filename, col, min_x)


def __deviation(x, y, cx, cy, kind):

    if 'abs' in kind:
        x_rel, y_rel = wicmath.absolute_deviation(x, y, cx, cy)
    else:
        x_rel, y_rel = wicmath.relative_deviation(x, y, cx, cy)

    return x_rel, y_rel


def __plot_alone(X, Y):
    f, ax_plt = plt.subplots(1, 1)

    for x_data, y_data, y_legend in zip(X['data'], Y['data'], Y['legend']):
        ax_plt.plot(x_data, y_data, label=y_legend)

    ax_plt.set_xscale(X['scale'])
    ax_plt.set_yscale(Y['scale'])
    ax_plt.set_xlabel(X['label'])
    ax_plt.set_ylabel(Y['label'])

    ax_plt.legend(loc='center', bbox_to_anchor=(0.5, -0.09),
                  ncol=3, fancybox=True, shadow=True)

    plt.tight_layout()

    return 1


def __plot_compare(X, Y):
    f, (ax_plt, ax_dev) = plt.subplots(2, 1, sharex='col')

    if len(X['data']) == len(X['cdata']):
        color = __colors(len(X['data']))

        zipped = zip(X['data'], Y['data'], X['cdata'], Y['cdata'], Y['legend'],
                     Y['clegend'])

        for x_data, y_data, x_cdata, y_cdata, y_legend, y_clegend in zipped:

            c = next(color)
            ax_plt.plot(x_data, y_data, color=c, label=y_legend)
            ax_plt.plot(x_cdata, y_cdata, '--', color=c, label=y_clegend)

            x_rel, y_rel = __deviation(x_data, y_data, x_cdata, y_cdata,
                                       Y['dev'])

            y_rel[abs(y_rel) > Y['rd_max']] = np.nan
            ax_dev.plot(x_rel, y_rel, color=c)

    else:

        for x_data, y_data, y_legend in zip(X['data'], Y['data'], Y['legend']):
            ax_plt.plot(x_data, y_data, label=y_legend)
            x_rel, y_rel = __deviation(x_data, y_data, X['cdata'], Y['cdata'],
                                       Y['dev'])

            y_rel[abs(y_rel) > Y['rd_max']] = np.nan
            ax_dev.plot(x_rel, y_rel)

        ax_plt.plot(X['cdata'], Y['cdata'], '--', label=Y['clegend'])

    ax_plt.set_xscale(X['scale'])
    ax_plt.set_yscale(Y['scale'])
    ax_plt.set_ylabel(Y['label'])

    ax_dev.set_xscale(X['scale'])
    ax_dev.set_yscale(Y['cscale'])
    ax_dev.set_xlabel(X['label'])
    ax_dev.set_ylabel(Y['clabel'])

    size = len(X['data'])

    # TODO: Tune better the factors controlling legend's position.
    ax_plt.legend(loc='center', bbox_to_anchor=(0.5, (-0.03 - 0.05 * size)),
                  ncol=2, fancybox=True, shadow=True, prop={'size': 10})

    plt.tight_layout(h_pad=1 + 1 * size)

    return 1


def __try_loadtxt(filename, min_x, usecols):
    """Try loading data. If variable labels were wrong, print the possible
    options and exit"""
    try:
        x, y = np.loadtxt(__filter_data(filename, usecols[0], min_x),
                          usecols=usecols, unpack=True)

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


def plot_columns(filename, x='', y='', min_x=0, x_scale='log', y_scale='log',
                 x_label='', y_label='', y_legend='', compare_with='',
                 compare_with_legend='', compare_with_scale='linear',
                 rd_max=np.inf, dev='rel', output=''):

    """Given a CLASS output file (or list of files) plot its selected variables
    x, y. If variable labels are not set or misspelled, print the possible
    ones"""
    # TODO: Add support to understand input like y='p_smg/rho_smg'
    X = {'var':  x, 'min': min_x, 'scale': x_scale}
    Y = {'var':  y, 'scale': y_scale, 'rd_max': rd_max, 'cscale':
         compare_with_scale, 'dev': dev}

    Y['clabel'] = 'rel. dev. [%]' if 'abs' not in dev else 'abs. dev.'

    X['data'], Y['data'], Y['legend'] = [], [], []

    if type(filename) is str:
        filename = [filename]
    if type(y_legend) is str:
        y_legend = [y_legend]

    y_legend.extend([''] * (len(filename) - len(y_legend)))
    # If subtraction <= 0 y_legend remains as it is.

    # TODO: Write function for both loops.

    for f, legend in zip(filename, y_legend):  # Zip clips element not paired.
        var_col_dic = chd.class_headers_to_dict(f)

        usecols = (var_col_dic[x], var_col_dic[y])  # If x/y is not a valid key, it
                                                    # returns an empty list. (see
                                                    # chd module)
        x_data, y_data = __try_loadtxt(f, min_x, usecols)
        X['data'].append(x_data)
        Y['data'].append(y_data)
        Y['legend'].append(legend or f.split('/')[-1])  # Label for legend

    X['label'] = x_label or vtd.var_title_dic[x]
    Y['label'] = y_label or vtd.var_title_dic[y]

    if not compare_with:
        __plot_alone(X, Y)
    else:
        # Both headers do not need to be equal, just have the plotting
        # variables.
        #
        # if not check_two_header_lists(filename[0], compare_with):
        #     sys.exit("File headers don't coincide")
        if type(compare_with) is str:
            compare_with = [compare_with]
        if type(compare_with_legend) is str:
            compare_with_legend = [compare_with_legend]

        compare_with_legend.extend([''] * (len(compare_with) -
                                           len(compare_with_legend)))

        for f, legend in zip(compare_with, compare_with_legend):
            var_col_dic = chd.class_headers_to_dict(f)

            usecols = (var_col_dic[x], var_col_dic[y])  # If x/y is not a valid key, it
                                                        # returns an empty list. (see
                                                        # chd module)
            x_data, y_data = __try_loadtxt(f, min_x, usecols)
            X['cdata'].append(x_data)
            Y['cdata'].append(y_data)
            Y['clegend'].append(legend or f.split('/')[-1])  # Label for legend

        __plot_compare(X, Y)

    __plot_close(output)


def plot_w(filename, x='z', w_add=0, min_x=0, x_scale='log', y_scale='log',
           x_label='', y_label='', y_legend='', rd_max=np.inf, dev='rel',
           compare_with_scale='linear', output='', theory=''):
    """Plot the equation of state from a background CLASS output. It compares
    the result of p_smg/rho_smg and p(phi)_smg/rho(phi)_smg"""
    X = {'var':  x, 'min': min_x, 'scale': x_scale}
    Y = {'var':  'w', 'scale': y_scale, 'label': 'w', 'rd_max': rd_max, 'dev':
         dev, 'cscale': compare_with_scale}

    X['label'] = x_label or vtd.var_title_dic[x]
    Y['clabel'] = 'rel. dev. [%]' if 'abs' not in dev else 'abs. dev.'

    X['data'], Y['data'], Y['legend'] = [], [], []
    X['cdata'], Y['cdata'], Y['clegend'] = [], [], []

    if type(filename) is str:
        filename = [filename]
    if type(y_legend) is str:
        y_legend = [y_legend]

    if any('background' not in f.split('/')[-1] for f in filename):
        raise ValueError("File must contain background quantities. If it is " +
                         "indeed a background file, include in its name " +
                         "'background' (e.g. quintessence_background.dat).")

    y_legend.extend([''] * (len(filename) - len(y_legend)))

    for f, legend in zip(filename, y_legend):

        var_col_dic = chd.class_headers_to_dict(f)
        usecols = (var_col_dic[x], var_col_dic['p_smg'], var_col_dic['rho_smg'])

        data = __filter_data(f, usecols[0], min_x)

        x_data, p_smg, rho_smg = np.loadtxt(data, usecols=usecols, unpack=True)
        w = np.divide(p_smg, rho_smg)

        X['data'].append(x_data)
        Y['data'].append(np.add(w, w_add))

        s = legend or f.split('/')[-1]
        Y['legend'].append('p/rho [' + s + ']')  # Label for legend

        data = __filter_data(f, usecols[0], min_x)

        if theory == 'quintessence':
            cw, clegend = quintessence.w(data, var_col_dic)

        Y['cdata'].append(np.add(cw, w_add))
        Y['clegend'].append(clegend + '\n [' + s + ']')

    X['cdata'] = X['data']

    if w_add:
        Y['var'] = str(w_add) + ' + ' + Y['var']
        Y['label'] = str(w_add) + ' + ' + Y['label']

    __plot_compare(X, Y)

    __plot_close(output)


def plot_w0_wa(filename, min_z=False, x_scale='linear', y_scale='linear',
               x_label='', y_label='', y_legend='', rd_max=np.inf, output='',
               theory=''):
    """Plot the evolution of wa with w0 where w(a) \simeq w0 + (a0-a) w'(a0)."""
    # TODO: Implement acceptance of various files. Or eliminate function if
    # useless...

    if 'background' not in filename.split('/')[-1]:
        raise ValueError("File must contain background quantities. If it is " +
                         "indeed a background file, include in its name " +
                         "'background' (e.g. quintessence_background.dat)")

    X = {'scale': x_scale, 'label': 'w0'}
    Y = {'scale': y_scale, 'label': 'wa', 'legend': filename.split('/')[-1],
         'rd_max': rd_max}

    var_col_dic = chd.class_headers_to_dict(filename)
    data = __filter_data(filename, var_col_dic['z'], min_z)

    if theory == 'quintessence':
        w0, wa = quintessence.w0_wa(data, var_col_dic)

    X['data'], Y['data'] = w0, wa

    __plot_alone(X, Y)

    __plot_close(output)
