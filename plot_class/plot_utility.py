#!/usr/bin/python

import class_headers_to_dict as chd
import matplotlib.pyplot as plt
import variable_title_dics as vtd
import wicmath
# from list_class_headers import check_two_header_lists
import quintessence
import model_independent_variables as miv
import numpy as np
import matplotlib.cm as cm
import sys


def __exit_variable_error(filename):
    title_var_dic = chd.class_file_title_var_dic(filename)
    title_var_str = ""
    for title, var in title_var_dic.iteritems():
        title_var_str += ''.join("'" + title + "'" + ': ' + "'" + var +
                                 "'" + '\n')

    my_err = ("You probably used wrong variable labels. See bellow " +
              "the avaibles ones for the present file:" +
              filename.split('/')[-1] + '\n' + title_var_str)

    sys.exit(my_err)


def __prepare_input_for_loop(filename, y_legend, IC={}):
    if type(filename) is str:
        filename = [filename]
    if type(y_legend) is str:
        y_legend = [y_legend]

    y_legend.extend([''] * (len(filename) - len(y_legend)))
    # If subtraction <= 0 y_legend remains as it is.

    if IC:
        if type(IC) is dict:
            IC = [IC]
        return filename, y_legend, IC
    else:
        return filename, y_legend


def __check_variables_in_files(filename, x, y=''):
    for f in filename:
        try:
            var_col_dic = chd.class_headers_to_dict(f)
            1 + var_col_dic[x]  # 1 + [] if x is not a key
            if y:
                1 + var_col_dic[y]

        except TypeError:
            __exit_variable_error(f)


def __init_dictionaries(var, x_boundaries, labels, scales, cscale, rd_max, dev,
                        xy_add=[0, 0], IC={}):

    x_add, y_add = xy_add

    X = {
        'var': var[0],
        'min': x_boundaries[0],
        'max': x_boundaries[1],
        'scale': scales[0],
        'label':  labels[0] or vtd.var_title_dic[var[0]],
        'add': xy_add[0]
    }

    Y = {
        'var': var[1],
        'scale': scales[1],
        'label':  labels[1] or vtd.var_title_dic[var[1]],
        'clabel': 'rel. dev. [%]' if 'abs' not in dev else 'abs. dev.',
        'cscale': cscale,
        'rd_max': rd_max,
        'dev': dev,
        'add': xy_add[1],
        'IC': IC
    }

    if x_add:
        X['label'] = str(x_add) + ' + ' + X['label']
    if y_add:
        Y['label'] = str(y_add) + ' + ' + Y['label']

    return X, Y


def __load_columns_file_loop(X, Y, filename, y_legend, func, compare=False):

    kdata = 'data'
    klegend = 'legend'

    x_add, y_add = X['add'], Y['add']

    if compare:
        kdata = 'c' + kdata
        klegend = 'c' + klegend

    X[kdata], Y[kdata], Y[klegend] = [], [], []

    for f, lgd in zip(filename, y_legend):  # Zip clips elements not paired.
        var_col_dic = chd.class_headers_to_dict(f)
        filter_col = var_col_dic[X['var']]
        data = __filter_data(f, filter_col, X['min'], X['max'])
        x_data, y_data, legend = func(X['var'], data, var_col_dic)

        X[kdata].append(np.add(x_add, x_data))
        Y[kdata].append(np.add(y_add, y_data))

        s = lgd or f.split('/')[-1]
        Y[klegend].append(legend + ' [' + s + ']')  # Label for legend

    return 1


def __load_var_data_legend_loop(X, Y, filename, y_legend, func, cfunc):

    X['data'], Y['data'], Y['legend'] = [], [], []
    X['cdata'], Y['cdata'], Y['clegend'] = [], [], []

    x_add, y_add = X['add'], Y['add']
    IC = Y['IC']

    if len(IC) != len(filename):
        sys.exit("There must be as much IC as input files.")

    for f, lgd, ic in zip(filename, y_legend, IC):

        var_col_dic = chd.class_headers_to_dict(f)
        filter_col = var_col_dic[X['var']]
        data = __filter_data(f, filter_col, X['min'], X['max'])

        x_data, y_data, legend = func(X['var'], data, var_col_dic)

        X['data'].append(np.add(x_add, x_data))
        Y['data'].append(np.add(y_add, y_data))

        s = lgd or f.split('/')[-1]
        Y['legend'].append(legend + ' [' + s + ']')  # Label for legend

        data = __filter_data(f, filter_col, X['min'], X['max'])

        y_cdata, clegend = cfunc(data, var_col_dic, ic)

        Y['cdata'].append(np.add(y_add, y_cdata))
        Y['clegend'].append(clegend + '\n [' + s + ']')

    X['cdata'] = X['data']

    return 1


def __filter_outranged_x(filename, x_col, x_min, x_max):
    """Filter values lower than x_min in x_col"""

    with open(filename, 'r') as f:
        for line in f:
            if '#' not in line:
                x = float(line.split()[x_col])
                if x >= x_min and x <= x_max:
                    yield line


def __filter_data(filename, col, x_min, x_max):
    """In case x_min is provided, filter all values lower than it"""
    def __data_gen(filename):
        with open(filename) as f:
            for line in f:
                yield line

    if x_min == -np.inf and x_max == np.inf:
        return __data_gen(filename)

    else:
        return __filter_outranged_x(filename, col, x_min, x_max)


def __deviation(x, y, cx, cy, kind):

    if 'abs' in kind:
        x_rel, y_rel = wicmath.absolute_deviation(x, y, cx, cy)
    else:
        x_rel, y_rel = wicmath.relative_deviation(x, y, cx, cy)

    return x_rel, y_rel


def __print_and_close(output, Y):
    ax_plt = Y['axes'][0]
    size = len(Y['data'])

    if output:
        # TODO: Is there a way not to repeat ourselves using ax_plt.legend?
        # TODO: Tune better the factors controlling legend's position.

        ax_plt.legend(loc='center', bbox_to_anchor=(0.5, (-0.3 - 0.06 * size)),
                      ncol=2, fancybox=True, shadow=True, prop={'size': 10})
        plt.tight_layout(h_pad=3 + 1.5 * size)
        plt.savefig(output, bbox_inches='tight', dpi=300)

    else:
        plt.show()

    plt.close()


def __plot_alone(X, Y):
    f, ax_plt = plt.subplots(1, 1)
    Y['axes'] = [ax_plt]  # For compatibility with __print_and_close

    for x_data, y_data, y_legend in zip(X['data'], Y['data'], Y['legend']):
        ax_plt.plot(x_data, y_data, label=y_legend)

    ax_plt.set_xscale(X['scale'])
    ax_plt.set_yscale(Y['scale'])
    ax_plt.set_xlabel(X['label'])
    ax_plt.set_ylabel(Y['label'])

    # TODO: Improve legend position management.
    ax_plt.legend(loc='center', bbox_to_anchor=(0.5, -0.09),
                  ncol=3, fancybox=True, shadow=True)

    plt.tight_layout()

    return 1


def __plot_compare_dirty_job(data, cdata, legends, color, Y):
    x_data, y_data = data
    x_cdata, y_cdata = cdata
    y_legend, y_clegend = legends
    ax_plt, ax_dev = Y['axes']

    c = next(color)
    ax_plt.plot(x_data, y_data, color=c, label=y_legend)

    if len(Y['cdata']) == len(Y['data']):  # TODO: Double check. __p_c and here.
        ax_plt.plot(x_cdata, y_cdata, '--', color=c, label=y_clegend)

    x_rel, y_rel = __deviation(x_data, y_data, x_cdata, y_cdata,
                               Y['dev'])
    y_rel[abs(y_rel) > Y['rd_max']] = np.nan
    ax_dev.plot(x_rel, y_rel, color=c)


def __plot_compare(X, Y):
    def __colors(length):
        color = iter(cm.rainbow(np.linspace(0, 1, length)))

        return color

    f, Y['axes'] = plt.subplots(2, 1, sharex='col')
    ax_plt, ax_dev = Y['axes']

    color = __colors(len(Y['data']))

    if len(X['data']) == len(X['cdata']):
        zipped = zip(X['data'], Y['data'], X['cdata'], Y['cdata'], Y['legend'],
                     Y['clegend'])

        for x_data, y_data, x_cdata, y_cdata, y_legend, y_clegend in zipped:
            __plot_compare_dirty_job([x_data, y_data], [x_cdata, y_cdata],
                                     [y_legend, y_clegend], color, Y)

    elif len(X['cdata']) == 1:

        x_cdata, y_cdata = X['cdata'][0], Y['cdata'][0]
        y_clegend = Y['clegend'][0]

        for x_data, y_data, y_legend in zip(X['data'], Y['data'], Y['legend']):
            __plot_compare_dirty_job([x_data, y_data], [x_cdata, y_cdata],
                                     [y_legend, y_clegend], color, Y)

        ax_plt.plot(x_cdata, y_cdata, '--', label=y_clegend)

    else:
        raise ValueError("compare_with must be of the same size as filename")

    ax_plt.set_xscale(X['scale'])
    ax_plt.set_yscale(Y['scale'])
    ax_plt.set_ylabel(Y['label'])

    ax_dev.set_xscale(X['scale'])
    ax_dev.set_yscale(Y['cscale'])
    ax_dev.set_xlabel(X['label'])
    ax_dev.set_ylabel(Y['clabel'])

    size = len(X['data'])
    # TODO: Tune better the factors controlling legend's position.
    ax_plt.legend(loc='center', bbox_to_anchor=(0.5, (-0.03 - 0.05 *
                                                      size)), ncol=2,
                  fancybox=True, shadow=True, prop={'size': 10})
    plt.tight_layout(h_pad=1 + 1 * size)

    return 1


def plot_columns(filename, x='', y='', x_add=0, y_add=0, x_min=-np.inf,
                 x_max=np.inf, x_scale='log', y_scale='log', x_label='',
                 y_label='', y_legend='', compare_with='',
                 compare_with_legend='', compare_with_scale='linear',
                 rd_max=np.inf, dev='rel', output=''):
    """Given a CLASS output file (or list of files) plot its selected variables
    x, y. If variable labels are not set or misspelled, print the possible
    ones"""
    # TODO: Add support to understand input like y='p_smg/rho_smg'
    filename, y_legend = __prepare_input_for_loop(filename, y_legend)
    __check_variables_in_files(filename, x, y)
    X, Y = __init_dictionaries([x, y], [x_min, x_max], [x_label, y_label],
                               [x_scale, y_scale], compare_with_scale, rd_max,
                               dev, xy_add=[x_add, y_add])

    columns = miv.columns(x, y)
    __load_columns_file_loop(X, Y, filename, y_legend, columns)

    if not compare_with:
        __plot_alone(X, Y)
    else:
        cw, cwl = __prepare_input_for_loop(compare_with, compare_with_legend)

        __load_columns_file_loop(X, Y, cw, cwl, columns, compare=True)

        __plot_compare(X, Y)

    __print_and_close(output, Y)


def plot_w0_wa(filename, z_min=-np.inf, z_max=np.inf, x_scale='linear',
               y_scale='linear', x_label='', y_label='', y_legend='',
               rd_max=np.inf, output='', theory='', IC={}):
    """Plot the evolution of wa with w0 where w(a) \simeq w0 + (a0-a) w'(a0)."""
    # TODO: Implement acceptance of various files. Or eliminate function if
    # useless...

    X = {'scale': x_scale, 'label': 'w0'}
    Y = {'scale': y_scale, 'label': 'wa', 'rd_max': rd_max}

    X['data'], Y['data'], Y['legend'] = [], [], []

    filename, y_legend = __prepare_input_for_loop(filename, y_legend)

    if theory == 'quintessence':
        w0_wa = quintessence.w0_wa

    for f, lgd in zip(filename, y_legend):  # Zip clips elements not paired.
        var_col_dic = chd.class_headers_to_dict(f)
        # If x/y is not a valid key, it returns an empty list. (see chd module)
        data = __filter_data(f, var_col_dic['z'], z_min, z_max)
        x_data, y_data = w0_wa(data, var_col_dic)

        X['data'].append(x_data)
        Y['data'].append(y_data)
        Y['legend'].append(lgd or f.split('/')[-1])  # Label for legend

    __plot_alone(X, Y)

    __print_and_close(output, Y)


def plot_check_w(filename, x='z', x_add=0, y_add=0, x_min=-np.inf,
                 x_max=np.inf, x_scale='log', y_scale='log', x_label='',
                 y_label='', y_legend='', rd_max=np.inf, dev='rel',
                 compare_with_scale='linear', output='', theory='', IC={}):
    """Plot the equation of state from a background CLASS output. It compares
    the result of p_smg/rho_smg and p(phi)_smg/rho(phi)_smg"""

    filename, y_legend, IC = __prepare_input_for_loop(filename, y_legend, IC)
    __check_variables_in_files(filename, x)

    y = 'w'
    y_label = y_label or 'w'
    X, Y = __init_dictionaries([x, y], [x_min, x_max], [x_label, y_label],
                               [x_scale, y_scale], compare_with_scale, rd_max,
                               dev, xy_add=[x_add, y_add], IC=IC)

    w = miv.w

    cw = quintessence.w(theory)

    __load_var_data_legend_loop(X, Y, filename, y_legend, w, cw)

    __plot_compare(X, Y)

    __print_and_close(output, Y)


def plot_check_alphaK(filename, x='z', x_add=0, y_add=0, x_min=-np.inf,
                      x_max=np.inf, x_scale='log', y_scale='log', x_label='',
                      y_label='', y_legend='', rd_max=np.inf, dev='rel',
                      compare_with_scale='linear', output='', theory='', IC={},
                      choice=1, kineticity_safe=0):
    """Plot alpha_k from a background CLASS output and compare it with other
    calculation of alpha_k"""

    filename, y_legend, IC = __prepare_input_for_loop(filename, y_legend, IC)
    __check_variables_in_files(filename, x)

    y = 'alpha_K'
    y_label = y_label or 'alpha_k'
    X, Y = __init_dictionaries([x, y], [x_min, x_max], [x_label, y_label],
                               [x_scale, y_scale], compare_with_scale, rd_max,
                               dev, xy_add=[x_add, y_add], IC=IC)

    alphaK = miv.alphaK(kineticity_safe)

    calphaK = quintessence.alphaK(choice, theory)

    __load_var_data_legend_loop(X, Y, filename, y_legend, alphaK, calphaK)

    __plot_compare(X, Y)

    __print_and_close(output, Y)


def plot_w(filename, x='z', x_add=0, y_add=0, x_min=-np.inf, x_max=np.inf,
           x_scale='log', y_scale='log', x_label='', y_label='', y_legend='',
           compare_with='', compare_with_legend='', rd_max=np.inf, dev='rel',
           compare_with_scale='linear', output='', IC={}):
    """Given a CLASS background output file (or list of files) plot the
    equation of state vs x (default 'z'). If x-variable label is misspelled,
    print the possible ones"""
    filename, y_legend, IC = __prepare_input_for_loop(filename, y_legend, IC)
    __check_variables_in_files(filename, x)

    y = 'w'
    y_label = y_label or 'w'
    X, Y = __init_dictionaries([x, y], [x_min, x_max], [x_label, y_label],
                               [x_scale, y_scale], compare_with_scale, rd_max,
                               dev, xy_add=[x_add, y_add], IC=IC)

    w = miv.w

    __load_columns_file_loop(X, Y, filename, y_legend, w)

    if not compare_with:
        __plot_alone(X, Y)
    else:
        cw, cwl = __prepare_input_for_loop(compare_with, compare_with_legend)

        __load_columns_file_loop(X, Y, cw, cwl, w, compare=True)

        __plot_compare(X, Y)

    __print_and_close(output, Y)


def plot_alphaK(filename, x='z', x_add=0, y_add=0, x_min=-np.inf, x_max=np.inf,
                x_scale='log', y_scale='log', x_label='', y_label='',
                y_legend='', rd_max=np.inf, dev='rel',
                compare_with_scale='linear', compare_with='',
                compare_with_legend='', output='', theory='', IC={},
                kineticity_safe=0):
    """Given a CLASS background output file (or list of files) plot alpha_k
    vs x (default 'z'). If x-variable label is misspelled, print the possible
    ones"""

    filename, y_legend, IC = __prepare_input_for_loop(filename, y_legend, IC)
    __check_variables_in_files(filename, x)

    y = 'alpha_K'
    y_label = y_label or 'alpha_k'
    X, Y = __init_dictionaries([x, y], [x_min, x_max], [x_label, y_label],
                               [x_scale, y_scale], compare_with_scale, rd_max,
                               dev, xy_add=[x_add, y_add], IC=IC)

    alphaK = miv.alphaK(kineticity_safe)

    __load_columns_file_loop(X, Y, filename, y_legend, alphaK)

    if not compare_with:
        __plot_alone(X, Y)
    else:
        cw, cwl = __prepare_input_for_loop(compare_with, compare_with_legend)

        __load_columns_file_loop(X, Y, cw, cwl, alphaK, compare=True)

        __plot_compare(X, Y)

    __print_and_close(output, Y)
