#!/usr/bin/python

from matplotlib import pyplot as plt
import numpy as np
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D


def __filter_outranged_x(filename, usecols, mins, maxs):
    """Filter values lower than x_min in x_col"""

    with open(filename, 'r') as f:
        for line in f:
            if '#' not in line:
                c = 0
                for col, x_min, x_max in zip(usecols, mins, maxs):
                    x = float(line.split()[col])
                    if x < x_min or x > x_max:
                        break
                    c += 1

                if c == len(usecols):
                    yield line


def plot(filename, x=0, y=1, x_abs=False, y_abs=False,
         x_scale='linear', y_scale='linear', x_label='x', y_label='y',
         scatter=False, size=12, title='', xmin=-np.inf, xmax=np.inf,
         ymin=-np.inf, ymax=np.inf):

    filtercols = (x, y)
    mins = (xmin, ymin)
    maxs = (xmax, ymax)
    data_filtered = __filter_outranged_x(filename, filtercols, mins, maxs)

    xs, ys = np.loadtxt(data_filtered, usecols=(x, y), unpack=True)

    if x_abs:
        xs = np.fabs(xs)
    if y_abs:
        ys = np.fabs(ys)

    f, ax = plt.subplots(1, 1)
    ax.set_xscale(x_scale)
    ax.set_yscale(y_scale)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)

    if title:
        plt.title(title)

    if scatter:
        ax.scatter(xs, ys, s=size)
    else:
        ax.plot(xs, ys, linewidth=size)

    plt.show()
    plt.close()


def plot_color(filename, x=0, y=1, c=2, x_abs=False, y_abs=False,
               x_scale='linear', y_scale='linear', x_label='col_0',
               y_label='col_1', z_label='None', c_label='', size=30, title='',
               vmax=None, vmin=None, cmap=1, xmin=-np.inf, xmax=np.inf,
               ymin=-np.inf, ymax=np.inf):

    filtercols = (x, y)
    mins = (xmin, ymin)
    maxs = (xmax, ymax)
    data_filtered = __filter_outranged_x(filename, filtercols, mins, maxs)

    fig, ax = plt.subplots()
    xs, ys, cs = np.loadtxt(data_filtered, usecols=(x, y, c), unpack=True)

    if x_abs:
        xs = np.fabs(xs)
    if y_abs:
        ys = np.fabs(ys)

    colormap = {
        1: cm.coolwarm,
        2: cm.plasma,
        3: cm.Blues,
        4: cm.afmhot,
    }

    cax = ax.scatter(xs, ys, s=size, c=cs, cmap=colormap[cmap], vmax=vmax,
                     vmin=vmin)

    ax.set_title(title)
    ax.set_xscale(x_scale)
    ax.set_yscale(y_scale)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)

    cbar = fig.colorbar(cax)
    cbar.set_label(c_label)

    plt.show()
    plt.close()


def plot_3d(filename, x=0, y=1, z=2, c=None, x_abs=False, y_abs=False,
            z_abs=False, x_scale='linear', y_scale='linear', z_scale='linear',
            x_label='col_0', y_label='col_1', z_label='col_2', c_label='',
            size=30, title='', vmax=None, vmin=None, cmap=1, xmin=-np.inf,
            xmax=np.inf, ymin=-np.inf, ymax=np.inf, zmin=-np.inf, zmax=np.inf,
            depthshade=0, xy_plot=False, line_ini=None, line_end=None,
            line_label=None):

    filtercols = (x, y, z)
    mins = (xmin, ymin, zmin)
    maxs = (xmax, ymax, zmax)
    data_filtered = __filter_outranged_x(filename, filtercols, mins, maxs)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    colormap = {
        1: cm.coolwarm,
        2: cm.plasma,
        3: cm.Blues,
        4: cm.afmhot,
    }

    if c is not None:
        cmap = colormap[cmap]
        xs, ys, zs, cs = np.loadtxt(data_filtered, usecols=(x, y, z, c),
                                    unpack=True)
    else:
        cmap = None
        xs, ys, zs = np.loadtxt(data_filtered, usecols=(x, y, z), unpack=True)

    if x_abs:
        xs = np.fabs(xs)
    if y_abs:
        ys = np.fabs(ys)
    if z_abs:
        zs = np.fabs(zs)

    if c is not None:
        cax = ax.scatter(xs, ys, zs, s=size, c=cs, cmap=cmap, vmax=vmax,
                         vmin=vmin, depthshade=depthshade)
        cbar = fig.colorbar(cax)
        cbar.set_label(c_label)
    else:
        cax = ax.scatter(xs, ys, zs, s=size)

    if xy_plot:
        zmin_aux = ax.set_zlim3d()[0]
        ax.scatter(xs, ys, zmin_aux, s=3, color='g')

    # Consider allowing any kind of line. ATM just straight lines.
    if (line_ini is not None) and (line_end is not None):
        line = zip(line_ini, line_end)
        ax.plot(line[0], line[1], line[2], linewidth=0.05 * size, c='black',
                label=line_label)
    elif ((line_ini is not None) or
          (line_end is not None)) and ((line_ini is None) or (line_end is
                                                              None)):
        print "Line not drawn. You must enter the initial and end points"

    ax.set_title(title)
    ax.set_xscale(x_scale)
    ax.set_yscale(y_scale)
    ax.set_zscale(z_scale)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_zlabel(z_label)

    ax.legend()

    plt.show()
    plt.close()
