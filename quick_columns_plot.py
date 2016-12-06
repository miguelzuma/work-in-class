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
         xscale='linear', yscale='linear', xlabel='x', ylabel='y',
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
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    if title:
        plt.title(title)

    if scatter:
        ax.scatter(xs, ys, s=size)
    else:
        ax.plot(xs, ys, linewidth=size)

    plt.axis('tight')
    plt.show()
    plt.close()


def plot_color(filename, x=0, y=1, c=2, x_abs=False, y_abs=False,
               xscale='linear', yscale='linear', xlabel='col_0',
               ylabel='col_1', zlabel='None', clabel='', size=30, title='',
               vmax=None, vmin=None, cmap=1, xmin=-np.inf, xmax=np.inf,
               ymin=-np.inf, ymax=np.inf, density=True, bins=None,
               gridsize=100):

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

    if density:
        cax = ax.hexbin(xs, ys, bins=bins, gridsize=gridsize, xscale=xscale,
                        yscale=yscale)
    else:
        cax = ax.scatter(xs, ys, s=size, c=cs, cmap=colormap[cmap], vmax=vmax,
                         vmin=vmin)
        ax.set_xscale(xscale)
        ax.set_yscale(yscale)

    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    cbar = fig.colorbar(cax)
    cbar.set_label(clabel)

    plt.axis('tight')
    plt.show()
    plt.close()


def plot_3d(filename, x=0, y=1, z=2, c=None, x_abs=False, y_abs=False,
            z_abs=False, xscale='linear', yscale='linear', zscale='linear',
            xlabel='col_0', ylabel='col_1', zlabel='col_2', clabel='',
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
        cbar.set_label(clabel)
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
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    ax.set_zscale(zscale)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)

    ax.legend()

    plt.show()
    plt.close()
