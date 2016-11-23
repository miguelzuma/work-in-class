#!/usr/bin/python

from matplotlib import pyplot as plt
import numpy as np
from matplotlib import cm


def plot(filename, x=0, y=1, x_abs=False, y_abs=False,
         x_scale='linear', y_scale='linear', x_label='x', y_label='y',
         scatter=False, size=12, title=''):

    x, y = np.loadtxt(filename, usecols=(x, y), unpack=True)

    if x_abs:
        x = np.fabs(x)
    if y_abs:
        y = np.fabs(y)

    f, ax = plt.subplots(1, 1)
    ax.set_xscale(x_scale)
    ax.set_yscale(y_scale)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)

    if title:
        plt.title(title)

    if scatter:
        ax.scatter(x, y, s=size)
    else:
        ax.plot(x, y, linewidth=size)

    plt.show()
    plt.close()


def plot_color(filename, x=0, y=1, c=2, x_abs=False, y_abs=False,
               x_scale='linear', y_scale='linear', x_label='col_0',
               y_label='col_1', z_label='None', c_label='', size=30, title='',
               vmax=None, vmin=None, cmap=1, xmin=None, xmax=None, ymin=None,
               ymax=None):

    fig, ax = plt.subplots()
    x, y, c = np.loadtxt(filename, usecols=(x, y, c), unpack=True)

    if x_abs:
        x = np.fabs(x)
    if y_abs:
        y = np.fabs(y)

    colormap = {
        1: cm.coolwarm,
        2: cm.plasma,
        3: cm.Blues,
        4: cm.afmhot,
    }

    cax = ax.scatter(x, y, s=size, c=c, cmap=colormap[cmap], vmax=vmax,
                     vmin=vmin)

    ax.set_title(title)
    ax.set_xscale(x_scale)
    ax.set_yscale(y_scale)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)

    cbar = fig.colorbar(cax)
    cbar.set_label(c_label)

    plt.axis([xmin, xmax, ymin, ymax])

    plt.show()
    plt.close()


def plot_3d(filename, x=0, y=1, z=2, c=None, x_abs=False, y_abs=False,
            z_abs=False, x_scale='linear', y_scale='linear', z_scale='linear',
            x_label='col_0', y_label='col_1', z_label='None', c_label='',
            size=30, title='', vmax=None, vmin=None, cmap=1, xmin=None,
            xmax=None, ymin=None, ymax=None, zmin=None, zmax=None,
            depthshade=0, xy_plot=False, line_ini=None, line_end=None,
            line_label=None):

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
        x, y, z, c = np.loadtxt(filename, usecols=(x, y, z, c), unpack=True)
    else:
        cmap = None
        x, y, z = np.loadtxt(filename, usecols=(x, y, z), unpack=True)

    if x_abs:
        x = np.fabs(x)
    if y_abs:
        y = np.fabs(y)
    if z_abs:
        z = np.fabs(z)

    cax = ax.scatter(x, y, z, s=size, c=c, cmap=cmap, vmax=vmax, vmin=vmin,
                     depthshade=depthshade)

    if xy_plot:
        zmin_aux = zmin or ax.set_zlim3d()[0]
        ax.plot(x, y, zmin_aux)

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

    cbar = fig.colorbar(cax)
    cbar.set_label(c_label)

    ax.set_xlim3d([xmin, xmax])
    ax.set_ylim3d([ymin, ymax])
    ax.set_zlim3d([zmin, zmax])

    ax.legend()

    plt.show()
    plt.close()
