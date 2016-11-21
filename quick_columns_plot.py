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


def plot_color(filename, x=0, y=1, c=2, x_abs=False, y_abs=False,
               x_scale='linear', y_scale='linear', x_label='col_0',
               y_label='col_1', c_label='', size=14, title='', vmax=None,
               vmin=None):

    x, y, c = np.loadtxt(filename, usecols=(x, y, c), unpack=True)

    if x_abs:
        x = np.fabs(x)
    if y_abs:
        y = np.fabs(y)

    fig, ax = plt.subplots()

    cax = ax.scatter(x, y, s=size, c=c, cmap=cm.coolwarm, vmax=vmax, vmin=vmin)
    ax.set_title(title)
    ax.set_xscale(x_scale)
    ax.set_yscale(y_scale)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)

    # Add colorbar, make sure to specify tick locations to match desired
    # ticklabels
    # cbar = fig.colorbar(cax, ticks=[-7, 0, 6])
    # cbar.ax.set_yticklabels(['-7', '0', '7'])

    cbar = fig.colorbar(cax)
    cbar.set_label(c_label)

    plt.show()
