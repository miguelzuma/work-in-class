#!/usr/bin/python

from matplotlib import pyplot as plt
import numpy as np


def plot(filename, x=0, y=1, x_abs=False, y_abs=False, x_scale='linear',
         y_scale='linear', x_label='x', y_label='y', scatter=False, size=10,
         title=''):

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
