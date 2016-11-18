#!/usr/bin/python

from matplotlib import pyplot as plt
import numpy as np


def plot(filename, x=0, y=1, x_abs=False, y_abs=False, x_scale='linear',
         y_scale='linear', x_label='x', y_label='y', scatter=False, size=10):

    x, y = np.loadtxt(filename, usecols=(x, y), unpack=True)

    if x_abs:
        x = np.fabs(x)
    if y_abs:
        y = np.fabs(y)

    f, ax = plt.subplots(1, 1)
    ax.x_scale(x_scale)
    ax.y_scale(y_scale)
    ax.x_label(x_label)
    ax.y_label(y_label)

    if scatter:
        ax.scatter(x, y, s=size)
    else:
        ax.plot(x, y, linewidth=size)

    plt.show()

