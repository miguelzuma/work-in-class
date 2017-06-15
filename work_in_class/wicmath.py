#!/usr/bin/python

import numpy as np
from scipy.interpolate import interp1d


def __deviation(x, y, cx, cy, kind):

    if len(x) == len(cx) and not (x - cx).all():
        abs_dev = np.subtract(y, cy)
        cY = cy
        X = x

    else:

        cf = interp1d(cx, cy)

        # Create an x-array contained by the interpolation boundaries
        X = x[x >= min(cx)]
        X = X[X <= max(cx)]

        # Use the y-values of the X-array
        b1 = np.where(x == X[0])[0][0]
        b2 = np.where(x == X[-1])[0][0]

        Y = y[b1:b2 + 1]
        cY = cf(X)

        abs_dev = np.subtract(Y, cY)

    if kind == 'abs':
        return X, abs_dev
    elif kind == 'rel':
        return X, np.divide(abs_dev, cY)
    else:
        raise ValueError("Deviation kind must be 'abs' or 'rel'")


def relative_deviation(x, y, cx, cy):  # c stands for compare
    """Return relative deviation between [cx, cy] and [x, y]."""

    X, rel_dev = __deviation(x, y, cx, cy, 'rel')

    return X, np.multiply(rel_dev, 100)


def absolute_deviation(x, y, cx, cy):
    """Return absolute deviation between [cx, cy] and [x, y]."""

    X, abs_dev = __deviation(x, y, cx, cy, 'abs')

    return X, abs_dev


def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx
