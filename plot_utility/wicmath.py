#!/usr/bin/python

import numpy as np
from scipy.interpolate import interp1d


def relative_deviation(x, y, cx, cy):  # c stands for compare
    """Return relative deviation between [cx, cy] and [x, y]."""

    if len(x) == len(cx) and not (x - cx).all:
        rel_dev = np.divide(np.substract(y, cy), cy)

    else:

        cf = interp1d(cx, cy)

        # Create an x-array contained by the interpolation boundaries
        X = x[x > max(x[0], cx[0])]
        X = X[X < min(x[-1], cx[-1])]

        # Use the y-values of the X-array
        b1 = np.where(x == X[0])[0][0]
        b2 = np.where(x == X[-1])[0][0]

        Y = y[b1:b2 + 1]
        cY = cf(X)

        rel_dev = np.divide(np.subtract(Y, cY), cY)

    return X, np.multiply(rel_dev, 100)
