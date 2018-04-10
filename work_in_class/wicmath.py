#!/usr/bin/python

import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit


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


def find_nearest(array, value):
    """
    Return index of the array element closest to value.
    """
    idx = (np.abs(array - value)).argmin()

    return idx


def diff(x, y):
    """
    Output X, dy/dx, where X are the intermediate values of the x array points.

    x = array of values which y change respect to.
    y = array of values to differentiate with respect to x.
    """
    dx = np.diff(x)
    dydx = np.diff(y) / dx
    X = x[:-1] + dx / 2.

    return X, dydx


def diff_three_points(x, y):
    """
    Output f'(x) = 1/2h [-3 f(x) + 4 f(x+h) - f(x+2h)]

    x = array of values which y change respect to. x should be equispaced.
    y = array of values to differentiate with respect to x.
    """

    h = np.diff(x)

    return x[:-3], 0.5/h[:-3] * (-3. * y[:-3] + 4. * y[1:-2] - y[2:])


def intermediate(array):
    """
    Return an array with the intermediate values of array
    """
    if type(array) is not np.ndarray:
        array = np.array(array)

    return array[:-1] + np.diff(array) / 2.


def reflect_data(data):
    """
    Reflect data at its lowest value.
    """
    data_rf = []
    for dat in data:
        x = np.sort(dat)

        x_rf = 2 * x[0] - x[x > x[0]]  # = x0 - (x - x0), x excluding x0

        X = np.concatenate([x_rf, dat])

        data_rf.append(X)

    return data_rf


"""
Fit padde wrapper functions by Emilio Bellini.
"""

def fit_pade(xdata, ydata, n_num, n_den):
    def pade_approx(x, coeff_num, coeff_den):  # the actual fit function
        num = 0.
        den = 1.
        for i, coeff in enumerate(coeff_num):
            num = num + coeff*(x**i)
        for i, coeff in enumerate(coeff_den):
            den = den + coeff*(x**(i+1))
        return num/den

    def wrapper(x, n_num, n_den, *args):
        coeff_num = list(args[0:n_num+1])
        coeff_den = list(args[n_num+1:n_num+n_den+1])
        return pade_approx(x, coeff_num, coeff_den)

    p0 = np.ones(n_num+n_den+1)
    popt, _ = curve_fit(lambda x, *p0: wrapper(xdata, n_num, n_den, *p0), xdata, ydata, p0=p0)
    yfit = wrapper(xdata, n_num, n_den, *popt)
    return popt, yfit

"""
End of fit Pade wrapper.
"""
