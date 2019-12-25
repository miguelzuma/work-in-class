#!/usr/bin/python

from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from scipy.special import factorial
from scipy import stats
import numpy as np
import warnings


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


def find_nearest(array, value, reldev=True):
    """
    Return index of the array element closest to value.
    """
    if reldev:
        idx = (np.abs(1 - array / value)).argmin()
    else:
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

    return x[:-2], 0.5/h[:-1] * (-3. * y[:-2] + 4. * y[1:-1] - y[2:])


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


def log_modulus(data):
    """
    Log10-modulus transform data.
    Inspired by: https://blogs.sas.com/content/iml/2014/07/14/log-transformation-of-pos-neg.html
    """

    return np.sign(data) * np.log10(np.abs(data) + 1)


def log_modulus_inverse(data):
    """
    Log10-modulus inverse function transform data.
    """

    return np.sign(data) * (10.**(np.abs(data)) - 1)

def log_complex(data, regularize=True):
    """
    Return log10(x) if x>0 and i*log(x) if x<0.

    If regularize is True, change 0 by 1e-100
    """
    dat = np.array(data)
    output = np.array(data, dtype=complex)

    if regularize:
        # TODO: In order not to oversize the positive values, consider
        # multiplying by a (-1)^n term, with n drawn from a binomial
        # distribution with even probability.
        warnings.warn("Changing 0's to 1e-100")
        output[dat == 0] = (1e-100 + 0j)

    output = np.log10(np.abs(output), dtype=complex)
    output[dat < 0] *= 1j

    return output

def log_complex_inverse(data, threshold=1e-3):
    """
    Return the inverse function of log_complex.
    """
    #TODO: Optimize this!!!
    flat_data = np.array(data).flatten()
    shape_data = np.shape(data)
    output = np.zeros(len(flat_data))

    for i, value in enumerate(flat_data):
        if np.abs(value.real/value.imag) < threshold:
            output[i] = -10.**value.imag
        elif np.abs(value.imag/value.real) < threshold:
            output[i] = 10. ** value.real
        else:
            raise ValueError('Log_complex transformation cannot be inverted ' + \
                             'since both parts of item #{} ({}) are allowed '.format(i, value) + \
                             'by threshold: try increasing it')

    return output.reshape(shape_data)


def sample_log_data(dat, bins):
    """
    Return data sampled in the logarithmic-scaled distribution of dat.
    """
    H = []
    for i in [np.log10(dat[dat > 0]), np.log10(-dat[dat < 0]), dat[dat == 0]]:
        h, x = np.histogram(i, bins=bins)
        H.append(stats.rv_histogram([h,x]))

    probs1 = len(dat[dat > 0]) / float(len(dat))
    probs2 = len(dat[dat < 0]) / float(len(dat))

    data = []

    for counter in dat:
        sign = stats.uniform.rvs()
        if sign < probs1:
            data.append(10**(H[0].rvs()))
        elif (sign >= probs1) & (sign <= (probs1 + probs2)):
            data.append(-10**(H[1].rvs()))
        else:
            data.append(0)

    return np.array(data)


def pade(an, m):
    """
    Return Pade approximation coefficients of a polynomial.
    Parameters
    ----------
    an : (N,) array_like
        Taylor series coefficients.
    m : int
        The order of the returned approximating polynomials.
    Returns
    -------
    p, q : array
        The pade approximation of the polynomial defined by `an` is
        the polynomial with coefficients `poly(p)/poly(q)`.

    Note
    --------
    This code is taken from scipy.misc.pade distributed under SciPy license (https://www.scipy.org/scipylib/license.html):

    Copyright © 2001, 2002 Enthought, Inc.
    All rights reserved.

    Copyright © 2003-2013 SciPy Developers.
    All rights reserved.

    Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

        Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
        Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
        Neither the name of Enthought nor the names of the SciPy Developers may be used to endorse or promote products derived from this software without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
    """
    from scipy import linalg
    an = np.asarray(an)
    N = len(an) - 1
    n = N - m

    if n < 0:
        raise ValueError("Order of q <m> must be smaller than len(an)-1.")
    Akj = np.eye(N+1, n+1)
    Bkj = np.zeros((N+1, m), 'd')
    for row in range(1, m+1):
        Bkj[row, :row] = -(an[:row])[::-1]
    for row in range(m+1, N+1):
        Bkj[row, :] = -(an[row-m:row])[::-1]
    C = np.hstack((Akj, Bkj))

    pq = linalg.solve(C, an)
    p = pq[:n+1]
    q = np.r_[1.0, pq[n+1:]]

    return p, q[1:]  # remove the denominator's 1 of (1 + ... )


"""
Fit padde wrapper functions by Emilio Bellini.
"""


def pade_approx(x, coeff_num, coeff_den):  # the actual fit function
    num = 0.
    den = 1.
    for i, coeff in enumerate(coeff_num):
        num = num + coeff*(x**i)
    for i, coeff in enumerate(coeff_den):
        den = den + coeff*(x**(i+1))
    return num/den


def _p0_calc(xdata, ydata, n_num, n_den):

    len_der = 2*(n_num+n_den)+1

    Y = ydata[:len_der]

    taylor_coeff = [Y[0]]

    for i in range(n_num+n_den):
        _, Y = diff_three_points(xdata[:len(Y)], Y)
        taylor_coeff.append(Y[0]/np.math.factorial(i+1))

    for i in range(n_den + 1):
        try:
            p, q = pade(taylor_coeff[:n_den-i + n_num + 1], n_den-i)
            break
        except np.linalg.LinAlgError:
            continue

    return p, np.concatenate([q, np.zeros(i)])


def _wrapper(x, n_num, n_den, *args):
    coeff_num = list(args[0:n_num+1])
    coeff_den = list(args[n_num+1:n_num+n_den+1])
    return pade_approx(x, coeff_num, coeff_den)


def fit_pade(xdata, ydata, n_num, n_den, p0=[], **kwards):

    if p0 == []:
        p0 = _p0_calc(xdata, ydata, n_num, n_den)
        p0 = np.concatenate(p0)

    popt, _ = curve_fit(lambda x, *p0: _wrapper(xdata, n_num, n_den, *p0), xdata, ydata, p0=p0, **kwards)
    yfit = _wrapper(xdata, n_num, n_den, *popt)
    return popt, yfit

"""
End of fit padde wrapper functions by Emilio Bellini.
"""


def Pade(x, c):
    """
    Pade. the first half of the coeffs will be in the numerator. The rest in the
    denominator.

    Coeffs are scaled by the factorial of the acompanying monomial magnitude order.
    """
    l = len(c)/2
    coeff_num = c[:l] / factorial(range(len(c[:l])), exact=True)
    coeff_den = c[l:] / factorial(range(len(c[l:])), exact=True)
    return pade_approx(x, coeff_num, coeff_den)


def Taylor_legacy(x, c):
    """
    Taylor expansion without factorial terms
    """
    s = 0
    for i, ci in enumerate(c):
        s += ci * x**i
    return s


def Taylor(x, c):
    s = 0
    for i, ci in enumerate(c):
        s += ci / np.math.factorial(i) * x**i
    return s

def Taylor_log(x, c):
    s = 0
    for i, ci in enumerate(c):
        s += np.sign(ci) * (10.**np.abs(1./ci) - 1) / np.math.factorial(i) * x**i
    return s

def Taylor_c0(x, c):
    return Taylor(x, np.concatenate([[0.], c]))

def Taylor_legacy_c0(x, c):
    return Taylor_legacy(x, np.concatenate([[0.], c]))

def Taylor_c1(x, c):
    return Taylor(x, np.concatenate([[1.], c]))

def Taylor_legacy_c1(x, c):
    return Taylor_legacy(x, np.concatenate([[1.], c]))

def bins_equally_spaced(x, c):
    h, x = np.histogram(x, bins=len(c))
    result = []

    for i, bin in enumerate(x[:-1]):
        result += [c[i]] * h[i]

    return np.array(result).flatten()


def _wrapper_fit(function):
    def wrapper_out(x, n_coeffs, *args):
        coeffs = list(args)
        return function(x, coeffs)
    return wrapper_out


def fit(function, xdata, ydata, n_coeffs, p0=[], **kargs):
    if p0 == []:
        p0 = np.ones(n_coeffs)
    popt, _ = curve_fit(lambda x, *p0: _wrapper_fit(function)(xdata, n_coeffs, *p0), xdata, ydata, p0=p0, **kargs)
    yfit = _wrapper_fit(function)(xdata, n_coeffs, *popt)
    return popt, yfit
