#!/usr/bin/python

#from classy import Class
from matplotlib import pyplot as plt
import numpy as np
import wicmath

######
# This program must be versatile. It will accept classy dictionaries to
# compare. It should be used together with classy.
######


class Compare():
    def __init__(self, x, *args):
        if self.__check_dicts(*args): self._dicts = args
        if self.__check_x(x): self._x= x
        self._adev = ()
        self._rdev = ()
        self.__refdict = self._dicts[0]

    def __check_x(self, *args):
        if args:
            x = args[0]
        else:
            x = self._x

        for adict in self._dicts:
            if x not in adict.viewkeys():
                raise ValueError(str(x) + " must be a key in ALL input dictionaries.\n \
                                 Use .update_x() to change x or update_dics to use different dictionaries")

        return True

    def __check_dicts(self, *args):
        dicts = args or self.__dicts

        for adict in dicts:
            if not dicts[0].viewkeys() == adict.viewkeys():
                print "Dictionaries do not have the same keys some functions may not work."
                # return False
        return True

    def __deviation(self, y, kind):
        ref = self.__refdict
        x = self._x
        keys = self.__select_key(y)

        if kind == 'rel':
            deviation = wicmath.relative_deviation
        elif kind == 'abs':
            deviation = wicmath.absolute_deviation
        else:
            raise ValueError("kind must be one of 'rel or 'abs'")

        output_dicts = []

        for adict in self._dicts[1:]:
            d = {}
            for akey in keys:
                d[akey] =  deviation(adict[x], adict[akey], ref[x], ref[akey])
            output_dicts.append(d)

        return tuple(output_dicts)

    def __plot(self, y, dicts, kind='', label='', xabs=False, yabs=False,
               xscale='linear', yscale='linear', xlabel='', ylabel='',
               title='', xlim=[None, None], ylim=[None, None], xsum=0, ysum=0):

        keys = self.__select_key(y)

        color = plt.cm.rainbow(np.linspace(0, 1, len(dicts)))

        for akey in keys:
            i = 0
            f, ax = plt.subplots(1, 1)

            for adict in dicts:
                if kind == "normal":
                    X, Y = adict[self._x], adict[akey]
                else:
                    X, Y = adict[akey]

                X = X + xsum
                Y = Y + ysum

                Y_pos = Y.copy()
                Y_neg = Y.copy()

                Y_pos[Y_pos < 0] = np.nan
                Y_neg[Y_neg >= 0] = np.nan

                label2 = label

                if np.any(np.isfinite(Y_pos)):
                    plt.plot(X, Y_pos, c=color[i], label=label.format(2+i, i))
                    label2 = ""  # This way second plt.plot won't produce label
                if np.any(np.isfinite(Y_neg)):
                    plt.plot(X, np.abs(Y_neg), '--', c=color[i], label=label2.format(2+i, i))
                i += 1

            ax.set_xscale(xscale)
            try:
                ax.set_yscale(yscale)
            except:
                print "Error logging-scale. Using linear scale for y axis"
                ax.set_yscale('linear')
                pass

            ax.set_xlabel(xlabel or str(self._x) + " + " + str(xsum))
            ax.set_ylabel(ylabel or akey + " [Dev] + " + str(ysum))
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)

            plt.title(title or kind)
            plt.legend(loc="best")
            plt.show()
            plt.close()

    def __remove_self(self, local):
        d = local.copy()
        del d['self']

        return d

    def __select_key(self, y):
        if y == 'all':
            keys = self.__refdict.iterkeys()
        elif hasattr(y, '__iter__'):
            keys = y
        else:
            keys = [y]

        return keys

    def relative_deviation(self, y='all'):
        """Returns the percentual relative deviation between the common
        dictionaries variables or the specified variable if y is given. The
        fist dictionary is taken as reference"""
        self._rdev = self.__deviation(y, 'rel')
        return self._rdev

    def absolute_deviation(self, y='all'):
        """Returns the absolute deviation between the common dictionaries
        variables or the specified variable if y is given"""
        self._adev = self.__deviation(y, 'abs')
        return self._adev

    def update_x(self, x):
        self.__check_x(x)
        self._x = x

    def update_dicts(self, *args):
        self.__check_dicts(*args)
        self._dicts = args
        self.__check_x()

    def plot_rdev(self, y='all', xabs=False, yabs=False, xscale='linear',
                  yscale='linear', xlabel='', ylabel='', title='',
                  xlim=[None, None], ylim=[None, None], xsum=0, ysum=0):
        dicts = self._rdev or self.relative_deviation()
        label = 'd1-d{0} rdev.'
        self.__plot(**self.__remove_self(locals()))
        return True

    def plot_adev(self, y='all', xabs=False, yabs=False, xscale='linear',
                  yscale='linear', xlabel='', ylabel='', title='',
                  xlim=[None, None], ylim=[None, None], xsum=0, ysum=0):
        dicts = self._adev or self.absolute_deviation()
        label = 'd1-d{0} adev.'
        self.__plot(**self.__remove_self(locals()))
        return True

    def plot(self, y='all', xabs=False, yabs=False, xscale='linear',
             yscale='linear', xlabel='', ylabel='', title='',
             xlim=[None, None], ylim=[None, None], xsum=0, ysum=0):
        kind = 'normal'
        dicts = self._dicts
        label = 'd{1}'
        self.__plot(**self.__remove_self(locals()))
        return True


if __name__ == "__main__":
    pass
