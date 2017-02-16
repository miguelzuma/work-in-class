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

    def __select_key(self, y):
        if y == 'all':
            keys = self.__refdict.iterkeys()
        else:
            keys = [y]

        return keys

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

    def plot(self, y="all"):
        pass

