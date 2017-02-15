#!/usr/bin/python

#from classy import Class
from matplotlib import pyplot as plt
import numpy as np

######
# This program must be versatile. It will accept classy dictionaries to
# compare. It should be use together with classy.
######


class Compare():
    def __init__(self, x, *args):
        if self.__check_dicts(*args): self._dicts = args
        if self.__check_x(x): self._x= x

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

    def update_x(self, x):
        self.__check_x(x)
        self._x = x

    def update_dicts(self, *args):
        self.__check_dicts(*args)
        self._dicts = args
        self.__check_x()

    def relative_deviation(self):
        pass

    def absolute_deviation(self):
        pass

    def plot(self, which="all"):
        pass



