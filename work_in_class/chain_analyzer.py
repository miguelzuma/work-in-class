#!/usr/bin/python

import numpy as np


class Chain():
    """
    Class to store and analyze MCMC chains.
    """
    def __init__(self):
        self.chains = []
        self.CosmoHammerArguments = {}

    def _readCosmoHammerOptions(self, fileArguments):
        """
        Fill CosmoHammerArguments with walkersRatio, burninIterations and
        sampleIterations values
        """
        with open(fileArguments) as f:
            for line in f:
                argument, value = line.split('=')
                self.CosmoHammerArguments[argument.strip()] = int(value)

    def readCosmoHammerChain(self, filepath, fileArguments, numberFreeParam, removeError=False):
        """
        Append individual walkers chains in self.chains array
        """
        self._readCosmoHammerOptions(fileArguments)

        if 'burn' in filepath:
            iterKey = 'burninIterations'
        else:
            iterKey = 'sampleIterations'

        with open(filepath) as f:
            walkers = self.CosmoHammerArguments['walkersRatio'] * numberFreeParam

            self.chains = [[] for i in range(walkers)]

            for iteration in range(self.CosmoHammerArguments[iterKey]):
                for walker in range(walkers):
                    line = f.next().split()
                    if ('nan' in line) and removeError:
                        print "Removing point with lkl=-np.inf. Be aware it can " +\
                               "make the rest of operations meaningless."
                        continue
                    self.chains[walker].append(np.array([float(x) for x in line]))

    def readChain(self, filepath):
        """
        Append chain to self.chains object
        """
        self.chains.append(np.loadtxt(filepath))

    def autocorrelation(self, T, paramIndex):
        """
        Compute autocorrelation following arXiv:1212.1721v2 Appendix C.
        """

        C = 0

        walkers = len(self.chains)

        for i in range(walkers):
            C += self.autocorrelationPerWalker(T, i, paramIndex)

        return C/walkers

    def autocorrelationPerWalker(self, T, walker, paramIndex):
        """
        Compute the autocorrelation for a walker
        """
        chain = zip(*self.chains[walker])[paramIndex]
        n = len(chain)
        mean = np.mean(chain)

        chaint = np.array(chain[:n-T])
        chainForward = np.array(chain[T:])

        return np.sum((chaint - mean)*(chainForward - mean))/(n-T)
