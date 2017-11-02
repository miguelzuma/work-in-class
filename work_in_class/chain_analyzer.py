#!/usr/bin/python

import numpy as np
from matplotlib import pyplot as plt


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

    def readChain(self, filepath, removeFirstCols=True):
        """
        Append chain to self.chains object
        """

        chain = []
        firstColumnToStore = 2

        if not removeFirstCols:
            firstColumnToStore = 0

        with open(filepath) as f:
            for line in f:
                columns = [float(x) for x in line.split()]
                multiplicity = columns[0] = 0
                columns[0] = 1  # Set multiplicity to 1 as they will be repeated
                for n in range(multiplicity):
                    chain.append(columns[firstColumnToStore:])

        self.chains.append(chain)

        # TODO: Think about removing first two columns or adding them to CH
        # (more problematic? Think if this program is going to be only for test
        # convergence so lkl is not necessarily required)

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

    def autocorrelationGoodmanWeare(self, T, paramIndex):
        """
        Compute the autocorrelation as Goodman and Weare (2010).

        It cannot be used if some step has been removed.
        """

        n = len(self.chains[0])

        Ft = [self._ensembleAveragePerIteration(iteration, paramIndex)
              for iteration in range(n-T)]

        Fforward = [self._ensembleAveragePerIteration(iteration, paramIndex)
                    for iteration in range(T, n)]

        Fmean = np.mean(np.concatenate((Ft, Fforward)))

        return np.sum((Ft - Fmean)*(Fforward - Fmean))/(n-T)

    def _ensembleAveragePerIteration(self, iteration, paramIndex):
        """
        Return the ensemble average per iteration
        """

        value = 0

        for walker in self.chains:
            value += walker[iteration][paramIndex]

        return value / len(self.chains)

    def plotAutocorrelation(self, paramIndex, GoodmanWeare=False):
        """
        Plot the autocorrelation function. It fails if some chains are removed.
        """
        n = len(self.chains[0])
        if type(paramIndex) is not list:
            paramIndex = [paramIndex]

        for index in paramIndex:
            if not GoodmanWeare:
                autocorrelation = [self.autocorrelation(T, index) for T in range(n)]
            else:
                autocorrelation = [self.autocorrelationGoodmanWeare(T, index) for T in range(n)]

            plt.plot(range(n), autocorrelation, label=index)
        plt.xlabel('Lag')
        plt.ylabel('Autocorrelation')
        plt.legend()
        plt.show()
        plt.close()

    def getGelmanRubin(self):
        """
        Return an array with the R-1 (original Gelman-Rubin criterion) for each
        parameter

        From Brooks & Gelman (2017)
        """
        varOfWalker, meanOfWalker = [], []

        for walker in self.chains:
            varOfWalker.append(np.var(walker, axis=0, ddof=1))
            meanOfWalker.append(np.mean(walker, axis=0))

        W = np.mean(varOfWalker, axis=0)
        b = np.var(meanOfWalker, axis=0, ddof=1)  # b = B/n
        n = len(self.chains[0])
        m = len(self.chains)

        Var = W - W/n + b * (1 + 1. / m)

        return np.sqrt(Var/W) - 1
