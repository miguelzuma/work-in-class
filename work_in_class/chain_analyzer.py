#!/usr/bin/python

import numpy as np
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
import re


class Chain():
    """
    Class to store and analyze MCMC chains.
    """
    def __init__(self):
        self.chains = np.array([])
        self.CosmoHammerArguments = {}
        self.paramsNames = []
        self.paramsScaleByType = {'cosmo': [],
                                  'nuisance': [],
                                  'derived': []}
        self.paramsByType = {'cosmo': [],
                             'nuisance': [],
                             'derived': []}
        self.paramsFixedByType = {'cosmo': {},
                                  'nuisance': {}}
        self.cosmoArguments = {}

    def empty(self):
        self.chains = np.array([])

    def _readCosmoHammerOptions(self, fileArguments):
        """
        Fill CosmoHammerArguments with walkersRatio, burninIterations and
        sampleIterations values
        """
        with open(fileArguments) as f:
            for line in f:
                argument, value = line.split('=')
                self.CosmoHammerArguments[argument.strip()] = int(value)

    def readMCMCParameters(self, logparamPath):
        """
        Read the name of MCMC parameters and distinguish between 'cosmo',
        'nuisance' and 'derived'.
        """
        with open(logparamPath) as f:
            for line in f:
                if ('=' not in line) or (line.strip()[0] == '#'):
                    continue

                if ('data.parameters' in line):
                    a = line.split("'")
                    paramType = a[-2]
                    paramName = a[1]

                    arguments = line.split("[")[-1].split(',')

                    if (float(arguments[3]) == 0) and (paramType != 'derived'):
                        value = float(arguments[0]) * float(arguments[4])
                        self.paramsFixedByType[paramType].update({paramName:
                                                                  value})
                    else:
                        self.paramsNames.append(paramName)
                        self.paramsByType[paramType].append(paramName)
                        self.paramsScaleByType[paramType].append(float(arguments[4]))

                elif ('data.cosmo_arguments' in line):
                    paramName = line.split("'")[1]
                    paramValue = line.split("=")[1]
                    try:
                        self.cosmoArguments.update({paramName: float(paramValue)})
                    except:
                        s = paramValue.strip().replace("'", '').replace('"', '')
                        self.cosmoArguments.update({paramName: s})

    def readCosmoHammerChain(self, chainPaths, fileArguments, numberFreeParam, removeError=False):
        """
        Append individual walkers chains in self.chains array
        """
        self._readCosmoHammerOptions(fileArguments)

        walkers = self.CosmoHammerArguments['walkersRatio'] * numberFreeParam

        chain = np.loadtxt(chainPaths[0])
        for cpath in chainPaths[1:]:
            chain = np.vstack((chain, np.loadtxt(cpath)))

        chainsWalker = np.array([chain[i::walkers] for i in range(walkers)])

        if removeError:
            chains = []
            for walker in chainsWalker:
                chains.append(walker[~np.isnan(walker).any(axis=1)])
            self.chains = np.array(chains)

        else:
            self.chains = chainsWalker

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
                multiplicity = int(columns[0])
                columns[0] = 1  # Set multiplicity to 1 as they will be repeated
                for n in range(multiplicity):
                    chain.append(columns[firstColumnToStore:])

        # self.chains.append(chain)
        self.chains = np.append(self.chains, chain, axis=0)

        # TODO: Think about removing first two columns or adding them to CH
        # (more problematic? Think if this program is going to be only for test
        # convergence so lkl is not necessarily required)

    def autocorrelation_acor(self, paramIndex, walker=None):
        """
        Returns the autocorrelation fucntion for paramIndex.

        walker = None or int. If int, return the autocorrelation function for
        paramIndex for walker.

        Uses acor
        """
        import acor

        if walker:
            acor.function(self.chains[walker, :, paramIndex])
        else:
            mean = 0
            walkers = len(self.chains)
            for i in range(walkers):
                mean += acor.function(self.chains[i, :, paramIndex])
            return mean/float(walkers)

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

    def plotAutocorrelation(self, paramIndex, acor=True, normalized=False, GoodmanWeare=False, withExponential=False, **kwards):
        """
        Plot the autocorrelation function. It fails if some chains are removed.
        """

        n = len(self.chains[0])
        if type(paramIndex) is not list:
            paramIndex = [paramIndex]

        for index in paramIndex:
            if acor:
                autocorrelation = self.autocorrelation_acor(index)
            elif not GoodmanWeare:
                autocorrelation = [self.autocorrelation(T, index) for T in range(n)]
            else:
                autocorrelation = [self.autocorrelationGoodmanWeare(T, index) for T in range(n)]

            if normalized and (not acor):
                autocorrelation = np.array(autocorrelation) / autocorrelation[0]

            plt.plot(range(n), autocorrelation, label=index)

            if withExponential:
                # It returns a list: [(tau, sigma)]
                tau = self.getExponentialAutocorrelationTime(index,
                                                             GoodmanWeare=GoodmanWeare,
                                                             **kwards)[0][0]
                plt.plot(range(n), autocorrelation[0] * np.exp(-np.arange(n)/tau),
                         '--', label='Fit {}'.format(index))

        plt.xlabel('Lag')
        plt.ylabel('Autocorrelation')
        plt.legend()
        plt.show()
        plt.close()

    def getGelmanRubin(self, MontePython=False):
        """
        Return an array with the R-1 (Gelman-Rubin criterion without the
        t-Student dimensional estimator correction) for each parameter as in
        MontePython.

        MontePython = Bool. If true, return R - 1 computed (almost) as in
        MontePython, i.e. R = Between / Within variances; except for the
        weighting for the chain length.

        """
        varOfWalker, meanOfWalker = [], []

        if self.CosmoHammerArguments:
            startingIteration = self.CosmoHammerArguments['burninIterations']
        else:
            startingIteration = 0

        for walker in self.chains:
            varOfWalker.append(np.var(walker[startingIteration:], axis=0, ddof=1))
            meanOfWalker.append(np.mean(walker[startingIteration:], axis=0))

        W = np.mean(varOfWalker, axis=0)
        b = np.var(meanOfWalker, axis=0, ddof=1)  # b = B/n
        n = float(len(self.chains[0]))
        m = float(len(self.chains))

        V = W - W/n + b * (1 + 1. / m)

        if MontePython:
            return V/W - 1

        return np.sqrt(V / W) - 1

    def getIntegratedAutocorrelationTime(self, paramIndex):
        """
        Return the integrated autocorrelation time.

        tau = 1/2 + sum C(T)/C(0)
        """

        C0 = self.autocorrelationGoodmanWeare(0, paramIndex)

        tau_sum_array = np.array([self.autocorrelationGoodmanWeare(T, paramIndex)
                                  for T in range(len(self.chains[0]))])/C0

        return 0.5 + tau_sum_array.sum()

    def getIntegratedAutocorrelationTime_acor(self, paramIndex):
        """
        Return the walkers average integrated autocorrelation time.
        It uses acor.
        """
        import acor

        tauMean = 0

        for walker in self.chains:
            paramChain = [row[paramIndex] for row in walker]
            tau, mean, sigma = acor.acor(paramChain)
            tauMean += tau

        return tauMean/len(self.chains)

    def getMaxIntegratedAutocorrelationTime(self, params=None, acor=False):
        """
        Return the greates integrated autocorrelation time
        """

        if params is None:
            params = len(self.chains[0][0])

        if acor:
            callAutocorrelationTime = self.getIntegratedAutocorrelationTime_acor
        else:
            callAutocorrelationTime = self.getIntegratedAutocorrelationTime

        tau = [callAutocorrelationTime(paramIndex) for paramIndex in
               range(params)]

        return max(tau)

    def _getExponentialAutocorrelationTime(self, paramIndex, tauGuess=1, GoodmanWeare=False):
        """
        Return the exponential autocorrelation time for the given paramIndex and
        its error.

        It is obtained fitting exp(-T/tau) to data up to C(T) < e^-1.
        """

        if GoodmanWeare:
            callAutocorrelation = self.autocorrelationGoodmanWeare
        else:
            callAutocorrelation = self.autocorrelation

        T = [0]
        C = [callAutocorrelation(0, paramIndex)]

        for t in range(1, len(self.chains[0])):
            c = callAutocorrelation(t, paramIndex)
            if c/C[0] > np.exp(-1):
                C.append(c)
                T.append(t)
            else:
                break

        tau, cov = curve_fit(lambda t, tau: np.exp(-t/tau), T, np.array(C)/C[0],
                             p0=(tauGuess))

        # tau, cov are np.arrays, let's return just the value
        return (T, C), (tau[0], np.sqrt(cov[0][0]))

    def getExponentialAutocorrelationTime(self, paramIndex, tauGuess=1,
                                          GoodmanWeare=False):
        """
        Return the exponential autocorrelation time for the given paramIndexes
        and its error.

        It is obtained fitting exp(-T/tau) to data up to C(T) < e^-1.
        """

        if type(paramIndex) is not list:
            paramIndex = [paramIndex]

        tau = []

        for index in paramIndex:
            tau.append(self._getExponentialAutocorrelationTime(index, tauGuess,
                                                               GoodmanWeare)[1])

        return tau

    def getStepFailed(self, cosmo=False):
        """
        Return steps for which computation failed (lkl = -inf).

        cosmo = False. If True return just the cosmological parameters.
        """
        result = [[]] * len(self.chains)

        if cosmo:
            selection = [param in self.paramsByType['cosmo'] for param in self.paramsNames]
        else:
            selection = [True] * len(self.chains[0][0])

        for i, walker in enumerate(self.chains):
            result[i] = np.array([step[selection] for step in walker if True in np.isnan(step)])

        return np.array(result)

    def getCosmoParamsStepFailed(self):
        """
        Return the classy dictionary to compute the steps that failed.
        """
        return self.getCosmoParamsFromWalkersSteps(self.getStepFailed(cosmo=True))

    def getCosmoParamsFromWalkersSteps(self, walkersSteps):
        """
        Return the classy dictionary to compute the steps inputted.

        walkersSteps must be an array with shape (#walkers, #steps, #params)
        """
        cosmoParams = walkersSteps

        params = [[]] * len(cosmoParams)

        for i, walker in enumerate(cosmoParams):
            for step in walker:
                params[i].append(self.getCosmoParamsFromOneStep(step))

        return params

    def getCosmoParamsFromOneStep(self, step):
        """
        Return the classy dictionary to compute the step inputted.
        """
        cosmoNames = self.paramsByType['cosmo']
        cosmoScales = self.paramsScaleByType['cosmo']
        cosmoFixedNames = [name for name in self.paramsFixedByType['cosmo'].iterkeys()]

        allNames = cosmoNames + cosmoFixedNames

        parameters_smgNames = filter(re.compile('parameters_smg').match, allNames)
        parameters_2_smgNames = filter(re.compile('parameters_2_smg').match, allNames)

        parameters_smgNames.sort()
        parameters_2_smgNames.sort()

        d = {key: val * scale for key, val, scale in zip(cosmoNames, step, cosmoScales)}
        d.update(self.paramsFixedByType['cosmo'])

        parameters_smg = []
        for name in parameters_smgNames:
            parameters_smg.append(d[name])
            del d[name]
        if parameters_smg:
            d['parameters_smg'] = str(parameters_smg).strip('[]')

        parameters_2_smg = []
        for name in parameters_2_smgNames:
            parameters_2_smg.append(d[name])
            del d[name]
        if parameters_2_smg:
            d['parameters_2_smg'] = str(parameters_smg).strip('[]')

        return d
