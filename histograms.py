#!/usr/bin/python

import numpy as np
import sys
from matplotlib import pyplot as plt
from matplotlib import cm


class Histograms():
    def __init__(self):
        self.fname = ''
        self.data = []
        self.histograms = []
        self.bins = []
        self.means = []
        self.sigmas = []

    def read_file(self, fname, header=True):
        """
        Read the file with the binned data.

        header = True or False. If True read extract bins from file header.
        """
        self.fname = fname

        self.data = np.loadtxt(fname, unpack=True)

        if header:
            with open(self.fname) as f:
                f.readline()
                header = f.readline()

            self.bins = map(float, header.strip('# ').split())

    def _compute_histograms(self, bins):
        """
        Compute the histograms of the data

        Output: [[hist, bin_edges], ...]
        """
        histograms = np.empty(len(self.data)).tolist()

        for i, data in enumerate(self.data):
            histograms[i] = np.histogram(data, bins=bins)  # histogram[i] = [hist, bin_edges]

        self.histograms = histograms

    def _compute_means(self):
        """
        Compute the mean for each bin.
        """
        self.means = np.mean(self.data, axis=1)

    def _compute_sigmas(self):
        """
        Compute the intervals that contains the 68% of the data around the means.
        """
        sigmas = np.empty(len(self.data)).tolist()

        for index, dataRaw in enumerate(self.data):
            data = dataRaw.copy()
            data.sort()
            mean = self.means[index]
            sigmaData = int(0.34*len(data) + 0.5)
            mean_index = np.abs(data-mean).argmin()
            try:
                sigma_up = data[mean_index + sigmaData] - mean
            except:
                sys.stderr.write("Upper 68% limmit not well defined for bin {}\n".format(index))
                if (mean_index + sigmaData - len(data)) == 0:
                    sys.stderr.write("Making an approximation\n")
                    sigma_up = data[mean_index + sigmaData - 1] - mean
                else:
                    sigma_up = np.nan
            try:
                sigma_down = data[mean_index - sigmaData] - mean
            except:
                sys.stderr.write("Lower 68% limmit not well defined for bin {}\n".format(index))
                if (mean_index - sigmaData) == -1:
                    sys.stderr.write("Making an approximation\n")
                    sigma_down = data[mean_index + sigmaData + 1] - mean
                else:
                    sigma_down = np.nan

            sigmas[index] = np.array([np.abs(sigma_down), sigma_up])  # Error are positive

        self.sigmas = sigmas

    def compute(self, bins='auto'):
        """
        Compute histograms, means and sigmas
        """
        self._compute_histograms(bins)
        self._compute_means()
        self._compute_sigmas()

    def plot_histogram(self, index, variable_binned='bin_i', xlabel='x', xlim=[None, None]):
        """
        Plot histogram for bin with index, index
        """
        data = self.data[index].copy()
        data.sort()
        histogram = self.histograms[index]
        mean = self.means[index]
        sigma = self.sigmas[index]

        bins = histogram[1][:-1]
        density = histogram[0]/histogram[0].sum(dtype=float)

        plt.bar(bins, density, width=np.diff(histogram[1]), color='b')
        plt.plot(bins, density, c='r')
        plt.axvline(mean, c='black', label='mean')
        plt.axvspan(mean - sigma[0], mean + sigma[1], alpha=0.5, color='grey', label='68% data')

        plt.xlim(xlim)
        plt.xlabel(xlabel)
        plt.legend(loc=0)
        plt.title('{} = {}'.format(variable_binned, self.bins[index]))
        plt.show()
        plt.close()

    def plot_evolution(self, xlabel='x', ylabel='y'): #, scale=['linear', 'linear']):
        """
        Plot evolution.
        """
        freq , y = zip(*self.histograms)
        prob = np.array(freq).T/np.sum(freq, dtype=float, axis=1)
        prob = prob.T
        y = np.delete(y, -1, 1)

        x = self.bins*np.ones(np.shape(y)).T
        x = x.T

        cax = plt.scatter(x, y, c=prob, s=2, cmap=cm.hot)

        cb = plt.colorbar(cax)
        cb.set_label('% data')
        plt.ylabel(ylabel)
        plt.xlabel(xlabel)
        plt.show()
        plt.close()

    def reset(self):
        """
        Reset the default values (erase everything)
        """
        self.__init__()
