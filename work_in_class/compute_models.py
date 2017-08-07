#!/usr/bin/pyhon

"""Module thought to be used when you want to obtain some quantities (e.g. the
background structures) for a model varying one of the parameters"""

from classy import Class
from collections import OrderedDict as od
from matplotlib import pyplot as plt
import inifile_parser as inip
import wicmath
import numpy as np


class Model():
    def __init__(self, cosmo=None):
        """
        Initialize the Model class. By default Model uses its own Class
        instance.

        cosmo = external Class instance. Default is None
        """
        if cosmo:
            self.cosmo = cosmo
        else:
            self.cosmo = Class()
        self.computed = {}

    def __set_scale(self, axes, xscale, yscale):
        """
        Set scales for axes in axes array.

        axes = axes array (e.g. f, ax = plt.subplots(2,2))
        xscale = linear array of xscale.
        yscale = linear array of yscale.

        Scales are set once axes is flatten. Each plot is counted from left to
        right an from top to bottom.
        """
        for i, ax in enumerate(axes.flat):
            ax.set_xscale(xscale[i])
            ax.set_yscale(yscale[i])

    def __set_label(self, axes, xlabel, ylabel):
        """
        Set labels for axes in axes array.

        axes = axes array (e.g. f, ax = plt.subplots(2,2))
        xlabel = linear array of xlabels.
        ylabel = linear array of ylabels.

        Labels are set once axes is flatten. Each plot is counted from left to
        right an from top to bottom.
        """
        for i, ax in enumerate(axes.flat):
            ax.set_xlabel(xlabel[i])
            ax.set_ylabel(ylabel[i])

    def compute_models(self, params, varied_name, index_variable, values,
                       back=[], thermo=[], prim=[], pert=[], trans=[],
                       extra=[], update=True, cosmo_msg=False):
        """
        Fill dic with the background structures for the model with given
        params, modifying the varied_name value with values.

        params = parameters to be set in Class. They must be in agreement with
                what is asked for.
        varied_name = the name of the variable you are modifying. It will be
                      used as key in dic assigned to its background structures.
        index_variable = variable's index in parameters_smg array.
        values = varied variable values you want to compute the cosmology for.
        back = list of variables to store from background. If 'all', store the
              whole dictionary.
        thermo = list of variables to store from thermodynamics. If 'all',
                  store the whole dictionary.
        prim = list of variables to store from primordial. If 'all', store the
               whole dictionary.
        pert = list of variables to store from perturbations. If 'all', store
               the whole dictionary.
        trans = list of variables to store from transfer. If 'all', store
                the whole dictionary. get_transfer accept two optional
                arguments: z=0 and output_format='class' (avaible options are
                'class' or 'camb'). If different values are desired, first
                item of trans must be {'z': value, 'output_format': value}.
        *args = any of the method or objects defined in cosmo, e.g. w0_smg().
                It will store {arg: cosmo.w0_smg()}
        """

        key = varied_name

        if (not update) or (key not in self.computed.keys()):
            self.computed[key] = od()

        for val in values:
            # key = "{}={}".format(varied_name, val)
            params["parameters_smg"] = inip.vary_params(params["parameters_smg"], [[index_variable, val]])

            d = self.computed[key][val] = {}

            self.cosmo.empty()
            self.cosmo.set(params)

            try:
                self.cosmo.compute()
            except Exception, e:
                print "Error: skipping {}={}".format(key, val)
                if cosmo_msg:
                    print e

                continue

            d['tunned'] = self.cosmo.get_current_derived_parameters(['tuning_parameter'])['tuning_parameter']

            for lst in [[back, 'back', self.cosmo.get_background],
                        [thermo, 'thermo', self.cosmo.get_thermodynamics],
                        [prim, 'prim', self.cosmo.get_thermodynamics]]:
                if lst[0]:
                    output = lst[2]()
                    if lst[0][0] == 'all':
                        d[lst[1]] = output
                    else:
                        d[lst[1]] = {}
                        for item in back:
                            if type(item) is list:
                                d[lst[1]].update({item[0]: output[item[0]][item[1]]})
                            else:
                                d[lst[1]].update({item: output[item]})

            if pert:
                # Perturbation is tricky because it can accept two optional
                # argument for get_perturbations and this method returns a
                # dictionary {'kind_of_pert': [{variable: list_values}]}, where
                # each item in the list is for a k (chosen in params).
                if type(pert[0]) is dict:
                    output = self.cosmo.get_perturbations(pert[0]['z'], pert[0]['output_format'])
                    if pert[1] == 'all':
                        d['pert'] = output
                else:
                    output = self.cosmo.get_perturbations()
                    if pert[0] == 'all':
                        d['pert'] = output

                if (type(pert[0]) is not dict) and (pert[0] != 'all'):
                    d['pert'] = {}
                    for subkey, lst in output.iteritems():
                        d['pert'].update({subkey: []})
                        for n, kdic in enumerate(lst):  # Each item is for a k
                            d['pert'][subkey].append({})
                            for item in pert:
                                if type(item) is list:
                                    d['pert'][subkey][n].update({item[0]: kdic[item[0]][item[1]]})
                                else:
                                    d['pert'][subkey][n].update({item: kdic[item]})

            for i in extra:
                exec('d[i] = self.cosmo.{}'.format(i))

            self.cosmo.struct_cleanup()

    def plot_4_vs_z(self, varied_name, y1name, y2name,
                    labelvaried_name, y1label, y2label, z_s=100,
                    yscale=['linear', 'linear', 'linear', 'linear'],
                    xscale=['log', 'log', 'log', 'log'],
                    exclude=[], scatter=False):
        """
        Plot a 2x2 figure, with two variables y1 and y2 vs z+1. The ones in
        the second column are a zoom of the ones in the first one.

        varied_name = varied variable's name
        y1name = output class name of the variable for Y axis 1
        y2name = output class name of the variable for Y axis 2
        labelvaried_name = label for varied_name
        y1label = label for Y1 axis
        y2label = label for Y2 axis
        z_s = greatest z to plot in zoomed figures
        xscale = list of scales for the x axis
        yscale = list of scales for the y axis
        exclude = list of the varied variable values to exclude from plotting.
        """
        fig, ax = plt.subplots(2, 2, figsize=(15, 10))
        xlabel = 'z+1'

        module = 'back'

        if scatter:
            ax00plot = ax[0, 0].scatter
            ax10plot = ax[1, 0].scatter
            ax01plot = ax[0, 1].scatter
            ax11plot = ax[1, 1].scatter
        else:
            ax00plot = ax[0, 0].plot
            ax10plot = ax[1, 0].plot
            ax01plot = ax[0, 1].plot
            ax11plot = ax[1, 1].plot

        for i, ba in self.computed[varied_name].iteritems():
            if (i in exclude) or (not ba):
                continue
            x = ba[module]['z'] + 1
            y1 = ba[module][y1name]
            y2 = ba[module][y2name]
            z_i = wicmath.find_nearest(ba[module]['z'], z_s)

            label_i = labelvaried_name + '={}'.format(i)

            ax00plot(x, y1, label=label_i)
            ax10plot(x, y2, label=label_i)
            ax01plot(x[z_i:], y1[z_i:], label=label_i)
            ax11plot(x[z_i:], y2[z_i:], label=label_i)

        self.__set_scale(ax, xscale, yscale)
        self.__set_label(ax, ['', '', xlabel, xlabel], [y1label, '', y2label,
                                                        ''])
        plt.legend(loc=0)
        plt.show()
        plt.close()

    def plot_fraction_density(self, varied_name, labelvaried_name, z_s=100,
                              yscale=['log', 'log'],
                              xscale=['log', 'log'], exclude=[], species=[]):
        """
        Plot Dark Energy fraction density. Second figure is plotted up to
        redshift z_s.

        varied_name = varied variable's name
        labelvaried_name = label for varied_name
        z_s = greatest z to plot in zoomed figures
        xscale = list of scales for the x axis
        yscale = list of scales for the y axis
        exclude = list of the varied variable values to exclude from plotting.
        species = list of extra species to plot, called as in hi_class tables.
        """
        xlabel = 'z+1'
        ylabel = r'$\Omega_{DE}$'

        f, ax = plt.subplots(1, 2, figsize=(15, 6))

        for i, ba in self.computed[varied_name].iteritems():
            if (i in exclude) or (not ba):
                continue

            z = ba['back']['z']
            z_i = wicmath.find_nearest(ba['back']['z'], z_s)

            OmegaDE = ba['back']['(.)rho_smg'] / ba['back']['(.)rho_crit']

            ax[0].plot(z + 1, OmegaDE, label=labelvaried_name + '={}'.format(i))
            ax[1].plot(z[z_i:] + 1, OmegaDE[z_i:], label=labelvaried_name + '={}'.format(i))

        for s in species:
            subindex_s = s.split('_')[-1]
            Omega_s = ba['back'][s] / ba['back']['(.)rho_crit']
            ax[0].plot(z + 1, Omega_s)
            ax[1].plot(z[z_i:] + 1, Omega_s[z_i:], label=r'$\Omega_{}$'.format(subindex_s))

        Omega0_Planck = np.ones(len(z)) * (1 - self.cosmo.Omega_m())
        ax[0].plot(z + 1, Omega0_Planck)
        ax[1].plot(z[z_i:] + 1, Omega0_Planck[z_i:], label=r'$1 - \Omega_{m,0}}$')

        self.__set_scale(ax, xscale, yscale)
        self.__set_label(ax, [xlabel, xlabel], [ylabel, ''])

        plt.legend(loc=0)
        plt.show()
        plt.close()

    def plot_density(self, varied_name, labelvaried_name, z_s=100,
                     yscale=['log', 'log', 'log', 'log'],
                     xscale=['log', 'log', 'log', 'log'],
                     exclude=[], species=[]):
        """
        Plot Dark Energy and critical density. Second figure is plotted until
        redshift z_s.

        varied_name = varied variable's name
        labelvaried_name = label for varied_name
        z_s = greatest z to plot in zoomed figures
        xscale = list of scales for the x axis
        yscale = list of scales for the y axis
        exclude = list of the varied variable values to exclude from plotting.
        species = list of extra species to plot, called as in hi_class tables.
        """
        xlabel = 'z+1'
        y1label = r'$\rho_c$'
        y2label = r'$\rho_{DE}$'

        f, ax = plt.subplots(2, 2, figsize=(15, 10))

        for i, ba in self.computed[varied_name].iteritems():
            if (i in exclude) or (not ba):
                continue

            z = ba['back']['z']
            z_i = wicmath.find_nearest(ba['back']['z'], z_s)

            rho = ba['back']['(.)rho_smg']
            rho_c = ba['back']['(.)rho_crit']

            ax[0, 0].plot(z + 1, rho_c, label=labelvaried_name + '={}'.format(i))
            ax[0, 1].plot(z[z_i:] + 1, rho_c[z_i:], label=labelvaried_name + '={}'.format(i))
            ax[1, 0].plot(z + 1, rho, label=labelvaried_name + '={}'.format(i))
            ax[1, 1].plot(z[z_i:] + 1, rho[z_i:], label=labelvaried_name + '={}'.format(i))

        rho0_Planck = np.ones(len(z)) * (1 - self.cosmo.Omega_m()) * rho_c[-1]
        ax[0, 0].plot(z + 1, rho0_Planck)
        ax[0, 1].plot(z[z_i:] + 1, rho0_Planck[z_i:])
        ax[1, 0].plot(z + 1, rho0_Planck)
        ax[1, 1].plot(z[z_i:] + 1, rho0_Planck[z_i:], label=r'$\rho_{DE}^{expected}$')

        for s in species:
            subindex_s = s.split('_')[-1]
            rho_s = ba['back'][s]
            ax[1, 0].plot(z + 1, rho_s)
            ax[1, 1].plot(z[z_i:] + 1, rho_s[z_i:], label=r'$\rho_{}$'.format(subindex_s))

        self.__set_scale(ax, xscale, yscale)
        self.__set_label(ax, ['', '', xlabel, xlabel], [y1label, '', y2label,
                                                        ''])

        ax[1, 0].set_xlabel('z+1')
        ax[1, 1].set_xlabel('z+1')
        ax[1, 0].set_ylabel(r'$\rho_{DE}$')
        ax[1, 1].set_ylabel(r'$\rho_{DE}$')
        ax[0, 0].set_ylabel(r'$\rho_c$')
        ax[0, 1].set_ylabel(r'$\rho_c$')
        plt.legend(loc=0)
        plt.show()
        plt.close()

    def plot(self, varied_name, x, y, labelvaried_name, xlabel, ylabel,
             scale=['linear', 'linear'], exclude=[], add=[0, 0], x_s=False, scatter=False):
        """
        Plot y vs x for all varied_names values.

        varied_name = varied variable's name
        x = list with ['variable name', 'module'] with  module = back, thermo...
            for X axis
        y = list with ['variable name', 'module'] with  module = back, thermo...
            for Y axis
        ylabelvaried_name = label for varied_name
        ylabel = label for Y1 axis
        scale = list of scale type for x and y axis e.g. ['linear', 'linear'].
        add = 2-item list with values to sum to x and y arrays respectively.
        x_s = greatest x to plot.
        exclude = list of the varied variable values to exclude from plotting.
        """
        fig, ax = plt.subplots(figsize=(10, 7))

        if scatter:
            axPlot = ax.scatter
        else:
            axPlot = ax.plot

        for i, ba in self.computed[varied_name].iteritems():
            if i in exclude:
                continue
            x1 = ba[x[1]][x[0]] + add[0]
            y1 = ba[y[1]][y[0]] + add[1]

            if x_s:
                x_i = wicmath.find_nearest(ba['back']['z'], x_s)

            label_i = labelvaried_name + '={}'.format(i)

            if x1[x_i] < x1[x_i + 1]:
                axPlot(x1[:x_i + 1], y1[:x_i + 1], label=label_i)
            else:
                axPlot(x1[x_i:], y1[x_i:], label=label_i)

            ax.set_xscale(scale[0])
            ax.set_yscale(scale[1])

            ax.set_ylabel(r'${}$'.format(ylabel))
            ax.set_xlabel(r'${}$'.format(xlabel))

        plt.legend(loc=0)
        plt.show()
        plt.close()
