#!/usr/bin/pyhon

"""Module thought to be used when you want to obtain some quantities (e.g. the
background structures) for a model varying one of the parameters"""

from classy import Class, CosmoSevereError
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

    def __store_cl(self, cl_dic):
        """
        Store cl's as (l*(l+1)/2pi)*cl, which is much more useful.
        """

        ell = cl_dic['ell'][2:]

        for cl, list_val in cl_dic.iteritems():
            list_val = list_val[2:]
            if (list_val == ell).all():
                cl_dic[cl] = list_val
                continue
            list_val = (ell * (ell + 1) / (2 * np.pi)) * list_val
            cl_dic[cl] = list_val  # Remove first two null items (l=0,1)

        return cl_dic

    def compute_models(self, params, varied_name, index_variable, values,
                       back=[], thermo=[], prim=[], pert=[], trans=[],
                       pk=[0.0001, 0.1, 100], extra=[], update=True,
                       cosmo_msg=False):
        """
        Fill dic with the hi_class output structures for the model with given
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
        pk = list with the minimum and maximum k values to store the present
             matter power spectrum and the number of points [k_min, k_max,
             number_points]. Default [10^-4, 10^1, 100].
        extra = list of any of the method or objects defined in cosmo, e.g.
                w0_smg().  It will store {'method': cosmo.w0_smg()}
        update = if True update old computed[key] dictionary elsewise create a
                 new one.  Default: True.
        cosmo_msg = if True, print cosmo.compute() messages. Default: False.
        """

        key = varied_name

        if (not update) or (key not in self.computed.keys()):
            self.computed[key] = od()

        for val in values:
            # key = "{}={}".format(varied_name, val)
            params["parameters_smg"] = inip.vary_params(params["parameters_smg"], [[index_variable, val]])

            # It might be after the try to not store empty dictionaries.
            # Nevertheless, I find more useful having them to keep track of
            # those failed and, perhaps, to implement a method to obtain them
            # with Omega_smg_debug.
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

            try:
                d['cl'] = self.__store_cl(self.cosmo.raw_cl())
            except CosmoSevereError:
                pass

            try:
                d['lcl'] = self.__store_cl(self.cosmo.lensed_cl())
            except CosmoSevereError:
                pass

            try:
                d['dcl'] = self.cosmo.density_cl()
            except CosmoSevereError:
                pass


            if ("output" in self.cosmo.pars) and ('mPk' in self.cosmo.pars['output']):
                k_array = np.linspace(*pk)
                pk_array = []
                for scale in k_array:
                    pk_array.append(self.cosmo.pk(scale, 0))

                d['pk'] = {'k': k_array, 'pk': np.array(pk_array)}

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

    def plot_4_density(self, varied_name, labelvaried_name, z_s=100,
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

            ba_s = ba

        rho0_Planck = np.ones(len(z)) * (1 - self.cosmo.Omega_m()) * rho_c[-1]
        ax[0, 0].plot(z + 1, rho0_Planck)
        ax[0, 1].plot(z[z_i:] + 1, rho0_Planck[z_i:])
        ax[1, 0].plot(z + 1, rho0_Planck)
        ax[1, 1].plot(z[z_i:] + 1, rho0_Planck[z_i:], label=r'$\rho_{DE}^{expected}$')

        for s in species:
            subindex_s = s.split('_')[-1]
            rho_s = ba_s['back'][s]
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

    def plot_cl(self, varied_name, labelvaried_name, cl=['tt', 'cl'], exclude=[], scale=['log', 'linear'], xlim=[]):
        """
        Plot angular power spectra: CMB raw or lensed power spectra.

        varied_name = varied variable's name
        labelvaried_name = label for varied_name
        cl = 2-item list whose first item is the desired spectra (e.g. 'tt',
             'ee'...) and whose second is 'cl' or 'lcl' for raw or lensed
             cl's. Default ['tt', 'cl'].
        exclude = list of the varied variable values to exclude from plotting.
        scale = list with scale for x and y axis. Default is ['log', 'linear']
        xlim = x limits to plot [x_min, x_max]. Default []
        """

        cl_label = r'$ [l (l+1) / 2\pi] C_l^{{{}}}$'.format(cl[0].upper())

        self.plot(varied_name, ['ell', 'cl'], cl, labelvaried_name, 'l', cl_label, exclude=exclude, scale=scale, xlim=xlim)

    def plot_pk(self, varied_name, labelvaried_name, exclude=[], scale=['log', 'log'], xlim=[]):
        """
        Plot present matter power spectrum.

        varied_name = varied variable's name
        labelvaried_name = label for varied_name
        exclude = list of the varied variable values to exclude from plotting.
        scale = list with scale for x and y axis. Default is ['log', 'log']
        xlim = x limits to plot [x_min, x_max]. Default []
        """

        self.plot(varied_name, ['k', 'pk'], ['pk', 'pk'], labelvaried_name, 'k', r'$P(k, z=0)$', exclude=exclude, scale=scale, xlim=xlim)

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

            ba_s = ba

        for s in species:
            subindex_s = s.split('_')[-1]
            Omega_s = ba_s['back'][s] / ba_s['back']['(.)rho_crit']
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
                     yscale=['log', 'log'], xscale=['log', 'log'], exclude=[],
                     species=[]):
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
        ylabel = r'$\rho_{DE}$'

        f, ax = plt.subplots(1, 2, figsize=(15, 6))

        for i, ba in self.computed[varied_name].iteritems():
            if (i in exclude) or (not ba):
                continue

            z = ba['back']['z']
            z_i = wicmath.find_nearest(ba['back']['z'], z_s)

            rho = ba['back']['(.)rho_smg']
            rho_c = ba['back']['(.)rho_crit']

            ax[0].plot(z + 1, rho, label=labelvaried_name + '={}'.format(i))
            ax[1].plot(z[z_i:] + 1, rho[z_i:], label=labelvaried_name + '={}'.format(i))

            ba_s = ba

        for s in species:
            subindex_s = s.split('_')[-1]
            rho_s = ba_s['back'][s]
            ax[0].plot(z + 1, rho_s)
            ax[1].plot(z[z_i:] + 1, rho_s[z_i:], label=r'$\rho_{}$'.format(subindex_s))

        rho0_Planck = np.ones(len(z)) * (1 - self.cosmo.Omega_m()) * rho_c[-1]
        ax[0].plot(z + 1, rho0_Planck)
        ax[1].plot(z[z_i:] + 1, rho0_Planck[z_i:])


        self.__set_scale(ax, xscale, yscale)
        self.__set_label(ax, [xlabel, xlabel], [ylabel, ''])

        plt.legend(loc=0)
        plt.show()
        plt.close()

    def plot(self, varied_name, x, y, labelvaried_name, xlabel, ylabel,
             scale=['linear', 'linear'], exclude=[], add=[0, 0], xlim=[],
             scatter=False):

        """
        Plot y vs x for all varied_names values.

        varied_name = varied variable's name
        x = list with ['variable name', 'module'] with  module = back, thermo...
            for X axis
        y = list with ['variable name', 'module'] with  module = back, thermo...
            for Y axis
        labelvaried_name = label for varied_name
        xlabel = label for x axis
        ylabel = label for y axis
        scale = list of scale type for x and y axis e.g. ['linear', 'linear'].
        add = 2-item list with values to sum to x and y arrays respectively.
        xlim = x limits to plot [x_min, x_max]. Default []
        exclude = list of the varied variable values to exclude from plotting.
        scatter = If True, plot scatter points.
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

            if xlim:
                indexes = []
                for index, x1_elem in enumerate(x1):
                    if (x1_elem > xlim[1]) or (x1_elem < xlim[0]):
                        indexes.append(index)

                x1 = np.delete(x1, indexes)
                y1 = np.delete(y1, indexes)

            label_i = labelvaried_name + '={}'.format(i)

            axPlot(x1, y1, label=label_i)

            ax.set_xscale(scale[0])
            ax.set_yscale(scale[1])

            ax.set_ylabel(ylabel)
            ax.set_xlabel(xlabel)

        plt.legend(loc=0)
        plt.show()
        plt.close()
