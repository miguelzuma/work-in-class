#!/usr/bin/pyhon

"""Module thought to be used when you want to obtain some quantities (e.g. the
background structures) for a model varying one of the parameters"""

from collections import OrderedDict as od
from matplotlib import pyplot as plt
import inifile_parser as inip
import wicmath


class Model():
    def __init__(self, cosmo):
        self.cosmo = cosmo
        self.computed = {}

    def compute_models(self, params, varied_name, index_variable, values,
                       back=[], thermo=[], prim=[], pert=[], trans=[],
                       extra=[]):
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
        self.computed[key] = od()

        for val in values:
            # key = "{}={}".format(varied_name, val)
            params["parameters_smg"] = inip.vary_params(params["parameters_smg"], [[index_variable, val]])

            d = self.computed[key][val] = {}

            self.cosmo.set(params)

            try:
                self.cosmo.compute()
            except:
                print "Error: skipping {}={}".format(key, val)
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
                    scale=[['linear', 'linear'], ['linear', 'linear']],
                    exclude=[]):
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
        scale = list of lists of scale type for x and y axis e.g. ['linear',
                'linear']. First item is for upper plots, and second for the
                lowers
        exclude = list of the varied variable values to exclude from plotting.
        """
        fig, ax = plt.subplots(2, 2, figsize=(15, 10))
        xlabel = 'z+1'

        module = 'back'

        for i, ba in self.computed[varied_name].iteritems():
            if i in exclude:
                continue
            x = ba[module]['z'] + 1
            y1 = ba[module][y1name]
            y2 = ba[module][y2name]
            z_i = wicmath.find_nearest(ba[module]['z'], z_s)

            label_i = r'${}={}$'.format(labelvaried_name, i)

            ax[0, 0].semilogx(x, y1, label=label_i)

            ax[0, 0].set_ylabel(r'${}$'.format(y1label))
            ax[0, 0].set_yscale(scale[0][0])

            ax[1, 0].semilogx(x, y2, label=label_i)
            ax[1, 0].set_ylabel(r'${}$'.format(y2label))
            ax[1, 0].set_xlabel(r'${}$'.format(xlabel))
            ax[1, 0].set_yscale(scale[1][0])

            ax[0, 1].semilogx(x[z_i:], y1[z_i:], label=label_i)
            ax[0, 1].set_yscale(scale[0][1])

            ax[1, 1].semilogx(x[z_i:], y2[z_i:], label=label_i)
            ax[1, 1].set_xlabel(r'${}$'.format(xlabel))
            ax[1, 1].set_yscale(scale[1][1])

        plt.legend(loc=0)
        plt.show()
        plt.close()

    def plot(self, varied_name, x, y, labelvaried_name, xlabel, ylabel,
             scale=['linear', 'linear'], exclude=[]):
        """
        Plot a 2x2 figure, with two variables y1 and y2 vs z+1. The ones in
        the second column are a zoom of the ones in the first one.

        varied_name = varied variable's name
        x = list with ['variable name', 'module'] with  module = back, thermo...
            for X axis
        y = list with ['variable name', 'module'] with  module = back, thermo...
            for Y axis
        ylabelvaried_name = label for varied_name
        ylabel = label for Y1 axis
        scale = list of scale type for x and y axis e.g. ['linear', 'linear'].
        exclude = list of the varied variable values to exclude from plotting.
        """
        fig, ax = plt.subplots(figsize=(10, 7))

        for i, ba in self.computed[varied_name].iteritems():
            if i in exclude:
                continue
            x1 = ba[x[1]][x[0]]
            y1 = ba[y[1]][y[0]]

            label_i = r'${}={}$'.format(labelvaried_name, i)

            ax.plot(x1, y1, label=label_i)

            ax.set_xscale(scale[0])
            ax.set_yscale(scale[1])

            ax.set_ylabel(r'${}$'.format(ylabel))
            ax.set_xlabel(r'${}$'.format(xlabel))

        plt.legend(loc=0)
        plt.show()
        plt.close()
