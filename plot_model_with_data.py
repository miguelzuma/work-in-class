#!/usr/bin/python
# -*- coding: utf-8 -*-

# Docstring:
"""This module provides an interface to interact directly with MontePython
Likelihoods."""

import sys
import os
from collections import OrderedDict

class CommandLine():
    """CommandLine class needed at likelihood initialization."""
    def __init__(self, folder):
       self.folder = folder


class Data():
    """Data class parallel to MontePython's defined only with the required
    parameters to load likelihoods.

    Once Data is initialized, likelihoods classes are accesible by
    self.lkl['Likelihood_name'].
    """
    def __init__(self, path_MP, experiments, path_data='data'):
        self.path = {}
        self.path['root'] = os.path.abspath(path_MP)
        self.path['MontePython'] = self.path['root'] + '/montepython'
        self.path['data'] = os.path.join(self.path['root'], path_data)
        self.log_flag = True # This makes Likelihood to not write a log file
        self.experiments = experiments
        self.cosmo_arguments = {}
        """
        Simple dictionary that will serve as a communication interface with the
        cosmological code. It contains all the parameters for the code that
        will not be set to their default values.

        :rtype:   dict

        Documentation from montepython/data.py (Modified)
        """
        self.mcmc_parameters = OrderedDict()
        """
        Ordered dictionary of dictionaries, it will contain the 'nuisance'
        parameters which are the only used from this variable in the Likelihood
        class.  Every parameter name will be the key of a dictionary, containing
        the current value and its scale.

        Example: OrderedDict(['A_planck', {'current': value, 'scale': value2}])

        :rtype: ordereddict

        Documentation from montepython/data.py (Modified)
        """
        self.cosmological_module_name = 'CLASS'
        self.command_line = CommandLine('/tmp')
        self.__initialise_likelihoods()
        os.remove(os.path.join(self.command_line.folder, 'log.param'))

    def __initialise_likelihoods(self):
        """
        Given an array of experiments, return an ordered dict of instances

        Copied from montepython/data.py

        """

        self.lkl = OrderedDict()
        # adding the likelihood directory to the path, to import the module
        # then, for each library, calling an instance of the likelihood.
        # Beware, though, if you add new likelihoods, they should go to the
        # folder likelihoods/yourlike/yourlike.py, and contain a yourlike.data,
        # otherwise the following set of commands will not work anymore.

        # For the logging if log_flag is True, each likelihood will log its
        # parameters

        # Due to problems in relative import, this line must be there. Until a
        # better solution is found. It adds the root folder of the MontePython
        # used as the first element in the sys.path
        sys.path.insert(0, self.path['root'])
        sys.path.insert(1, self.path['MontePython'])

        for elem in self.experiments:

            folder = os.path.abspath(os.path.join(
                self.path['MontePython'], "likelihoods", "%s" % elem))
            # add the folder of the likelihood to the path of libraries to...
            # ... import easily the likelihood.py program
            try:
                exec "from likelihoods.%s import %s" % (
                    elem, elem)
            except ImportError as message:
                raise io_mp.ConfigurationError(
                    "Trying to import the %s likelihood" % elem +
                    " as asked in the parameter file, and failed."
                    " Please make sure it is in the `montepython/"
                    "likelihoods` folder, and is a proper python "
                    "module. Check also that the name of the class"
                    " defined in the __init__.py matches the name "
                    "of the folder. In case this is not enough, "
                    "here is the original message: %s\n" % message)
            # Initialize the likelihoods. Depending on the values of
            # command_line and log_flag, the routine will call slightly
            # different things. If log_flag is True, the log.param will be
            # appended.
            try:
                exec "self.lkl['%s'] = %s('%s/%s.data',\
                    self, self.command_line)" % (
                    elem, elem, folder, elem)
            except KeyError as e:
                if e.find('clik') != -1:
                    raise io_mp.ConfigurationError(
                        "You should provide a 'clik' entry in the dictionary "
                        "path defined in the file default.conf")
                else:
                    raise io_mp.ConfigurationError(
                        "The following key: '%s' was not found" % e)

    def get_mcmc_parameters(self, table_of_strings):
        """
        Returns an ordered array of parameter names filtered by
        `table_of_strings`.

        Parameters
        ----------
        table_of_strings : list
            List of strings whose role and status must be matched by a
            parameter. For instance,

            >>> data.get_mcmc_parameters(['varying'])
            ['omega_b', 'h', 'amplitude', 'other']

            will return a list of all the varying parameters, both
            cosmological and nuisance ones (derived parameters being `fixed`,
            they wont be part of this list). Instead,

            >>> data.get_mcmc_parameters(['nuisance', 'varying'])
            ['amplitude', 'other']

            will only return the nuisance parameters that are being varied.

        Copied from montepython/data.py

        """
        table = []
        for key, value in self.mcmc_parameters.iteritems():
            number = 0
            for subvalue in value.itervalues():
                for string in table_of_strings:
                    if subvalue == string:
                        number += 1
            if number == len(table_of_strings):
                table.append(key)
        return table

    def update_params(self, params):
        """
        Given a classy-type parameters dictionary update the
        data.cosmo_arguments dictionary.

        Paramters_smg can be set as in classy or as in MontePython.
        """

        self.cosmo_arguments.update(params)

        if 'parameters_smg' in params:
            parameters_smg = self.cosmo_arguments['parameters_smg'].split(',')

            for i, value in enumerate(parameters_smg):
                self.cosmo_arguments['parameters_smg__{}'.format(i+1)] = value

            del self.cosmo_arguments['parameters_smg']

    def update_nuisance(self, nuisance):
        """Update the nuisance parameters in the mcmc_parameters dictionary.
        Input a dictionary or list of list containing the nuisance parameter
        name and its value. All scales are put to 1."""

        if type(nuisance) is dict:
            nuisance = nuisance.iteritems()

        for key, value in nuisance:
            self.mcmc_parameters[key] = {'current': value, 'scale': 1, 'role':
                                          'nuisance'}

    def set_params(self, params):
        """
        Given a classy parameters dictionary fill the data.cosmo_arguments
        variable.

        Paramters_smg can be set as in classy or as in MontePython.

        If cosmo = Class(), you can use cosmo.pars as input parameter.
        """

        if self.cosmo_arguments:
            self.cosmo_arguments.clear()

        self.cosmo_arguments = params.copy()

        if 'parameters_smg' in params:
            parameters_smg = self.cosmo_arguments['parameters_smg'].split(',')

            for i, value in enumerate(parameters_smg):
                self.cosmo_arguments['parameters_smg__{}'.format(i+1)] = value

            del self.cosmo_arguments['parameters_smg']

    def set_nuisance(self, nuisance):
        """Set the nuisance parameters in the mcmc_parameters dictionary. Input a
        dictionary or list of list containing the nuisance parameter name and
        its value. All scales are put to 1."""

        if self.mcmc_parameters:
            self.mcmc_parameters.clear()

        if type(nuisance) is dict:
            nuisance = nuisance.iteritems()

        for key, value in nuisance:
            self.mcmc_parameters[key] = {'current': value, 'scale': 1, 'role':
                                         'nuisance'}

