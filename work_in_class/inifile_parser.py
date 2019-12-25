#!/usr/bin/python


# Docstring
"""Set of functions helpful to deal with ini files."""

def inifile_parser(filename):
    """Return a dictionary with key the name of the parameter
    and value, its value."""
    params = {}
    with open(filename) as f:
        for line in f:
            if not ('#' in line) and ('=' in line):
                var, key = line.split("=")
                params[var.strip()] = key.strip()
    return params


def parameters_smg(parameters_str):
    """Return a list of floats with the parameters of parameters_smg"""
    return map(float, parameters_str.split(", "))


def vary_params(parameters_str, new):
    """
    Return a modified parameters_smg string with the values given in the 'new' array.
    params = old params to modify

    new = [[index_of_value1, value1], [index_of_value2, value2],...]
    """

    params = parameters_str.split(',')
    for i in new:
        params[i[0]] = str(i[1])
    return ','.join(params)
