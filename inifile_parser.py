#!/usr/bin/python

def inifile_parser(filename):
    """Return a dictionary with key the name of the parameter
    and value, its value."""
    params = {}
    with open(filename) as f:
        for line in f:
            if '=' in line:
                var, key = line.strip().split(" = ")
                params[var] = key
    return params


def parameters_smg(parameters_str):
    """Return a list of floats with the parameters of parameters_smg"""
    return map(float, parameters_str.split(", "))
