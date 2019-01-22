#!/usr/bin/python
import numpy as np
from classy import CosmoSevereError
import scipy.integrate as integrate

#####
# NOTE: Emilio's functions from classy.pyx. self should not be used like this, but let's do an exception.
#####

def growthrate_at_k_and_z(cosmo, k, z):
    dz = 0.005
    f=-0.5*(1+z)*(-3.*np.log(cosmo.pk(k,z))+4.*np.log(cosmo.pk(k,z+dz))-np.log(cosmo.pk(k,z+2.*dz)))/(2.*dz)
    if (f==0.):
        raise CosmoSevereError(
            "The step in redshift is too small to compute growth rate")
    return f

def growthrate_at_z(cosmo, z):
    k_fid = 0.01
    return growthrate_at_k_and_z(cosmo, k_fid, z)
#########

def fprime(OmegaDEw, OmegaM):
    """
    Return function df/dz(z, f).

    OmegaDEw and OmegaM must be functions of z
    """
    def wrap(z, f):
        try:
            #output = 1.5 * OmegaM(z) - 0.5 * (1 - 3*w(z)*OmegaDE(z)) * f - f**2
            output = ((0.5 * (1 - 3*OmegaDEw(z))) * f + f**2 - 1.5 * OmegaM(z)) / (1+z)
        except Exception as e:
            raise e
        return output
    return wrap
def angular_distance(z, H, z_points=None):
    """
    Return the angular distance.
    z, H must be CLASS-like arrays
    """
    DA_fit = []
    if z_points is None:
        z_points = z

    for i in z_points:
        #DA_fit.append(1/(1+i)*integrate.trapz(1/H_fit[zTMP<=i][::-1], zTMP[zTMP<=i][::-1]))
        DA_fit.append(1/(1+i)*integrate.simps(1/H[z <= i][::-1], z[z <= i][::-1], even='last'))

    return np.array(DA_fit)
