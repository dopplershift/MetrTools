'''
A collection of utility functions for doing drop-size distribution
calculations.
'''

import numpy as np
from constants import density_water
from scipy.constants import milli
from scipy.special import gamma as gamma_func

__version__ = 0.7

mp_N0 = 8.0e3 / milli # m^-3 mm^-1 / m mm^-1 -> m^-4

def constrained_gamma_shape(lam):
    '''Calculates the shape factor (mu) for the constrained gamma relation
       as given by Zhang et al. (2001), using the slope factor *lam*.'''
    return -0.016 * lam**2 + 1.213 * lam - 1.957

def mp_slope_3rd(lwc):
    '''Gives the slope for the Marshall-Palmer drop size distribution
       given the liquid water content density in kg/m^3. Returns
       slope in m^-1.'''
    return (np.pi * density_water * mp_N0 / lwc)**0.25

def mp_from_lwc(d, lwc):
    '''Returns the Marshall-Palmer distribution weights corresponding to the
       given diameters using the given liquid water content to calculate the
       slope parameter.  All quantities should be in MKS.'''
    return marshall_palmer(d, mp_slope_3rd(lwc))

def marshall_palmer(d, lam):
    '''Returns the Marshall-Palmer distribution weights corresponding to the
       given diameters using the given slope parameter.  All quantities should
       be in MKS.'''
    return exponential(d, lam, mp_N0)

def exponential(d, lam, N0):
    '''Returns the exponential distribution weights corresponding to the
       given diameters using the given slope and intercept parameters.  All
       quantities should be in MKS.'''
    return modified_gamma(d, lam, N0, 0.0)

def modified_gamma(d, lam, N0, mu):
    '''Returns the modifed gamma distribution weights corresponding to the
       given diameters using the given slope, intercept, and shape parameters.
       All quantities should be in MKS.'''
    return N0 * d**mu * np.exp(-d * lam)

def volume_gamma(d, d0, N, nu=-0.8):
    '''Returns the gamma distribution weights corresponding to the
       given diameters using the given parameters.
       All quantities should be in MKS.'''
    rat = d / d0
    return ((N * (3 *(nu + 1) ** (nu + 1)) / (d0 * gamma_func(nu + 1)))
        * (rat**(3*nu + 2) * np.exp(-(nu+1) * rat**3)))

def gamma_d0_from_lwc(lwc, N):
    '''Returns the gamma distribution d0 parameter given the value of liquid
       water content and number concentration.
       All quantities should be in MKS.'''
    return np.power(6 * lwc / (np.pi * density_water * N), 1. / 3.)

def rain_fallspeed(d):
    '''Returns the raindrop fallspeed in m s^-1 given the drop diameter in m'''
    d = d / milli # Formulas need diameter in mm.  Assumed passed in as m.
    #From Brandes et al. 2002
    #This formula does not work above 1.2 cm
    #vt = -0.1021 + 4.932*d - 0.9551*d**2 + 0.07934*d**3 - 0.002362*d**4
    vt = -0.1021 + d*(4.932 + d*(-0.9551 + d*(0.07934 - 0.002362*d)))
    return vt

def lwc(d, dsd):
    '''Returns the liquid water content in kg/m^3 given the number density
    and diameters. These should be in MKS.'''
    return np.trapz(d**3 * (np.pi * density_water / 6.) * dsd, x=d, axis=0)

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from scipy.constants import kilo as g_per_kg

    lwc = 19.0 / g_per_kg
    d = np.linspace(0.01, 10.0, 100) * milli # Up to 10mm in meters
    mp = mp_from_lwc(d, lwc)
    plt.semilogy(d / milli, mp * milli, 'b')
    plt.xlabel('Diameter (mm)')
    plt.ylabel(r'Concentration (# m$^{-3}$ mm$^{-1}$)')
    plt.title(r'Marshall-Palmer Distribution for %.1f g kg$^{-1}$ LWC'
        % (lwc * g_per_kg))
    plt.show()
