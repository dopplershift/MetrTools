import numpy as N
from constants import density_water, mm_per_m

mp_N0 = 8.0e3 * mm_per_m # m^-3 mm^-1 * mm/m -> m^-4

def mp_slope_3rd(lwc):
    '''Gives the slope for the Marshall-Palmer drop size distribution
       given the liquid water content density in kg/m^3. Returns
       slope in m^-1.'''
    return (N.pi * density_water * mp_N0 / lwc)**0.25

def mp_from_lwc(d, lwc):
    return marshall_palmer(d, mp_slope_3rd(lwc))

def marshall_palmer(d, lam):
    return exponential(d, lam, mp_N0)

def exponential(d, lam, N0):
    return gamma(d, lam, N0, 0.0)

def gamma(d, lam, N0, mu):
    return N0 * d**mu * N.exp(-lam * d)

if __name__ == '__main__':
    import matplotlib.pyplot as P
    from constants import g_per_kg
    lwc = 19.0 / g_per_kg 
    d = N.linspace(0.01, 10.0, 100) / mm_per_m # Up to 10mm in meters
    mp = mp_from_lwc(d, lwc)
    P.semilogy(d * mm_per_m, mp / mm_per_m, 'b')
    P.xlabel('Diameter (mm)')
    P.ylabel(r'Concentration (# m$^{-3}$ mm$^{-1}$)')
    P.title(r'Marshall-Palmer Distribution for %.1f g kg$^{-1}$ LWC'
        % (lwc * g_per_kg))
    P.show()
