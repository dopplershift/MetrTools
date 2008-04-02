import numpy as N

mp_N0 = 8.0e3 # m^-3 mm^-1
density_water = 1.0e3 # Kg m^-3

def mp_slope_3rd(lwc):
    '''Gives the slope for the Marshall-Palmer drop size distribution
       given the liquid water content mixing ratio in g/kg. Returns
       slope in mm^-1.'''
    return (N.pi * density_water * mp_N0 * 1.0e-6 / lwc)**0.25

def mp_from_lwc(d, lwc):
    return marshall_palmer(d, mp_slope_3rd(lwc))

def marshall_palmer(d, lam):
    return exponential(d, lam, mp_N0)

def exponential(d, lam, N0):
    return gamma(d, lam, N0, 0.0)

def gamma(d, lam, N0, mu):
    return N0 * d**mu * N.exp(-lam * d)

if __name__ == '__main__':
    import pylab as P
    lwc = 19.0
    d = N.linspace(0.01, 10.0, 100)
    mp = mp_from_lwc(d, lwc)
    P.semilogy(d, mp, 'b')
    P.xlabel('Diameter (mm)')
    P.ylabel(r'Concentration (# m$^{-3}$)')
    P.title(r'Marshall-Palmer Distribution for %.1f g kg$^{-1}$ LWC' % lwc)
    P.show()
