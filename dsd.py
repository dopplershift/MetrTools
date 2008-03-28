import numpy as N

mp_N0 = 8.0e3
density_water = 1.0e3

def mp_slope_3rd(lwc):
    '''Gives the slope for the Marshall-Palmer drop size distribution
       given the liquid water content mixing ratio in g/kg'''
    return (N.pi * density_water * mp_N0 * 1.0e6 / lwc)**0.25

def mp_from_lwc(d, lwc):
    return marshall_palmer(d, mp_slope_3rd(lwc))

def marshall_palmer(d, lam):
    return exponential(d, lam, mp_N0)

def exponential(d, lam, N0):
    return gamma(d, lam, N0, 0.0)

def gamma(d, lam, N0, mu):
    return (N0 * (d[N.newaxis,:]**mu) * N.exp(-lam * d)).squeeze()

if __name__ == '__main__':
    import pylab as P
    lwc = 19.0
    d = N.linspace(0.01, 1.0, 100)
    mp = mp_from_lwc(d, lwc)
    P.semilogy(d, mp)
    P.xlabel('Diameter (cm)')
    P.ylabel('Concentration (#/m^3)')
    P.title('Marshall-Palmer Distribution for %.2f g/kg LWC' % lwc)
    P.show()
