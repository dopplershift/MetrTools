import numpy as N
from numpy import fft, fftshift
from matlab import xcorr

def periodogram(time_series, taper=None, L=None):
    'Implements the periodogram spectral estimator.'
    #If we don't have a taper, just use a rectangular one (array of ones)
    if not taper:
        taper = N.ones_like(time_series)
    #Other wise, make sure the sizes are compatible
    elif taper.size != time_series.size and taper.size != time_series.shape[-1]:
        raise ValueError('Data and Taper must have the same size')
        
    #Zero pad if necessary
    if L is None:
        L = time_series.shape[-1]
    elif L > time_series.shape[-1]:
        num_zeros = L - time_series.shape[-1]
        zeropad = N.zeros(time_series.shape[:-1] + (num_zeros,),
            dtype=time_series.dtype)
        time_series = N.concatenate((time_series, zero_pad), axis=-1)
    elif L < time_series.shape[-1]:
        time_series = time_series[...,:L]

    Z = fftshift(fft(taper * time_series))
    return N.abs(Z)*N.abs(Z)/time_series.size

def correlogram(time_series, scale='biased'):
    '''Implements the correlogram spectral estimator. This is often known as the
    Blackman-Tukey estimator with a rectangular window.'''
    R = xcorr(time_series, scale=scale)
    return N.abs(fftshift(fft(R)))

