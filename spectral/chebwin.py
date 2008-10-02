import numpy as N
import scipy.signal as sig
import scipy.special as spec
import matplotlib.pyplot as P

def mycheb(num, atten):
    gammainv = 10.**(atten/20)
    M = (num - 1)/2
    beta = N.cosh(N.arccosh(gammainv)/(2*M))
    n = N.arange(-M, M+1).reshape(-1,1)
    k = N.arange(1, M+1)
#    T_2M = spec.chebyt(2*M)
    T_2M = myT(2*M, beta*N.cos(k*N.pi/(2*M+1)))
    sumterm = (T_2M*N.cos(2*n*k*N.pi/(2*M+1))).sum(axis=-1)
    w = (gammainv + 2 * sumterm) / (2 * M + 1)
    return w

def mycheb2(M, at, sym=1):
    """Dolph-Chebyshev window.

    INPUTS:

      M : int
        Window size
      at : float
        Attenuation (in dB)
      sym : bool
        Generates symmetric window if True.

    """
    if M < 1:
        return array([])
    if M == 1:
        return ones(1,'d')

    odd = M % 2
    if not sym and not odd:
        M = M+1

    # compute the parameter beta
    beta = N.cosh(1.0/(M-1.0)*N.arccosh(10**(abs(at)/20.)))
#    print beta
#    print M
    k = N.r_[0:M]*1.0
    x = beta*N.cos(N.pi*k/M)
#   print x
    # find the window's DFT coefficients
    # using the Chebyshev polynomial
#    print spec.chebyt(M - 1)
#    p = N.polyval(spec.chebyt(M - 1), x)
    p = myT(M - 1, x)
    #    print p, x
    # Appropriate IDFT and filling up
    # depending on even/odd M
    if M % 2:
        w = N.real(N.fft.fft(p * 1.0));
        n = (M + 1) / 2;
        w = w[:n] / w[0];
        w = N.concatenate((w[n - 1:0:-1], w))
    else:
        p = p * N.exp(1.j*N.pi / M * N.r_[0:M])
        w = N.real(N.fft.fft(p));
        n = M / 2 + 1;
        w = w / w[1];
        w = N.concatenate((w[n - 1:0:-1], w[1:n]));
    if not sym and not odd:
        w = w[:-1]
    return w

def myT(order, x):
    retval = N.zeros_like(x)
    retval[x > 1] = N.cosh(order*N.arccosh(x[x>1]))
    retval[x < -1] = N.cosh(order*N.arccosh(-x[x<-1]))*((-1)*(order%2))
    retval[N.abs(x)<=1] = N.cos(order*N.arccos(x[N.abs(x)<=1]))
    return retval
