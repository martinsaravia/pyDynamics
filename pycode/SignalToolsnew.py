# -*- coding: utf-8 -*-
"""
Created on Thu Apr 21 23:54:45 2016

@author: root
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.fftpack
import os
from mslib import msutil as mu


    
def fft( y, t):
    n = len(y) # length of the signal
    Fs= 1 / ( t[1] - t[0] )#sampling Rate
    k = np.arange(n)
    T = n/Fs
    frq = k/T # two sides frequency range
    frq = frq[range(n/2)] # one side frequency range
    Y = np.fft.fft(y)/n # fft computing and normalization
    Y = Y[range(n/2)]
      
    plt.plot(frq,abs(Y),'r') # plotting the spectrum

    


def fft2( y, t, lim, exp, name):

    # Number of samplepoints
    N = len(y)
    # sample spacing
    T = ( t[1] - t[0] )
    x = np.linspace(0.0, N*T, N)
    
    yf = scipy.fftpack.fft(y)
    xf = np.linspace(0.0, 1.0/(2.0*T), N/2)
    
    fig, ax = plt.subplots()
    ax.plot(xf, 2.0/N * np.abs(yf[:N/2]))
    plt.xlabel('Frequency (Hz)',   fontsize=14)
    plt.ylabel('Voltage (V)',fontsize=14)
    plt.xlim(lim)
    plt.show()
    if exp == 'yes':
        mu.figexp(name, plt)
    
def ReadNational(self, name):
    
    cpath = os.getcwd()
    self.jname = name
    self.inpdir = cpath + '\\' + 'data' + '\\'
    self.ifile = open(self.inpdir + self.jname + '.dat' )      
    print ( ('-> Running National Instruments Reader for input file:' + self.jname))  
    signal = np.loadtxt(self.ifile, dtype='float', skiprows=2)
    return signal

#        chan = 3
#        for i, l in enumerate(self.ifile):
#            if i > 1:
#                print (float(linea[1]))
#                slist.append( [float(linea[0]), float(linea[1]), float(linea[2]) ] )
#        #Convert list to array
#        signal = np.zeros((len(slist),chan))
#        for i, line in enumerate(slist):
#            for c in range(chan):
#                signal[i, c] =  line[i]
#        
#        return signal
         
       
def SavGol(y, window_size, order, deriv=0):
#   ## mooth (and optionally differentiate) data with a Savitzky-Golay filter.
#    The Savitzky-Golay filter removes high frequency noise from data.
#    It has the advantage of preserving the original shape and
#    features of the signal better than other types of filtering
#    approaches, such as moving averages techhniques.
#    
#    This code has been taken from http://www.scipy.org/Cookbook/SavitzkyGolay
#    Parameters
#    ----------
#    y : array_like, shape (N,)
#        the values of the time history of the signal.
#    window_size : int
#        the length of the window. Must be an odd integer number.
#    order : int
#        the order of the polynomial used in the filtering.
#        Must be less then `window_size` - 1.
#    deriv: int
#        the order of the derivative to compute (default = 0 means only smoothing)
#    Returns
#    -------
#    ys : ndarray, shape (N)
#        the smoothed signal (or it's n-th derivative).
#    Notes
#    -----
#    The Savitzky-Golay is a type of low-pass filter, particularly
#    suited for smoothing noisy data. The main idea behind this
#    approach is to make for each point a least-square fit with a
#    polynomial of high order over a odd-sized window centered at
#    the point.
#    Examples
#    --------
#    t = np.linspace(-4, 4, 500)
#    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
#    ysg = savitzky_golay(y, window_size=31, order=4)
#    import matplotlib.pyplot as plt
#    plt.plot(t, y, label='Noisy signal')
#    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
#    plt.plot(t, ysg, 'r', label='Filtered signal')
#    plt.legend()
#    plt.savefig('images/golay.png')
#    #plt.show()
#    References
#    ----------
#    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
#       Data by Simplified Least Squares Procedures. Analytical
#       Chemistry, 1964, 36 (8), pp 1627-1639.
#    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
#       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
#       Cambridge University Press ISBN-13: 9780521880688

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError:
        print('SavGol filter Error')
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv]
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    #    print('convolving')
    return np.convolve( m, y, mode='valid')
        
        
        