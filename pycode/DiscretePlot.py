# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 21:48:24 2015

@author: root
"""
import os; clear = lambda: os.system('cls'); clear() 
import matplotlib.pyplot as plt
import numpy as np


def PlotSpectrum(signal):
    
    sr = 1.0 / (signal[1] - signal[0]) # Sampling Rate
    ts = len(signal) # length of the signal
    k = np.arange(ts)
    T = ts / sr
    frq = k / T # two sides frequency range
    frq = frq[range(ts/2)] # one side frequency range
    
    Y = np.fft.fft(y)/n # fft computing and normalization
    Y = Y[range(n/2)]
     
     
    plt.plot(frq, abs(Y) )
 