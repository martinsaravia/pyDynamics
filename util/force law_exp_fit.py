# -*- coding: utf-8 -*-
"""
Created on Tue Jun 06 18:12:01 2017

@author: Root
"""

#import os; clear = lambda: os.system('cls'); clear() 
import matplotlib.pyplot as plt, numpy as np
import matplotlib as mpl
from scipy.optimize import curve_fit
import time

t0 = time.time()
mpl.use('Qt4Agg') 
font = {'family' : 'Times New Roman','weight' : 'normal','size'   : 12}
mpl.rc('font', **font)


unitconv = 1.0 / 1000.0
gr2N = 1.0 / 100.0

def func(x, a, b, c):
   return a * np.exp(-b * x) + c

# Experimento de la curva para caida libre con altura 32.2 mm (medido con buena alineacion, usando el timon de la fresadora)
xp038_experimental_damp_iden = unitconv *  np.array(  [ 0.0, 5.0, 7.5, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0, 31.0    ])   # 12.0 (  0.584566 ) 
Fx038_experimental_damp_iden = gr2N *      np.array(  [ 0.0, 0.03, 0.05, 0.1,  0.2,  0.33, 0.5,   0.8,  1.2,  1.7,  2.4,  3.4,  4.9,  6.9,  9.9, 14.1, 20.7, 30.4, 44.9   ]) 

xp = xp038_experimental_damp_iden
Fx = Fx038_experimental_damp_iden

popt, pcov = curve_fit(func, xp, Fx)
p1 = plt.plot(xp, Fx, 'b-', label='data')
p2 = plt.plot(xp, func(xp, *popt), 'r-', label='fit')
plt.legend(handles=[p1, p2]) 



# Plots

#p1 = plt.plot(xp, Fx,      'r*', label='Finite Element')
##plt.plot(xp, Fx_3reg, 'b--', label='Cubic reg')
#
#p2 = plt.plot(xp0, Fx_3reg, 'b', label='Cubic reg')
#
#p3 = plt.plot(xp0, Fx_5reg, 'k', label='Quintic reg')
#
#p4 = plt.plot(xp0, Fx_7reg, 'r', label='Sept reg')
#
#
#
#plt.legend(handles=[p1, p2, p3, p4]) 
#
#plt.show()