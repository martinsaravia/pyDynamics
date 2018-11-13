# -*- coding: utf-8 -*-
"""
Created on Tue Jun 06 18:12:01 2017

@author: Root
"""

#import os; clear = lambda: os.system('cls'); clear() 
import matplotlib.pyplot as plt, numpy as np
import matplotlib as mpl
from scipy.optimize import curve_fit
from scipy import interpolate
import time

t0 = time.time()
mpl.use('Qt4Agg') 
font = {'family' : 'Times New Roman','weight' : 'normal','size'   : 12}
mpl.rc('font', **font)


unitconv = 1.0 / 1000.0
gr2N = 1.0 / 100.0



# Experimento de la curva para caida libre con altura 32.2 mm (medido con buena alineacion, usando el timon de la fresadora)
xp038_experimental_damp_iden = unitconv *  np.array(  [ -0.5, 1.0, 2.0, 5.0, 14, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0, 31.0, 40.0   ])   # 12.0 (  0.584566 ) 
Fx038_experimental_damp_iden = gr2N *      np.array(  [ 0.0, 0.0, 0.0, 0.0,  0.0, 0.1,  0.2,  0.33, 0.5,   0.8,  1.2,  1.7,  2.4,  3.4,  4.9,  6.9,  9.9, 14.1, 20.7, 30.4, 44.9, 400.0  ]) 


xp038_experimental_damp_iden_2 = unitconv * (-9.8 + np.array(  [  0.0, 5.0,  10.0, 28.0, 29.0, 30.0, 31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0, 40.0, 41.0, 42.0 ])   )# 12.0 (  0.584566 ) 
Fx038_experimental_damp_iden_2 = gr2N *             np.array(  [  0.0,  0.0,  0.0,  0.1,  0.2,  0.4,  0.7,  1.0,  1.5,  2.1,  3.0,  4.2,  5.9,  8.4, 12.1, 17.5, 25.4, 37.5 ]) 

xp = xp038_experimental_damp_iden_2
Fx = Fx038_experimental_damp_iden_2

xint = np.arange(-0.01, unitconv*35.0, unitconv*0.1)
xint = np.flipud(-xint)
xp = np.flipud(-xp)
Fx = np.flipud(-Fx)



MM = np.zeros((2, len(xp)))
MM[0,:] = xp
MM[1,:] = Fx
np.savetxt('workfile.dat',  MM.T, fmt='%1.3e', delimiter=',')   # X is an array



# ---------   Spline interpiolation of experimental curve   ---------------

Fint = interpolate.spline(xp, Fx, xint, order=3, kind='smoothest', conds=None)
tck = interpolate.splrep(xp, Fx)
Fint2 = interpolate.splev(xint, tck, der=0, ext=0)
Fint2_der = interpolate.splev(xint, tck, der=1, ext=1)

finter= interpolate.spline(xp, Fx, 0., order=3, kind='smoothest', conds=None)

print(interpolate.splev(0.001, tck, der=1, ext=1))

# ---------   Exponential fitting curve   ---------------

def fitFunc(x, a, b, c):
    return a*np.exp(-b*x) + c

cfp = curve_fit(fitFunc, xp, Fx)
Fint_exp = fitFunc(xint, cfp[0][0], cfp[0][1], cfp[0][2])

plt.figure(1)
p1 = plt.plot(xp, Fx, 'b*', label='Force - Experiment')
#p2 = plt.plot(xint, Fint, 'k', label='Force - Spline Fit')
p3 = plt.plot(xint, Fint2, 'r', label='Force - Spline Fit from Rep')
#p4 = plt.plot(xint, Fint_exp , 'g-', label='Force - Exponential Fit from Rep')
plt.legend(handles=[p1[0], p2[0], p3[0], p4[0] ]) 

plt.figure(2)
#p5 = plt.plot(xint, Fint2 / xint, 'r',    label='Stiffness from Splines Fit')
#p6 = plt.plot(xint, Fint_exp / xint, 'k', label='Stiffness from Exponential Fit')
tol = 1.0E-5
p7 = plt.plot(xint, (Fint2) / (xint+tol), 'g', label='Stiffness from Splines Fit with tol')
#p8 = plt.plot(xint, (Fint_exp) / (xint+tol), 'y', label='Stiffness from Exponential Fit with tol')
plt.legend(handles=[ p5[0], p6[0], p7[0], p8[0]]) 
plt.show()
# ----------------------------------------------------------------------------




