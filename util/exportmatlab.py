# -*- coding: utf-8 -*-
"""
Created on Tue Oct 25 17:32:27 2016

@author: Martin
"""
import os; clear = lambda: os.system('cls'); clear() 
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as so 
import csv
import time
import msalg as ma

t0 = time.time()


name = 'SNLWT100A2_NWP113_2.mat'
name = 'SNLWT100A2_NWP113_damp005.mat'

var = ['U', 'Ui',  'Vel', 'Acc', 'T', 'CTable', 'X', 'N', 'OME']
MAT = so.loadmat(name, variable_names=var)
n0 = 45 - 1 # Base Node
udof = MAT['X']['CTable'][0][0][n0,0:3]
rdof = MAT['X']['CTable'][0][0][n0,3:6]
tt = MAT['X']['T'][0][0][1]
stps = len(tt) 

# Initialize output data
X0 = np.zeros((stps, 3))   # Base Point position
V0 = np.zeros((stps, 3))   # Base Point velocity (not used?)
A0 = np.zeros((stps, 3))   # Base Point acceration 
VO = np.zeros((stps, 3))   # Base Point angular
WO = np.zeros((stps, 3))   # Base Point angular
O0 = np.zeros((stps, 3))   # Base Point incremental rotation vector


#IC = np.zeros((stps, 13))
#IC[:,0] = tt
udof = udof - 1
rdof = rdof - 1
for idx, t in enumerate(tt):
    X0[idx] = MAT['N'][n0,1:4] + MAT['U'][udof, idx]
#    V0[idx] = MAT['Vel'][udof, idx]
    A0[idx] = MAT['Acc'][udof, idx]
    WO[idx] = MAT['OME'][rdof, idx]
    VO[idx] = MAT['Vel'][rdof, idx]
    O0[idx] = MAT['Ui'] [rdof, idx]
    

#IC[:,1:4] = A0
#IC[:,4:7] = O0
#IC[:,7:10] = VO
#IC[:,10:13] = X0
#
#
#plt.plot(tt.T, A0[:,2])
#plt.plot(tt.T, O0[:,0])
#
#
#    
#f = open('SNL100Tip_3.dat', 'wb')   
#f.write('Finite Motion Formated Data \n') 
#f.write('Acceleration, Incremental Rotation Vector, Rotation Vector Velocity, Absolute Position  \n') 
#writer = csv.writer(f, delimiter = '\t',quotechar =',',quoting=csv.QUOTE_MINIMAL)
#for idx, t in enumerate(tt):
#    writer.writerow(IC[idx,:])
#    
#print(time.time()-t0)
from scipy.interpolate import RegularGridInterpolator as rgi
CONS = {}
CONS['T0'] = 0  
CONS['TF'] = 100        
CONS['DT'] = 0.001
CONS['TIME'] = np.arange(CONS['T0'], CONS['TF']+CONS['DT'], CONS['DT']) 
CONS['LP'] = np.arange(0, 1, 1.0/ len(CONS['TIME']))        
CONS['TSTP'] = len(CONS['TIME']) 


X0i = np.zeros((CONS['TSTP'], 3))   # Base Point position
V0i = np.zeros((CONS['TSTP'], 3))   # Base Point velocity (not used?)
A0i = np.zeros((CONS['TSTP'], 3))   # Base Point acceration 
VOi = np.zeros((CONS['TSTP'], 3))   # Base Point angular
WOi = np.zeros((CONS['TSTP'], 3))   # Base Point angular
O0i = np.zeros((CONS['TSTP'], 3))   # Base Point incremental rotation vector


for idx, t in enumerate(CONS['TIME']):
   X0i[idx,0] =  np.interp(t, tt, X0[:,0]) 
   X0i[idx,1] =  np.interp(t, tt, X0[:,1])
   X0i[idx,2] =  np.interp(t, tt, X0[:,2])
   VOi[idx,0] =  np.interp(t, tt, VO[:,0]) 
   VOi[idx,1] =  np.interp(t, tt, VO[:,1])
   VOi[idx,2] =  np.interp(t, tt, VO[:,2])
   A0i[idx,0] =  np.interp(t, tt, A0[:,0]) 
   A0i[idx,1] =  np.interp(t, tt, A0[:,1])
   A0i[idx,2] =  np.interp(t, tt, A0[:,2])   


fac = 10 # Multiplicity between time steps
for idx in range(stps-1):
    do = O0[idx+1,:] / fac
    for jdx in np.arange(1,fac+1):
        i = (idx * fac) + jdx      
        O0acum = do * jdx
        O0i[i, :] = do # Rotacion incremental
        DT = ma.tangmap(O0acum)
        WOi[i, :] = np.dot(DT, VOi[i, :] )


   
plt.plot(tt, O0[:,0])
plt.plot(CONS['TIME'],O0i[:,1])

#plt.plot(tt, WO[:,0])
#plt.plot(CONS['TIME'],WOi[:,0])
#
#plt.plot(tt, X0[:,0])
#plt.plot(CONS['TIME'],X0i[:,0])
#
#plt.plot(tt, A0[:,0])
plt.plot(CONS['TIME'],A0i[:,:])
#
#plt.plot(tt, VO[:,0])
#plt.plot(CONS['TIME'],VOi[:,0])



IC = np.zeros((CONS['TSTP'], 13))
IC[:,0] = CONS['TIME']   
    
IC[:,1:4] = A0i
IC[:,4:7] = O0i
IC[:,7:10] = VOi
IC[:,10:13] = X0i




    
f = open('SNL100Tip_damp005_interp.dat', 'wb')   
f.write('Finite Motion Formated Data \n') 
f.write('Acceleration, Incremental Rotation Vector, Rotation Vector Velocity, Absolute Position  \n') 
writer = csv.writer(f, delimiter = '\t',quotechar =',',quoting=csv.QUOTE_MINIMAL)
for idx, t in enumerate(CONS['TIME']):
    writer.writerow(IC[idx,:])
#    
print(time.time()-t0)