# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 17:14:09 2016

@author: Martin
"""
import numpy as np
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import csv
dt = 1
r0 = 100.0
tf = 100.0
time = np.arange(0.0,tf,dt)

DO = np.array([0.0, 0.0, 0.2])
w = DO[2] / dt



X0 = np.zeros((len(time), 3))
V0 = np.zeros((len(time), 3))
A0 = np.zeros((len(time), 3))
VO = np.zeros((len(time), 3))


signal = np.zeros((len(time), 13))
signal[:,0] = time


fig = plt.figure()
ax = fig.gca(projection='3d')
for idx, t in enumerate(time):
    X0[idx, :] = np.array([r0 * np.cos(w*t), r0 * np.sin(w*t),  0.0])
    V0[idx, :] = np.array([-r0 * w * np.sin(w*t), r0 * w* np.cos(w*t),  0.0])
    A0[idx, :] = np.array([-r0 * w**2 * np.cos(w*t), -r0 * w**2 * np.sin(w*t),  0.0])
    VO[idx, :] = np.array([0.0, 0.0, w])
    signal[idx,1:4] = A0[idx, :]
    signal[idx,4:7] = DO
    signal[idx,7:10] = VO[idx, :]
    signal[idx,10:13] = X0[idx, :]
    ax.plot([X0[idx, 0], X0[idx, 0]], [X0[idx, 1], X0[idx, 1]] , [0,5.0])
    
    
## ----------  Export the file   ---------------- 
#f = open('SNL100TipSimple.dat', 'wb')   
#f.write('Finite Motion Formated Data \n') 
#f.write('Acceleration, Incremental Rotation Vector, Rotation Vector Velocity, Absolute Position Velocity \n') 
#writer = csv.writer(f, delimiter = '\t',quotechar =',',quoting=csv.QUOTE_MINIMAL)
#for idx, t in enumerate(time):
#    writer.writerow(signal[idx,:])
    
 
#np.savetxt('SNL100Tip.dat', signal, delimiter=',') 
#    
#plt.plot(time,X0)
    
#B  = np.loadtxt('circmotion.dat', dtype='float', skiprows=2)