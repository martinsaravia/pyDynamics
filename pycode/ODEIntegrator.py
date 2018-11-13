# -*- coding: utf-8 -*-
"""
Created on Mon Jun 08 14:28:33 2015

@author: root
"""
import numpy as np

class ODEIntegrator(object):
    def __init__(self, DS):
        DS.CONS['TS'] = 0 
        pass
        
    
class LinearNewmark(ODEIntegrator):
    def __init__(self, DS):
        ODEIntegrator.__init__(self, DS)
        # Bathe Variant of Newmark
        dt = DS.CONS['DT']
        bet = 0.25
        gam = 0.50 
        a0 = 1.0 / (bet * dt**2)
        a1 = gam / (bet * dt)
        a2 = 1.0 / (bet * dt)
        a3 = 1.0 / (2.0 * bet) - 1.0
        a4 = gam / bet - 1.0
        a5 = 0.5 * dt * (gam / bet - 2.0)
        a6 = dt * (1.0 - gam)
        a7 = gam * dt  
        # Form Effective Stiffness
        DS.Tangent()
        KE = DS.KK + a0 * DS.MM + a1 * DS.CC 
        DS.ForceVector()
        for ts, ti in enumerate(DS.CONS['TIME'][0:-1]):
            RE = DS.F[ts+1] + np.dot(DS.MM, a0 * DS.U[ts] + a2 * DS.V[ts] + a3 * DS.A[ts]) + np.dot(DS.CC, a1 * DS.U[ts] + a4 * DS.V[ts] + a5 * DS.A[ts])
            if DS.ndof is 1:
                DS.U[ts+1] = RE / KE 
            else:
                DS.U[ts+1] = np.linalg.solve(KE, RE) 
            DS.A[ts+1] = a0 * (DS.U[ts+1] - DS.U[ts]) - a2 * DS.V[ts] - a3 * DS.A[ts]
            DS.V[ts+1] = DS.V[ts] + a6 * DS.A[ts] + a7 * DS.A[ts+1]
            DS.CONS['TS'] += 1

        
class RungeKutta4(ODEIntegrator):
    def __init__(self, DS):
        ODEIntegrator.__init__(self, DS)
        DS.SV = np.zeros( (DS.CONS['TDIM'] , DS.ssdof) ) # State Vector
        dt = DS.CONS['DT']
        for ts in range(DS.CONS['TSTP'] ):
            k1 = dt * DS.StateFunction( DS.SV[ts,:] )
            k2 = dt * DS.StateFunction( DS.SV[ts,:] + 0.5 * k1)
            k3 = dt * DS.StateFunction( DS.SV[ts,:] + 0.5 * k2)
            k4 = dt * DS.StateFunction( DS.SV[ts,:] + k3)
            DS.SV[ts+1,:] = DS.SV[ts,:] + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
            
            DS.U[ts+1,:] = DS.SV[ts+1, DS.suloc]
            DS.V[ts+1,:] = DS.SV[ts+1, DS.svloc]
            try:
                DS.I[ts+1,:] = DS.SV[ts+1, DS.siloc]    # Electrical DOF
            except:
                pass
            # Converence Activities
            DS.ConvergeStep(ts)
            
            DS.CONS['TS'] += 1
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        