# -*- coding: utf-8 -*-
"""
Created on Mon Jun 08 14:28:33 2015

@author: root
"""
import numpy as np

class ODEIntegrator(object):
    def __init__(self, DS):
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
        KE = DS.KK + a0 * DS.MM + a1 * DS.CC 
        DS.CONS['TS'] = 0
        for ts, ti in enumerate(DS.CONS['TIME'][0:-1]):
#            print('Running Time Step: ' + str(ts+1))
            DS.ForceVector()
            RE = DS.F[ts+1] + np.dot(DS.MM, a0 * DS.U[ts] + a2 * DS.V[ts] + a3 * DS.A[ts]) + np.dot(DS.CC, a1 * DS.U[ts] + a4 * DS.V[ts] + a5 * DS.A[ts])
            DS.U[ts+1] = np.linalg.solve(KE, RE)
            DS.A[ts+1] = a0 * (DS.U[ts+1] - DS.U[ts]) - a2 * DS.V[ts] - a3 * DS.A[ts]
            DS.V[ts+1] = DS.V[ts] + a6 * DS.A[ts] + a7 * DS.A[ts+1]
            DS.CONS['TS'] += 1


class EigenSolution(object):
    def __init__(self, DS):
        eigsys = np.dot( np.linalg.inv(DS.MM), DS.KK )
        sol = np.linalg.eig(eigsys)
        DS.FREQ = np.sqrt( sol[0]  )  / (2.0 * np.pi) 
        

        
#class RKF45(ODEIntegrator):
#    def __init__(self, DS):
#        ODEIntegrator.__init__(self)
#        dt = DS.CONS['TSTEP']
#        # Meirovitch Variant of RKF45
#        fxt = np.array([ U[ts], V[ts] ])
#        g1 = dt * fxt
#        g2 = dt * 
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        