# -*- coding: utf-8 -*-
"""
Created on Mon Jun 08 10:47:45 2015
@author: root
"""

# FULL VEHICLE SEVEN DOF MODEL

from __future__ import print_function
import sys, os, numpy as np, csv, time
from mslib import msalg as ma, msutil as mu
import scipy as sp
import scipy.integrate as spint

class DiscreteSystem(object):
    def __repr__(self): return 'Discrete System Object'
    def __init__(self, XM, syskeys, snode, stif, damp, mass, solu, force):
        # Pass Element Objects to System Object
        self.CONS = {}
#        self.name = syskeys[0]
        self.name = XM.jname
        self.nodes = snode 
        self.nnode = len(self.nodes)
#        self.CONS['NDIM'] = self.nnode - 1
        self.CONS['SOLVER'] = solu[2]
        if syskeys[2] is None:
            self.CONS['FORMULATION'] = 'ABSOLUTE'
        else:
            self.CONS['FORMULATION'] = syskeys[2]
        self.CONS['MOTION'] = syskeys[3] # Motion flag ('finite' or 'local')
        self.Stif = []
        self.Damp = []
        self.Fric = []
        self.Mass = []
        self.Force = force
        self.Spring = []
        self.Damper = []
        
        for i in stif:
            self.Spring.append(XM.ELEM[i])
            self.Stif.append( XM.ELEM[i].stif )
        for i in damp:
            self.Damper.append(XM.ELEM[i])
            self.Damp.append( XM.ELEM[i].damp )
            self.Fric.append( XM.ELEM[i].fric )
        for i in mass:
            self.Mass.append( XM.ELEM[i].mass )
        
        # Form Time Table
        self.CONS['T0'] = solu[3][0]  
        self.CONS['TF'] = solu[3][1]        
        self.CONS['DT'] = solu[3][2]
        self.CONS['TIME'] = np.arange(self.CONS['T0'], self.CONS['TF'] + self.CONS['DT'], self.CONS['DT']) 
        lpstep = 1.0/ len(self.CONS['TIME'])
        self.CONS['LP'] = np.arange(0.0, 1.0 + lpstep, lpstep)        
        self.CONS['TSTP'] = len(self.CONS['TIME'])  - 1
        self.CONS['TDIM'] = len(self.CONS['TIME'])
        
        # Vectors Initialization
        self.U0 = np.zeros( (self.CONS['TDIM'], self.ndim)  ) # Base Displacement
        self.V0 = np.zeros( (self.CONS['TDIM'], self.ndim)  ) # Base Velocity
        self.A0 = np.zeros( (self.CONS['TDIM'], self.ndim)  ) # Base Acceleration
        self.X0 = np.zeros( (self.CONS['TDIM'], self.ndim)  ) # Base Position ( only for animation in relative formulation)
        
        # Current Variables
        self.Coor = np.array( [XM.NODE[i] for i in snode] )
        self.x = np.zeros( (self.CONS['TDIM'], self.nnode, 3) )    # Current Position
        self.o = np.zeros( (self.CONS['TDIM'], self.nnode, 3) ) # Current Incremental Rotation Vector
        self.e = np.zeros( (self.CONS['TDIM'], self.nnode, 3, 3) )    # Current Triad
        self.U = np.zeros( (self.CONS['TDIM'], self.nmdof) )       # Displacement
        self.V = np.zeros( (self.CONS['TDIM'], self.nmdof) )       # Velocity
        self.A = np.zeros( (self.CONS['TDIM'], self.nmdof) )       # Acceleration
        self.F = np.zeros( (self.CONS['TDIM'], self.nmdof) )       # Force
        self.FA = np.zeros( (self.CONS['TDIM'], self.nmdof) )       # Force
        self.FO = np.zeros( (self.CONS['TDIM'], self.nmdof) )       # Force
        self.P = np.zeros( self.CONS['TDIM'] ) # ???
        
        # Matrix Initialization
        self.KK = np.zeros((self.npdof, self.npdof))
        self.CC = np.zeros((self.npdof, self.npdof))
        self.MM = np.zeros((self.npdof, self.npdof))
        
        # Runge-Kutta Initialization
        self.ssdof = 2 * self.npdof
        self.SV = np.zeros( (self.CONS['TDIM'] , self.ssdof) ) # Runge-Kutta State Vector
        self.FS = np.zeros( (self.CONS['TDIM'] , self.ssdof) )                # Runge-Kutta Force ???
        
        # Initialize the Force Vector 
        self.InitForceVector()



    def StateFunction(self, sv):
    
        ts = self.CONS['TS']
        # Initialize State Matrices
        self.KS = np.zeros( (2*self.npdof,2*self.npdof) )    

        # Build state Space Matrix
        du = sv[self.suloc]
        dv = sv[self.svloc]        

        try:
            self.MagneticInduction(sv)
        except:
            pass
        
        # Update the force Vector
        self.UpdateForceVector(ts, sv)

        # Build State Space Submatrices
        KK = self.StiffnessMatrix(du)
        CC = self.DampingMatrix(dv)   
        FF = self.F[ts+1]
        MMinv = np.linalg.inv(self.MM) 
        MMKK = np.dot(MMinv, KK)
        MMCC = np.dot(MMinv, CC)
        MMFF = np.dot(MMinv, FF)
        
        # Assemble State Space Matrix
        ite = range(self.npdof)
        for i in ite:
            self.KS[i, i+self.npdof] = 1.0
            for j in ite:
                self.KS[i+self.npdof, j] = -MMKK[i,j]
                self.KS[i+self.npdof, j+self.npdof] = -MMCC[i,j]
                # Update the RK Force Vector (FS) F is the physical force vector
        
        self.FS[ts+1,self.npdof: 2*self.npdof] = MMFF
        
  
        # Stability of the method
        eig = np.real(np.linalg.eigvals(self.KS))
        dt = self.CONS['DT']
        ra = np.zeros(4)
        for idx, i in enumerate(eig):
            mu = i  * dt
            ra[idx] = 1.0 + mu + 0.5 * mu**2 + (1/6.0) * mu**3 + (1/24.0) * mu**4
        
        return (np.dot(self.KS, sv) + self.FS[ts+1])          


    def InitForceVector(self):     
        for f in self.Force:
            if f[1] == 'BASEMOTION':
                if self.CONS['FORMULATION'] == 'RELATIVE':
                    print('----- *BASE MOTION is still not compatible with RELATIVE formulation -----')
                print('   WARNING: Magnetic damping is not considered in the force term of *BASE MOTION cases -----')                
                if f[2] == 'SWEEP':
                    self.wf = np.zeros(self.CONS['TSTP'])
                    for ts in range(self.CONS['TSTP']-1):
                        for para in f[3]:
                            floc = para[0] - 1
                            A0 = para[1]                     
                            ti = self.CONS['TIME'][ts+1]
                            w0 = para[2] * 6.283185
                            wf = para[3] * 6.283185    
                            tf = self.CONS['TF']
                            t0 = self.CONS['T0']
                            aa =  ( wf - w0 ) / (tf - t0)
                            bb =  aa * t0 + w0
                            sin_a0_ti2 = np.sin(aa * ti**2 + bb * ti)
                            cos_a0_ti2 = np.cos(aa * ti**2 + bb * ti)
                            self.U0[ts+1,floc] += A0 * sin_a0_ti2
                            self.V0[ts+1,floc] += A0 * ( 2 * aa * ti + bb) * cos_a0_ti2 # ojo que wi = a * t, por eso aparece el 2
                            self.A0[ts+1,floc] += A0 * ( 2 * aa * cos_a0_ti2 - (2 * aa * ti + bb)**2 * sin_a0_ti2 )
                        elasticforce = self.Spring[0].force(self.U0[ts+1,floc])
                        dampingforce = (self.Damp[0] ) * self.V0[ts+1,floc] 
                        self.F [ts+1,floc] += elasticforce + dampingforce
            
            if f[1] == 'BASEACCELERATION':
            
                if self.CONS['FORMULATION'] == 'ABSOLUTE':
                    print('----- *BASE ACCELERATION is still not compatible with ABSOLUTE formulation -----')
                    sys.exit('*BASE ACCELERATION is still not compatible with ABSOLUTE formulation')
                if f[2] == 'SWEEP':
                    self.wf = np.zeros(self.CONS['TSTP'])
                    for ts in range(self.CONS['TSTP']-1):
                        for para in f[3]:
                            floc = (para[0] * 3 - 4) + para[1] 
                            A0 = para[2]                     
                            ti = self.CONS['TIME'][ts+1]
                            w0 = para[3] * 6.283185
                            wf = para[4] * 6.283185    
                            tf = self.CONS['TF']
                            t0 = self.CONS['T0']
                            aa =  ( wf - w0 ) / (tf - t0)
                            bb =  aa * t0 + w0         
                            self.A0[ts+1,floc] += A0 * np.sin(0.5 * aa * ti**2 + bb * ti)
                            self.wf[ts+1] = np.abs( ( aa * ti + bb) / 6.283185 ) 
                
                if f[2] == 'MEASUREMENT':
                    # Update the Time Table
                    from code import SignalTools as st
                    signal = st.SignalTools()
                    table = signal.ReadNational(f[4])
                    achannel = f[3][0][2]
                    floc = f[3][0][0] - 1
                    asignal = table[:, achannel]
                    tsignal = table[:,0]
                    factor = f[3][0][1]
                    for ts in range(self.CONS['TSTP']-1):
                        ti = self.CONS['TIME'][ts+1]
                        self.A0[ts+1, floc ] += np.interp(ti, tsignal, asignal) * factor

            if f[1] == 'BASEFINITEMOTION':
                if f[2] == 'MEASUREMENT':
                    # Finite kinematic excitation for arbitraty motion REC
                    from code import SignalTools as st
                    signal = st.SignalTools()       
                    table = signal.ReadNational(f[4])
                    floc = f[3][0][0] - 1
                    tt = table[:, 0]
                    # Check if time vector coincide (it is neccesary because rotation cannot be interpolated)
                    if (tt[1]-tt[0]) != self.CONS['DT']:
                        errmsg = '\n\n----- ERROR: Finite Motion Excitation time vector must be IDENTICAL TO the Discrete System solution vector -----'
                        print(errmsg)
                        raise Exception(errmsg)
                    self.A0 = table[:, 1:4] # Ojo, se sobreescribe el vector A0
                    self.O0 = table[:, 4:7]
                    self.VO = table[:, 7:10]
                    self.X0 = table[:, 10:13]
                                      
                    
            if f[1] == 'ACCELERATION':
                if f[2] == 'GRAVITY':
                    for acce in f[3]:
                        for iacc in self.A0:
                            iacc -=  acce[1:4]   # Le RESTO (ojo) la gravedad a la aceleracion de la base
                else:
                    for acce in f[3]:
                        for iacc in self.A0:
                            iacc +=  acce[1:4]   # Le suma la aceleracion impuesta a la aceleracion de la base

#         self.V0[:,floc] = sp.integrate.simps( self.CONS['TIME'], self.A0[ts+1,floc]) #NUMERICAL INTEGRATION OF ACCELERATION
#         self.U0[:,floc] = sp.integrate.simps( self.CONS['TIME'], self.V0[ts+1,floc]) 
            for idx, dof in enumerate(self.V0[0,:]):
                self.V0[:,idx] = spint.cumtrapz( self.A0[:,idx], self.CONS['TIME'], initial=0.0) #NUMERICAL INTEGRATION OF ACCELERATION
                self.U0[:,idx] = spint.cumtrapz( self.V0[:,idx], self.CONS['TIME'], initial=0.0)  


    def ConvergeStep(self,ts):
        pass




class REC3(DiscreteSystem):  # Recolector
    def __repr__(self): 
        return (self.name + '_REC3')

    def __init__(self, XM, para, snode, stif, damp, mass, coil, magn,  load, solu, force):
        self.nmdof =  1   # Mechanical DOFS
        self.npdof =  2   # MultiPhysical DOFS
        self.suloc = [0] # State space displacement location
        self.svloc = [1] # State space velocity location
        self.sqloc = [2] # State space charge location
        self.siloc = [3] # State space current location
        self.ndim = 3        
        
        DiscreteSystem.__init__(self, XM, para, snode, stif, damp, mass, solu, force)
        
        self.ELEM = {}
        elem = stif + damp + mass+ coil + magn + load
        for key in elem: 
            self.ELEM[key] = XM.ELEM[key]
        
        # Initialize Electrical Vectors 
        self.W = {}  # Power Element dictionary
        self.E = {}  # Voltage Element dictionary
        self.H = {}  #  Coil Flux density Element dictionary (flux derivative on every coil segment)
        for el in self.ELEM.keys():
            self.W[el] = np.zeros( self.CONS['TDIM'] ) # Electrical Power
            self.E[el] = np.zeros( self.CONS['TDIM'] ) # Voltage
                       
        self.I = np.zeros( self.CONS['TDIM'] ) # Current
        self.E['dphi'] = np.zeros( self.CONS['TDIM'] ) # Total induced voltage calculated from dphi
        self.dphi = np.zeros( self.CONS['TDIM'] )
        
#       Initialization of position and triads (Note that only the FIRST VECTOR of the triad is needed)
        for idx, node in enumerate(self.nodes): # Set the initial position
            self.x[0, idx, :] = XM.NODE[node] 
        V1 = self.x[0, 1, :] - self.x[0, 0, :] 
        self.le = ma.nrm(V1)
        self.E1 = V1 / self.le
        self.e[0] [0] [0, :] = self.E1 # Triad of the Base
        self.e[0] [1] [0, :] = self.E1 # Triad of the stacks is the same of the base       
        self.InertiaMatrix()
        # Initialize constant Matrices (in case of a new object with non constant matrices add return to the method definiton)
        self.Coil = []
        self.Magn = []
        self.Load = []
        for i in coil:
            self.Coil.append( XM.ELEM[i])
            self.H[i] = np.zeros( (XM.ELEM[i].nseg, self.CONS['TDIM'] ) )
        for i in magn:
            self.Magn.append( XM.ELEM[i] )
        for i in load:
            self.Load.append( XM.ELEM[i] )
        
            
        # System Electrical Properties    
        self.R = 0.0 # Total resistance in the system 
        self.L = 0.0 # Total inductance in the system 
        self.C = 0.0 # Total reactance  in in the system 
        self.S = 0.0 # Total elastance of the system
        for icoil in self.Coil:
            self.R += icoil.resi
            self.L += icoil.indu         
        for iload in self.Load:
            self.R += iload.resi
            self.L += iload.indu
            self.C += iload.capa
            self.S += iload.elas

        # System Properties       
        self.we = (self.S/self.L)**0.5 / (2.0 * np.pi)
        self.wn = (self.Stif[0]/self.Mass[0])**0.5 / (2.0 * np.pi)
        self.us = self.Spring[0].roots(self.Mass[0] * 9.8 )  # Static displacement in meters  
        print('   INFO: The resonant frequency of the electrical system is  ', self.we, ' Hz ')     
        print('   INFO: The resonant frequency of the mechanical system is  ', self.wn, ' Hz ') 
        # Stable time step of the electrical problem    
        sl = 2.78 / (self.R / self.L)
        
        if sl <= self.CONS['DT'] :
            print('   WARNING: The stable time step for RK4 is ', sl, ' but current step is ', self.CONS['DT'])
 
        
    def InertiaMatrix(self):
        skwE1 = ma.skw(self.E1)
        self.JJ = np.dot( skwE1, skwE1)
    
        
    def MagneticInduction(self, sv):
        Magn = self.Magn[0] 
        ts = self.CONS['TS']
        self.cm = 0.0         # EM Damping
        self.dphi[ts+1] = 0.0 # Integral of Magnetix flux derivative

        if Magn.flaw is 'NUMERIC':
            xm = sv[self.suloc] # magnet displacement (relative to coil since FORMULATION=RELATIVE)
            for Coil in self.Coil:
                for seg in Coil.table: 
                    xs = seg[2] # coil segment midpoint location
                    dflx = np.interp(xs, Magn.dflx[:,0] + xm, Magn.dflx[:,1]) # Derivative of flux at current segment location (flux function is shifted with magnet position)
                    sdphi = Coil.wind * (seg[3] * seg[4])  * dflx # Total segment integrated derivative
                    self.dphi[ts+1] += sdphi
                    Coil.flux()
            self.cm = sdphi**2 / self.R      # Magnetic Damping
             
            self.dphi_i = self.dphi[ts+1]
            self.cm_i = self.cm

    def StateFunction(self, sv):
        ts = self.CONS['TS']
        # Initialize State Matrices and Build state Space Matrix
        self.KS = np.zeros( (2*self.npdof,2*self.npdof) )    
        du = sv[self.suloc]
        dv = sv[self.svloc]        
        # Magnetics
        try:
            self.MagneticInduction(sv)
        except:
            pass
        # -- Damping parameters
        ff = 0.0
        ftol = 1.0E-5  # Friction Tolerance
        dvabs = np.absolute(dv)
        if ftol <= dvabs:
            ff = self.Fric[0] /  dvabs  # Friction Coefficient
        # -- Damping Matrix
        damp = self.Damp[0] + ff + self.cm_i
        self.KS[0,1] =   1.0
        self.KS[1,0] = - self.Spring[0].kk(du) / self.Mass[0]
        self.KS[1,1] = - damp / self.Mass[0]
        self.KS[2,3] =   1.0
        self.KS[3,1] = -  self.dphi_i / self.L
        self.KS[3,2] = - self.S / self.L
        self.KS[3,3] = - self.R / self.L
        # Update the force Vector
        self.UpdateForceVector(ts, sv)
        self.FS[ts+1,self.svloc] = -self.F[ts+1, 0] / self.Mass[0]
#       Stability of the method
        eig = np.real(np.linalg.eigvals(self.KS))
        dt = self.CONS['DT']
        ra = np.zeros(4)
        for idx, i in enumerate(eig):
            msu = i  * dt
            ra[idx] = 1.0 + msu + 0.5 * msu**2 + (1/6.0) * msu**3 + (1/24.0) * msu**4

        return (np.dot(self.KS, sv) + self.FS[ts+1]) 


    def UpdateForceVector(self, ts, sv):       
        FM = self.CONS['MOTION']  #Add support for harvester without finite motion
        if FM == 'FINITE':
            A0 = self.A0[ts+1, :]  # Base Acceleration
            O0 = self.O0[ts+1, :]  # Base Incremental Rotation Vector
            VO = self.VO[ts+1, :]  # Base Rotation Vector Velocity
            DT = ma.tangmap(O0)    # Incremental Tangential Map
            DR = ma.expmap (O0)    # Incremental Rotation tensor
            VW = np.dot( DT, VO )  # Angular Velocity
            um = sv[self.suloc]    # Magnet Displacement
            JJu = (self.le ) * self.JJ
            self.e[ts+1, 0, 0, :] = np.dot(DR, self.e[ts, 0, 0, :])
            self.e[ts+1, 1, 0, :] = np.dot(DR, self.e[ts, 0, 0, :])
            self.F[ts+1, 0]  = self.Mass[0] * np.dot(A0, self.e[ts+1, 0, 0, :]) + np.dot(VW, np.dot(JJu, VW))
            self.FA[ts+1, 0] = self.Mass[0] * np.dot(A0, self.e[ts+1, 0, 0, :]) 
            self.FO[ts+1, 0] = self.Mass[0] * np.dot(VW, np.dot(JJu, VW))
            
        elif FM == 'LOCAL':
            self.F[ts+1, 0] = self.Mass[0] * np.dot(self.A0[ts+1], self.E1) # Direction is not changed, so the initial director is used.
        


    def ConvergeStep(self, ts):   #Only for relative formulation     
        # Update the Kinematic Configuration
        sv = self.SV[ts+1]      # Converged state vector
        dsv = self.StateFunction(sv)  # Converged derivatives
        self.x[ts+1][0,:] = self.X0[ts+1][:]
        self.x[ts+1][1,:] = self.x[ts+1][0,:] + (self.le + self.U[ts+1, self.suloc] ) * self.e[ts+1][0][0,:]
        # Update Electrical variables     
        for key, el in self.ELEM.iteritems():
            q   = sv[self.sqloc]
            dq  = sv[self.siloc]
            ddq = dsv[self.siloc]
            try:
                self.W[key][ts] = el.power(q, dq, ddq)
                self.E[key][ts] = el.dvolt(q, dq, ddq)
            except:
                pass
            self.E['dphi'][ts]  = self.dphi[ts+1] * self.V[ts+1]
        
            



