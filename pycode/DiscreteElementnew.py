# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 18:19:51 2015

@author: root
"""

from __future__ import print_function
import sys, os, numpy as np, scipy as sp, csv, time
from mslib import msalg as ma, msutil as mu
from code import SignalTools as st

class DiscreteElement(object):
    def __repr__(self): return 'Discrete_Element'
    def __init__(self, XM, ename, etype, edata):
        self.name = ename
        self.type = etype

class Set1(DiscreteElement):
    def __init__(self, XM, ename, etype, edata):
        DiscreteElement.__init__(self, XM, ename, etype, edata)
        self.LL = edata[0]
        self.TD = edata[1]
        self.TT = edata[2]
        self.CG = edata[3]
    
class SpringL1(DiscreteElement):
    def __init__(self, XM, ename, etype, edata):
        DiscreteElement.__init__(self, XM, ename, etype, edata)
        self.stif = edata[0]
        self.beta = 0.0
    def force(self, dx): 
        return self.stif * dx 
             
class SpringRM(DiscreteElement):
    def __init__(self, XM, ename, etype, edata):
        DiscreteElement.__init__(self, XM, ename, etype, edata)
        self.krea = edata[0]
        self.rmov = edata[1]
        self.coef = edata[2]
        self.stif = self.krea * self.coef * self.rmov**2.0
        self.beta = 0.0
    def force(self, dx): 
        return self.stif * dx              
             
class SpringNL1(DiscreteElement):
    def __init__(self, XM, ename, etype, edata):
        DiscreteElement.__init__(self, XM, ename, etype, edata)
        self.stif = edata[0]
        self.beta = edata[1]
    def force(self, dx): 
        return self.stif * ( dx + self.beta * dx**3.0 ) 
    def kk(self, du):
        return self.stif * (1.0 + self.beta * du**2.0 )
         
class SpringNL2(DiscreteElement):
    def __init__(self, XM, ename, etype, edata):
        DiscreteElement.__init__(self, XM, ename, etype, edata)
        self.stif = edata[0]
        self.beta = edata[1]
        self.gama = edata[2]
    def force(self, dx): 
        return self.stif * ( dx + self.beta * dx**3.0 + self.gama * dx**5.0 ) 
    def kk(self, du):
        return self.stif * (1.0 + self.beta * du**2.0 + self.gama * du**4.0)
    def roots(self, force):
        coefs = [self.gama, 0.0, self.beta, 0.0,  self.stif, -force]
        root = np.roots(coefs )
        return root
        
class DamperL1(DiscreteElement):
    def __init__(self, XM, ename, etype, edata):
        DiscreteElement.__init__(self, XM, ename, etype, edata)
        self.damp = edata[0]
        self.fric = edata[1]
    def force(self, dv):
        return self.damp * dv

class DamperF1(DiscreteElement):
    def __init__(self, XM, ename, etype, edata):
        DiscreteElement.__init__(self, XM, ename, etype, edata)
        self.damp = edata[0]     
        self.fric = edata[1]  

class DamperNL1(DiscreteElement):
    def __init__(self, XM, ename, etype, edata):
        DiscreteElement.__init__(self, XM, ename, etype, edata)
        self.damp = edata[0]

class Mass(DiscreteElement):
    def __init__(self, XM, ename, etype, edata):
        DiscreteElement.__init__(self, XM, ename, etype, edata)
        self.mass = edata[0]

class Inertia(DiscreteElement):
    def __init__(self, XM, ename, etype, edata):
        DiscreteElement.__init__(self, XM, ename, etype, edata)
        self.mass = edata[0]

        
class Coil(DiscreteElement):  #Continuos Coil
    def __init__(self, XM, ename, etype, edata, ewind, egage):
        DiscreteElement.__init__(self, XM, ename, etype, edata)
        Wire = AWGWire(egage)
        self.wire = Wire
        self.node = edata[0] # Location of coil
        self.ntrn = edata[1] # Number of coil turns
        self.diam = edata[2] # Coil mean diameter
        self.radi = self.diam / 2.0
        self.area = np.pi * self.radi**2.0
        self.ires = np.pi * self.diam * self.wire.dres
        self.resi = self.ntrn * self.ires # Calculated total resistance
        self.resi = edata[3] # Measured total resistance
        self.leng = edata[4]        
        self.indu = 4.0 * np.pi *1.0E-7 * (self.area * self.ntrn**2.0) / self.leng
        self.elas = 0.0

        self.x0 = XM.NODE[self.node][0] # Note that only the x coordinate is read, this means the coil formulation is relative
        nseg = 20
        self.nseg = nseg
        dtrn = self.ntrn / nseg #Sgment Turns
        dres = self.resi / nseg #Segment Resistance
        xbot = self.x0 - self.leng / 2.0
        xtop = self.x0 + self.leng / 2.0
        dx = self.leng / nseg
        self.table = np.zeros( (nseg, 7) )
        if ewind == 'POSITIVE':  # WINDING DIRECTION OF THE COIL
            self.wind = 1.0
        elif ewind == 'NEGATIVE':
            self.wind = -1.0
        for i in range(nseg):
            self.table[i,0] = xbot + i * dx  
            self.table[i,1] = xbot + (i + 1) * dx
            self.table[i,2] = xbot + i * dx + dx/2.0  # Mid coordinate of segment
            self.table[i,3] = dtrn 
            self.table[i,4] = self.area
            self.table[i,5] = dres
            
    def power(self, q, dq, ddq):
        power = ( self.elas * q + self.resi * dq + self.indu * ddq) * dq
        return power
    def dvolt(self, q, dq, ddq):
        dvolt =  self.elas * q + self.resi * dq + self.indu * ddq
        return dvolt
        
class Magnet(DiscreteElement):
    def __init__(self, XM, ename, etype, eflux, efile, edata):
        DiscreteElement.__init__(self, XM, ename, etype, edata)
        self.node = edata[0]
        self.x0 = XM.NODE[self.node][0] # Note that only the x coordinate is read, this means the Magnet formulation is relative
        if eflux == 'NUMERIC':
            filename = XM.inpdir + efile
            self.flaw = 'NUMERIC'
            self.data = np.genfromtxt (filename, delimiter=",", skip_header=3) # Avergarage Flux
            self.aflx = np.zeros( (len(self.data[:,0]), 2) )
            self.aflx[:,0] = self.data[:,0] /1000.0 
            self.aflx[:,1] = st.SavGol(self.data[:,1], 101, 2, deriv=0)    
            self.dflx = np.zeros((len(self.aflx[:,0]),2)) # Derivative of average flux
            self.dflx[:,0] = self.aflx[:,0] #LENGTH COORDINATE
            temp = st.SavGol(self.data[:,1], 101, 2, deriv=1)
            dx = self.aflx[1,0] - self.aflx[0,0] # Necesito el dx para dividir la derivada porqeu sale calculada tomando dx=1 por defecto
            self.dflx[:,1] = -st.SavGol(temp, 101, 2, deriv=0) / dx # Ojo que tira la derivada cambiada de signo
        else:
            self.flaw = 'ANALYTIC'
            self.perm = 4.0E-7 * 3.14159 
            self.mmom = edata[1]
            self.diam = edata[2]
            self.leng = edata[3]


class Resistance(DiscreteElement):
    def __init__(self, XM, ename, etype, edata):
        DiscreteElement.__init__(self, XM, ename, etype, edata)
        self.node = edata[0]   # Unused data
        self.resi = edata[1]

        
class ELoad(DiscreteElement):  # Electrical Load
    def __init__(self, XM, ename, etype, edata):
        DiscreteElement.__init__(self, XM, ename, etype, edata)
        self.node = int(edata[0]) # Location of coil
        self.resi = edata[1] # Resistance
        self.indu = edata[2] # Inductance
        self.capa = edata[3] # Capacitance
        if self.capa != 0.0: # Elastance
            self.elas = 1.0 / self.capa
        else:
            self.elas = 0.0
        
    def power(self, q, dq, ddq):
        power = ( self.elas * q + self.resi * dq + self.indu * ddq) * dq
        return power
    def dvolt(self, q, dq, ddq):
        dvolt =  self.elas * q + self.resi * dq + self.indu * ddq
        return dvolt
        
        
class AWGWire(object):
    def __init__(self, gauge):
#                    diam(m),   area(m2),   weight(N/m), break(N), resistance(Ohm/m)
        AWG = {'42':[0.0635E-3, 0.00321E-6, 0.1764E-3,   0.784,    5.00]}  # Resistance is measured
        self.gage = int(gauge)
        self.diam = AWG[gauge][0]
        self.area = AWG[gauge][1]
        self.wght = AWG[gauge][2]
        self.brkm = AWG[gauge][3]       
        self.dres = AWG[gauge][4]         
