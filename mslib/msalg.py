# -*- coding: utf-8 -*-
"""
Created on Tue Jul 08 20:43:22 2014

@author: gasmgpu1
"""

import numpy as np
# NORM OF A 3D VECTOR
def nrm(v):
    return (v[0]**2 + v[1]**2 + v[2]**2)**0.5    

# CROSS PRODUCT OF 3D VECTORS
def crs(v1, v2):
    x = ((v1[1] * v2[2]) - (v1[2] * v2[1]))
    y = ((v1[2] * v2[0]) - (v1[0] * v2[2]))
    z = ((v1[0] * v2[1]) - (v1[1] * v2[0]))
    v = np.array( [x,y,z] )
    return v

# DOT PRODUCT OF 3D VECTORS    
def dot(v1, v2):
    s = (v1[0] * v2[0]) + (v1[1] * v2[1])  + (v1[2] * v2[2]) 
    return s

# 3x3 IDENTITY MATRIX    
def eye():    
    eye = np.array( [[1.0, 0.0, 0.0],[0.0, 1.0, 0.0],[0.0, 0.0, 1.0]] )
    return eye    
    
def dia(n):    
    dia = np.array( [[n, 0.0, 0.0],[0.0, n, 0.0],[0.0, 0.0, n]] )
    return dia 

def do3(m1, m2, m3):
    return np.dot(m1,np.dot(m2, m3))

def do4(m1, m2, m3, m4):
    return np.dot( m1, np.dot(m2, np.dot(m3, m4) ) )

def skw(v):
    skw = np.array( [ [ 0.0, -v[2], v[1] ],
                      [ v[2], 0.0, -v[0] ],
                      [-v[1], v[0], 0.0 ] ] ) 
    return skw
# ANTICLOCKWISE NORMAL OF SEGMENTS   
def normal(p1, p2):
    dp =  p2 - p1
    sv = dp / nrm( dp )     # versor tangente 
    nvt = np.array( [ 0.0, sv[2], -sv[1] ])
    return (nvt/nrm(nvt))
    
def tangen(p1, p2):
    tv = (p2 - p1)
    return (tv/nrm(tv))

def versor(v):    
    return  (v / ( (v[0]**2 + v[1]**2 + v[2]**2)**0.5 )) 

def unitx():
    return np.array( [1.0, 0.0, 0.0] )

def unity():
    return np.array( [0.0, 1.0, 0.0] )

def unitz():
    return np.array( [0.0, 0.0, 1.0] )
    
    
#=================================================================
#              ROTATION MATRICES AND SO3 ALGEBRA
#=================================================================
def dT(dphi, phi):
    """Time derivative of T according to Makinen (50) - 
                  NO ANDA, O BIEN ANDA PARA EL OJETE """
    mphi = 1.0E-25 + nrm(phi)    
    C1 = ( mphi * np.cos(mphi) - np.sin(mphi) ) / (mphi**3)
    C2 = ( mphi * np.sin(mphi) + 2.0 * np.cos(mphi) - 2.0 ) / (mphi**4)
    C3 = ( 3.0 * np.sin(mphi) - 2.0 * mphi - mphi * np.cos(mphi) ) / (mphi**5)
    C4 = ( np.cos(mphi) - 1.0 ) / (mphi**2)
    C5 = ( mphi - np.sin(mphi) ) / (mphi**3)    

    dT = (    C1 * np.dot( phi, dphi) * eye() 
            - C2 * np.dot( phi, dphi) * skw(phi) 
            + C3 * np.dot( phi, dphi) * np.outer(phi, phi)
            + C4 * skw(dphi)  
            + C5 * ( np.outer(dphi, phi) + np.outer(phi, dphi) )  )
    return dT
    
def dTV(phi, V):
    """dT * V according to Cardona Thesis 6.15"""
    mphi = 1.0E-25 + nrm(phi)
    mV = 1.0E-25 + nrm(V)
    nphi = phi / mphi
    cphi = np.cos(phi)
    sphi = np.sin(phi)
    c1 = np.dot(nphi, V) / mphi
    c2 = sphi / mphi
    c3 = np.sin(0.5 * mphi) / (0.5 * mphi) 

    dTV = (   ( 1.0 + cphi - 2.0 * c2) * c1 * V 
              + ( (3.0 * c2 - cphi - 2.0) * c1**2 + (1.0 - c2) * (mV / mphi)**2  ) * phi
              + (c3**2 - c2 ) * c1 * np.cross(phi, V) )
    return dTV    
    
def xiT(V, phi):
    """Linearization of V in the direction of phi"""
    mphi = 1.0E-25 + nrm(phi)    
    C1 = ( mphi * np.cos(mphi) - np.sin(mphi) ) / (mphi**3)
    C2 = ( mphi * np.sin(mphi) + 2.0 * np.cos(mphi) - 2.0 ) / (mphi**4)
    C3 = ( 3.0 * np.sin(mphi) - 2.0 * mphi - mphi * np.cos(mphi) ) / (mphi**5)
    C4 = ( np.cos(mphi) - 1.0 ) / (mphi**2)
    C5 = ( mphi - np.sin(mphi) ) / (mphi**3)    
    Sphi = skw(phi)  
    SV   = skw(V)    
    xiT = (   C1 * np.outer( V, phi) 
            - C2 * np.outer( np.dot( Sphi, V), phi) 
            + C3 * np.outer( ( np.dot(phi, V) * phi), phi) 
            - C4 * SV 
            + C5 * ( np.dot(phi, V) * dia(1.0)  + np.outer(phi,V ) )  ) 
    return xiT

def tangmap(phi):

    phi_x = phi[0]       
    phi_y = phi[1]    
    phi_z = phi[2]
    penal = 1e-25         
    mphi = penal + nrm(phi);
    a1 = np.sin(mphi) / mphi;      
    a2 = (1 - np.cos(mphi)) / (mphi**2)      
    a3 = (mphi - np.sin(mphi)) / (mphi**3)    
    # Matriz T de Cardona (Ritto-Correa traspuesta)
    T = np.array( [ [          a1 + a3 * phi_x**2,           a3 * phi_x * phi_y - a2 * phi_z,          a2 * phi_y + a3 * phi_x * phi_z  ] ,  
                    [   a3 * phi_x * phi_y + a2 * phi_z,           a1 + a3 * phi_y**2       ,         -a2 * phi_x + a3 * phi_y * phi_z  ] ,
                    [  -a2 * phi_y + a3 * phi_x * phi_z,     a2 * phi_x + a3 * phi_y * phi_z,                  a1 + a3 * phi_z**2       ] ] )
    
    return T.T   

def expmap(phi):
    phi_x = phi[0]       
    phi_y = phi[1]    
    phi_z = phi[2]
    penal = 1e-25         
    mphi = penal + nrm(phi)
        
    R = np.array( [ [   1.0 + (2.0*(-(phi_y)**2.0 - (phi_z)**2.0)*np.sin(mphi/2.0)**2.0)/mphi**2.0,           (2.0*(phi_x)*(phi_y)*np.sin(mphi/2.0)**2.0)/mphi**2.0 - ((phi_z)*np.sin(mphi))/mphi   ,      (2.0*(phi_x)*(phi_z)*np.sin(mphi/2.0)**2.0)/mphi**2.0 + ((phi_y)*np.sin(mphi))/mphi ],
                    [ (2.0*(phi_x)*(phi_y)*np.sin(mphi/2.0)**2.0)/mphi**2.0 + ((phi_z)*np.sin(mphi))/mphi,     1.0 + (2.0*(-(phi_x)**2.0 - (phi_z)**2.0)*np.sin(mphi/2.0)**2.0)/mphi**2.0            ,       (2.0*(phi_y)*(phi_z)*np.sin(mphi/2.0)**2.0)/mphi**2.0 - ((phi_x)*np.sin(mphi))/mphi],
                    [ (2.0*(phi_x)*(phi_z)*np.sin(mphi/2.0)**2.0)/mphi**2.0 - ((phi_y)*np.sin(mphi))/mphi,    (2.0*(phi_y)*(phi_z)*np.sin(mphi/2.0)**2.0)/mphi**2.0 + ((phi_x)*np.sin(mphi))/mphi   ,       1.0 + (2.0*(-(phi_x)**2.0 - (phi_y)**2.0)*np.sin(mphi/2.0)**2.0)/mphi**2.0]]    ) 
    return R
    
    
    
#=================================================================
#                      INTERPOLATION ROUTINES
#=================================================================
def inter2(v, l):
    """ Linear Interpolation of Matrix or Vector v """
    vi = 0.5 * (v[0] + v[1])
    vd = (v[1] - v[0]) / l
    return vi, vd

#def interp(v, x):
#    """ Linear Interpolation of Matrix or Vector v """
#    
#    vi = 0.5 * (v[0] + v[1])
#    vd = (v[1] - v[0]) / l
#    return vi

def nint1d(pts, mm, ja):
    """ 1D Numerical Integration """
    # 1 Point Weight and Factor 
    if pts == 1:
        ips = 0.0
        wfs = 2.0        
    # 2 Point Weights and Factors        
    if pts == 2:
        ips = np.array( [ -(1.0/3.0)**0.5,  (1.0/3.0)**0.5 ] )
        wfs = np.array( [      1.0,             1.0        ] )
    # Matrix Integration    
    for i in np.arange(len(ips)):
        if mm.ndim == 1: # Vector Integration            
            x = mm.shape
            mi = np.zeros(mshape)
            sf = np.zeros(mshape)
            for j in np.arange( len(mm)/2 ):
                sf[j, j  ] = ( 1.0 - ips[i] ) / 2.0                     
                sf[j, j+6] = ( 1.0 + ips[i] ) / 2.0    
#            print mm
#            print sf.T
#            print np.dot(mm, sf)
#            print ja
#            print wfs[i]                  
            mi += wfs[i] * np.dot(sf.T, np.dot(sf, mm)) * ja                     
        if mm.ndim == 2: # Matrix Integration
            for j in np.arange( len(mm)/2 ):
                sf[j, j  ] = ( 1.0 - ips[i] ) / 2.0     
                sf[j, j+6] = ( 1.0 + ips[i] ) / 2.0                      
            mi += wfs[i] * np.dot(sf.T, np.dot(mm, sf)) * ja                      
    return mi # Integrated matrix