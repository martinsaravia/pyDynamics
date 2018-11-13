# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 00:04:44 2015

@author: root
"""
#from __future__ import print_function
import os, time
from multiprocessing import Pool
import pickle
from mslib import  msutil as mu
from pycode import DiscreteElement, DiscreteSystem, ODEIntegrator


class DiscreteModel(object):
    def __repr__(self): return 'PostRig_Model'
    def __init__(self, name):
        cpath = os.getcwd()
        self.t0 = time.time()
        self.jname = name
        self.inpdir = cpath + '\\' + 'input' + '\\'
        self.outdir = cpath + '\\' + 'output' + '\\'
        self.datadir = cpath + '\\' + 'data' + '\\'
        self.fname = self.inpdir + self.jname + '.pr7'
        self.ifile = open(self.inpdir + self.jname + '.dat' )
        self.ofile = open(self.outdir + self.jname + '.out', 'w')
        print ( ' ') 
        print ( ' ===================================================================== ')
        print ( ' ||               Discrete System Simulator v0.4.0                  || ')
        print ( ' ===================================================================== \n')        
        print ( ' ---> Initializing Discrete Model...')  
        self.CONS = {}
        self.CTRL = {}
        self.DSYS = {}
        self.NODE = {}
        self.ELEM = {}
        self.EXCI = {}
        self.ParseInputFile()
        
        
    def Solve(self, **kwargs):
        processing = 'Serial'
        if processing == 'Serial':
            for DS in self.DSYS.values():
                print (' ---> Solving system: ' + DS.name)  
                if DS.CONS['SOLVER'] == 'NEWMARK':
                    ODEIntegrator.LinearNewmark(DS)
                if DS.CONS['SOLVER'] == 'RK4':
                    ODEIntegrator.RungeKutta4(DS)
    #            DiscreteSolver.EigenSolution(DS)
    
                # Save object to dict
                if 'save' in kwargs:
                    if kwargs['save'] == 'YES':
                        filehandler = open(DS.name + '_1-40Hz_025g_aprox_damp' + '.obj', 'wb') 
                        pickle.dump(DS, filehandler) 
        else:
            with Pool(8) as p:              
                p.map(ODEIntegrator.RungeKutta4, self.DSYS.values())            
    
     
                
        
    def ParseInputFile(self):
        print ( ('-> Running DSS Parser for Input File:' + self.fname))   
        readflag = 'Succesfully'
        inpfile = self.ifile
        linea = mu.lread(inpfile)
        # Start processing
        while  linea[0] != '*END':
            # Process Nodes
            if linea[0] == '*NODE':
                linea = mu.lread(inpfile)
                while linea[0][0] != '*':    
                    import numpy as np
                    self.NODE[int(linea[0])] = np.array( [ float(i) for i in linea[1:] ] )                   
                    linea = mu.lread(inpfile)
            
            if linea[0] == '*ELEMENT':
                coman = linea
                ename = mu.kget(coman, 'NAME')
                etype = mu.kget(coman, 'TYPE')
                edata = []
                linea = mu.lread(inpfile)
                while linea[0][0] != '*':                   
                    edata.append( [ float(i) for i in linea ] )                   
                    linea = mu.lread(inpfile)
                if etype == 'SET1':
                    self.ELEM[ename] = DiscreteElement.Set1(self, ename, etype, edata[0])
                if etype == 'SPRINGL1':
                    self.ELEM[ename] = DiscreteElement.SpringL1(self, ename, etype, edata[0])      
                if etype == 'SPRINGPOLY':
                    self.ELEM[ename] = DiscreteElement.SpringPoly(self, ename, etype, edata[0])   
                if etype == 'SPRINGTABLE':
                    self.ELEM[ename] = DiscreteElement.SpringTable(self, ename, etype, edata)  
                if etype == 'DAMPERL1':
                    self.ELEM[ename] = DiscreteElement.DamperL1(self, ename, etype, edata[0])
                if etype == 'DAMPERF1':
                    self.ELEM[ename] = DiscreteElement.DamperF1(self, ename, etype, edata[0])
                if etype == 'DAMPERNL1':
                    self.ELEM[ename] = DiscreteElement.DamperNL1(self, ename, etype, edata[0])    
                if etype == 'MASS':
                    self.ELEM[ename] = DiscreteElement.Mass(self, ename, etype, edata[0])
                if etype == 'INERTIA':
                    self.ELEM[ename] = DiscreteElement.Inertia(self, ename, etype, edata[0])  
                if etype == 'RESISTANCE':
                    self.ELEM[ename] = DiscreteElement.Resistance(self, ename, etype, edata[0]) 
                if etype == 'ELOAD':
                    self.ELEM[ename] = DiscreteElement.ELoad(self, ename, etype, edata[0])  
                if etype == 'COIL':
                    ewind = mu.kget(coman, 'WINDING')
                    egage = mu.kget(coman, 'AWG')
                    self.ELEM[ename] = DiscreteElement.Coil(self, ename, etype, edata[0], ewind, egage)   
                if etype == 'MAGNET':
                    eflux = mu.kget(coman, 'FLUX')
                    efile = mu.kget(coman, 'FILE')
                    self.ELEM[ename] = DiscreteElement.Magnet(self, ename, etype, eflux, efile, edata[0])  


            # Process Discrete Systems
            elif linea[0] == '*DISCRETESYSTEM':
                sforce = []
#                sysname = mu.kget(linea, 'NAME')
                sysname = self.jname
                systype = mu.kget(linea, 'TYPE')
                sysform = mu.kget(linea, 'FORMULATION')  
                sysmotn = mu.kget(linea, 'MOTION')
                sysdamp = mu.kget(linea, 'DAMPING')
                syskeys = [sysname, systype, sysform, sysmotn, sysdamp]       
                ssets = [] # defined empty to avoid error
                if systype == 'CAR7' or systype == 'CAR5' or systype == 'DOF1' or systype == 'REC1' or systype == 'REC2'or systype == 'REC3':
                    snode = [int(i) for i in mu.lread(inpfile) ]
                    sstif = mu.lread(inpfile)
                    sdamp = mu.lread(inpfile)
                    smass = mu.lread(inpfile)
                    if systype == 'REC1' or systype == 'REC2'or systype == 'REC3':
                        scoil = mu.lread(inpfile)
                        smagn = mu.lread(inpfile)  
                        sload = mu.lread(inpfile)   # System Load
                    if systype == 'CAR7':
                        ssets = mu.lread(inpfile)
                    linea = mu.lread(inpfile)   
                while  linea[0] != '*ENDSYSTEM': 

                    if linea[0] == '*SOLUTION':
                        sname = mu.kget(linea, 'NAME')
                        stype = mu.kget(linea, 'TYPE')
                        salgo = mu.kget(linea, 'ALGORITHM')
                        spara = [ float(i) for i in mu.lread(inpfile) ]
                        ssolu = [sname, stype, salgo, spara]
                        linea = mu.lread(inpfile)
                        
                    elif linea[0] == '*BASEMOTION':
                        fname = mu.kget(linea, 'NAME')
                        fnatu = 'BASEMOTION'
                        ftype = mu.kget(linea, 'TYPE')
                        fpara = []
                        linea = mu.lread(inpfile)                        
                        while linea[0][0] != '*':
                            fpara.append( [int(linea[0]), float(linea[1]), float(linea[2]), float(linea[3]) ] )
                            linea = mu.lread(inpfile) 
                        # Update the force data
                        sforce.append( [fname, fnatu, ftype, fpara] )
                    
                    elif linea[0] == '*BASEFINITEMOTION':
                        fname = mu.kget(linea, 'NAME')
                        fnatu = 'BASEFINITEMOTION'
                        ftype = mu.kget(linea, 'TYPE')
                        ffile = mu.kget(linea, 'FILE')
                        fpara = []
                        linea = mu.lread(inpfile)                        
                        while linea[0][0] != '*':
                            fpara.append( [int(linea[0])] )
                            linea = mu.lread(inpfile) 
                        # Update the force data
                        sforce.append( [fname, fnatu, ftype, fpara, ffile] )
                    
                    elif linea[0] == '*BASEACCELERATION':
                        fname = mu.kget(linea, 'NAME')
                        fnatu = 'BASEACCELERATION'
                        ftype = mu.kget(linea, 'TYPE')
                        ffile = mu.kget(linea, 'FILE')
                        fpara = []
                        linea = mu.lread(inpfile)                        
                        while linea[0][0] != '*':
                            fpara.append( [int(linea[0]), int(linea[1]), float(linea[2]), float(linea[3]), float(linea[4]) ] )
                            linea = mu.lread(inpfile) 
                        # Update the force data
                        sforce.append( [fname, fnatu, ftype, fpara, ffile] )
                        
                    elif linea[0] == '*FORCE':
                        fname = mu.kget(linea, 'NAME')
                        fnatu = 'FORCE'
                        ftype = mu.kget(linea, 'TYPE')
                        fpara = []
                        linea = mu.lread(inpfile)                        
                        while linea[0][0] != '*':
                            fpara.append( [int(linea[0]), float(linea[1])] )
                            linea = mu.lread(inpfile) 
                        # Update the force data
                        sforce.append( [fname, fnatu, ftype, fpara] )
                        
                    elif linea[0] == '*GRAVITY':    
                        fname = mu.kget(linea, 'NAME')
                        fnatu = 'ACCELERATION'
                        ftype = 'GRAVITY'
                        fpara = []
                        linea = mu.lread(inpfile)     
                        while linea[0][0] != '*':
                            fpara.append( [int(linea[0]), float(linea[1]), float(linea[2]), float(linea[3])] )
                            linea = mu.lread(inpfile) 
                        # Update the force data
                        sforce.append( [fname, fnatu, ftype, fpara] )
                    elif linea[0][0:2] == '**':
                        linea = mu.lread(inpfile)

            elif linea[0] == '*ENDSYSTEM':
                self.CreateSystem(sysname, systype, syskeys, snode, sstif, sdamp, smass, scoil, smagn, sload, ssolu, sforce, ssets )
                linea = mu.lread(inpfile)
            
            elif linea[0] == '*PARAMETER':
                ptype = mu.kget(linea, 'TYPE')
                linea = mu.lread(inpfile)
                if ptype == 'SWEEP':
                    w0 = float(linea[0])
                    wf = float(linea[1])
                    dw = float(linea[2])
#                        np = int(linea[3])     # Number of sampling points
#                        ta = float(linea[4])   # tau                     
                    freqs = np.arange(w0, wf, dw)
                    
                    for fr in freqs:
                        sforce[0][3][0][3] = fr
                        sforce[0][3][0][4] = fr
                        sname = sysname + 'par=' + str(fr)
                        self.CreateSystem(sname, systype, syskeys, snode, sstif, sdamp, smass, scoil, smagn, sload, ssolu, sforce, ssets )
                linea = mu.lread(inpfile)

        
                    
            elif linea[0][0:2] == '**':
                linea = mu.lread(inpfile)
            # EOF  
            elif linea[0] == '*END':
                print('Parser has reached EOF ' + readflag)
                break  
            # Abort if unknown command or blank line
            else:                
                if not linea[0]:
                    print('    ERROR: Blank line found.' + '\n' + str(linea))
                    readflag = 'Unsuccesfully'
                    linea = mu.lread(inpfile)
                else:
                    print( '    ERROR: Unidentified Global Keyword: ' + linea[0])
                    readflag = 'Unsuccesfully'
                    linea = mu.lread(inpfile)
            
        print('-> Parser has reached EOF ' + readflag)
      
      
    def CreateSystem(self, sysname, systype, syskeys, snode, sstif, sdamp, smass, scoil, smagn, sload, ssolu, sforce, ssets ):   
          
        if systype == 'CAR7':  # PostRig
            self.DSYS[sysname] = DiscreteSystem.CAR7(self, syskeys, snode, sstif, sdamp, smass, ssets, ssolu, sforce)
        elif systype == 'CAR5': # Midget Start Model
            self.DSYS[sysname] = DiscreteSystem.CAR5(self, syskeys, snode, sstif, sdamp, smass, ssolu, sforce)
        elif systype == 'DOF1': # 1 Degree of Freedom Model (do not work  because 1x1 matrix size)
            self.DSYS[sysname] = DiscreteSystem.DOF1(self, syskeys, snode, sstif, sdamp, smass, ssolu, sforce)
        elif systype == 'REC3': #Recolector
            self.DSYS[sysname] = DiscreteSystem.REC3(self, syskeys, snode, sstif, sdamp, smass, scoil, smagn, sload, ssolu, sforce)
        self.DSYS[sysname].name = sysname