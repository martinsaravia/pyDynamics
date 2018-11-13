# -*- coding: utf-8 -*-
"""
Created on Thu Aug 28 18:15:42 2014

@author: gasmgpu1
"""
import re
import numpy as np
# Get keyword 
def kget(line, keyword):
    "Get keyword from line. Return word after *KEYWORD="
    # kstr = '(?<=' + keyword + '=)\w+'
    # Same as above but match the point in order to read numbers after strings
    kstr = '(?<=' + keyword + '=)[a-zA-Z0-9._]+'
    for i in line:
        match = re.search(kstr, i)
        if match is not None:
            break
    if match is None:
        print ('   WARNING: Unable to find ', keyword, ' keyword in line', line)
        return match
    else:
        return match.group()

#Get float from keyword 
def fget(line, keyword):
    "Get keyword from line. Return float after *KEYWORD="
    kstr = '(?<=' + keyword + '=)[a-zA-Z0-9._]+'
    for i in line:
        match = re.search(kstr, i)
        if match is not None:
            break
    if match is None:
        print ('   WARNING: Unable to find ', keyword, ' keyword in line', line)
        return match
    else:
        return float(match.group())

   
def lread(fobj): 
    "Line Read From File Object. Return line string without spaces"            
    linea = fobj.readline()
    linea = linea.replace(" ", "")
    linea = linea.replace("\n", "")
    linea = linea.replace("\t", "")
    linea = linea.split(',')
    return linea
    
def fread(fobj):
    "Float Read From File Object. Return list of floats"
    linea = lread(fobj)
    linealist = map(float, linea.split(','))
    return linealist

def iread(fobj):
    "Integer Read From File Object. Return list of Integers"
    linea = lread(fobj)
    linealist = map(int, linea.split(','))
    return linealist
    
def tab_read(fobj):
    import itertools
    "Read VABS tab delimited file"
    linea = fobj.readline()
#    linea = linea.replace(" ", "")
    linea = linea.replace("\n", "")
    linea = linea.split(' ')
    linea = [ i.split('\t') for i in linea ]
    linea = list(itertools.chain(*linea))
    linea = filter(None, linea)
    return linea     

def saveobj(obj, filename):
    import pickle
    with open(filename, 'wb') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)

def openobj(filename):
    import pickle
    with open(filename, 'rb') as input:
        return pickle.load(input)
        

def figexp(name, pltobj):
    pltobj.savefig(name+'.eps', format='eps', dpi=1200)
    pltobj.savefig(name+'.pdf', format='pdf', dpi=1200)
    pltobj.savefig(name+'.png', format='png', dpi=1200)