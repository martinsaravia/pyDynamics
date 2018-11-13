# -*- coding: utf-8 -*-
"""
Created on Sat Nov 05 18:47:06 2016

@author: root
"""


import os; clear = lambda: os.system('cls'); clear() 
import mslib.msutil as mu
import matplotlib.pyplot as plt, numpy as np
import matplotlib as mpl
import time
t0 = time.time()

import sys
sys.path.append(cpath[:-4] )
from code import DiscreteModel as DM
from code import SignalTools as st

fname = cpath[:-4] + 'output\\' + 'msharv023FMlag.pkl'


mu.openobj(fname)