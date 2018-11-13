
# ------------------------------------------------------------------------
#                           pyDynamics
#
#     GIMAP-CONICET 1D Discrete Dynamics Solver for Energy Harvesting
# ------------------------------------------------------------------------
 
#  License
#  This file is part of pyDynamics.
#
#  pyDynamics is free software: you can redistribute it and/or modify it
#  under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  pyDynamics is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#  for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with pyDynamics.  If not, see <http://www.gnu.org/licenses/>.
# ------------------------------------------------------------------------




import os; clear = lambda: os.system('cls'); clear() 
import  pickle 
from pycode.DiscreteModel import DiscreteModel 
from mslib import msutil as mu
import matplotlib as mpl
import matplotlib.pyplot as plt, numpy as np
#mpl.use('Qt4Agg') 
#font = {'family' : 'Times New Roman','weight' : 'normal','size': 11}
#mpl.rc('font', **font)
## -----------------------------------------------------------------------------


# ====================================================
name = 'yourFileName'
DM = DiscreteModel(name)
DM.Solve()
F_a = -DM.DSYS[name].cm * DM.DSYS[name].V[:,0]
uap = DM.DSYS[name].U[:,0]


name = 'yourFileName'
DM = DiscreteModel(name)
DM.Solve()
F_e = ( DM.DSYS[name].dphi  ) * DM.DSYS[name].SV[:,3]
uex = DM.DSYS[name].U[:,0]

t = DM.DSYS[name].CONS['TIME']

#plt.plot(t, 1000*DM.DSYS[name].U)   


plt.xlabel('Time (s)',   fontsize=11)
plt.ylabel('Force (N)',fontsize=11)
p1, = plt.plot(t, F_a, 'k--', label='Approximate' )
p2, = plt.plot(t, F_e, 'k-', label='Present' )
plt.legend(handles=[p2,p1]) 
plt.xlim([0.54,0.74])
plt.show()
mu.figexp('yourFileName', plt)    