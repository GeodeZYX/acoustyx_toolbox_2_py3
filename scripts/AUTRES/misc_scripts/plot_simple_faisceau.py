# -*- coding: utf-8 -*-
"""
Éditeur de Spyder

Ceci est un script temporaire.
"""

import raytrace as rt
import SSP
import scipy.interpolate

Z,C = SSP.munk(6000)
ZC = np.loadtxt('/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/exemple/SSP/SSP_NOAA_dep5781_20050827000000')
Z,C = ZC[:,0] , ZC[:,1]

Tetalon = np.arange(0,10,1)

Xt_stk = []
Zt_stk = []


for t in Tetalon:    
    Xt_stk.append(list())
    Zt_stk.append(list())


Tstk = []
Xstk = []
Zstk = []

for ang in np.arange(0,90,2):
    OUT = rt.raytrace_SnellDesc(Z,C,ang,5000)
    
    X = np.cumsum(OUT[2])
    T = np.cumsum(OUT[3])
    Z = OUT[0]
    Tstk.append(T)
    plt.plot(X,-Z)
    
    IX = scipy.interpolate.interp1d(T,X,bounds_error=0)
    IZ = scipy.interpolate.interp1d(T,Z,bounds_error=0)
    
    for it, t in enumerate(Tetalon):
        Xt_stk[it].append(IX(t))
        Zt_stk[it].append(-IZ(t))
        


for it , t  in enumerate(Tetalon):
    if 0:
        continue
    else:
        plt.plot(Xt_stk[it] , Zt_stk[it] ,'r.-' )
    
        
    
plt.axis('equal')