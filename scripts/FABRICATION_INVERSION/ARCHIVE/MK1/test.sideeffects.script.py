# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 14:45:46 2015

@author: psakicki
"""

import acouclass as acls
import raytrace as rt
import numpy as np
import scipy.optimize as optimize

zm,cm = rt.munk(2000)
SSPmunk=acls.SSP(zm,cm,e=0)
ssf_munk = acls.make_SSF3D_from_SSP(SSPmunk)
Zextract , Cextract = ssf_munk.make_a_SSP()
ssf_munk.Interp_C([0,0,6500])

iposi = [0,0,300]
iangle = [80,90]
h = 500
smax = 5000
t = 5
resty = 'rkdp'

XYZ , _ , _ , _ , _  = acls.raytrace_ODE_2or3d(ssf_munk,iposi,iangle,smax,h,resotype=resty,return_mode='friendly')
BIGTUP   = rt.raytrace_SnellDesc2(Zextract,Cextract,zstart=iposi[-1],zmax=max(Zextract),tmax=t,theta=-iangle[0])

Z_clasik = np.cumsum(BIGTUP[0])
X_clasik = np.cumsum(BIGTUP[2])

plt.clf()
plt.plot(X_clasik,-Z_clasik,'rx-')
plt.plot(XYZ[:,1],-XYZ[:,2],'+')





