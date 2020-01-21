# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 16:18:08 2016

@author: psakicki
"""

import acouclass as acls
import raytrace as rt
import numpy as np
import scipy.optimize as optimize
#from megalib import *
import time
import copy

reload(acls)

# ==================== zone de test ====================

z0=0
r0=0
theta0 = 80

theta_4_integ =  theta0
theta_4_clasik = - theta0


Zbazik = np.arange(0,5000,1)
Cbazik1 = 1500 + 20 * np.cos(Zbazik  * (0.001))
Cbazik2 = 1500 + 20 * np.sin(Zbazik  * (0.001))

SSPbazik1 = acls.SSP(Zbazik,Cbazik1,e=0)
SSPbazik2 = acls.SSP(Zbazik,Cbazik2,e=10)
zm,cm = rt.munk(6000)

SSPmunk=acls.SSP(zm,cm,e=0)




#SSPbazik1.plot()
#SSPbazik2.plot()

#ssf1 = acls.make_SSF3D_from_SSP(SSPbazik1,xmin,xmax,ymin,ymax,xstep,ystep)
#ssf2 = acls.make_SSF3D_from_SSP(SSPbazik2,xmin,xmax,ymin,ymax,xstep,ystep)
ssf_munk = acls.make_SSF3D_from_SSP(SSPmunk)

OUT = acls.raytrace_ODE_stop(ssf_munk, np.array([0,0]), (theta_4_clasik,) , 
                       stop_max= 6 , path_step = 1 ,
                       stop_mode = "t" , resotype = 'rk4',  
                       theta_orient_hor = True ,
                       return_mode = 'short' , zmax_inp = 0 , 
                       s_init = 10000 , tolerance = 10**-8) 
                       
#                       
#BIGTUP_EIKO = acls.raytrace_ODE_2or3d(ssf_munk,np.array([0,0]),
#                              (theta_4_clasik,), 5500,
#                              h=1,resotype='rk4',
#                              return_mode='short',
#                              zmax_inp = 5000 )

Zextract , Cextract = ssf_munk.make_a_SSP(obj=False)
BIGTUP = rt.raytrace_SD3_frontend(Zextract,Cextract,0,np.max(Zextract),tmax=6,
                                  theta=-theta_4_clasik,cumsum = 1,plotflag=0)
X_clasik , Z_clasik , T_clasik = BIGTUP
X_clasik[-1] , Z_clasik[-1]