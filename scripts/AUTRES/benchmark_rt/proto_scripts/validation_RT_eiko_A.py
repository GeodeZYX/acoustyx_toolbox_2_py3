# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 14:57:34 2016

@author: psakicki
"""

import acouclass as acls
import raytrace as rt
import numpy as np
import scipy.optimize as optimize
#from megalib import *
import time
import copy



# ==================== zone de test ====================

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

ssf_munk.Z


#fig = mlab.figure()
#ssf1.plot()
#time.sleep(2)
##mlab.figure()
#mlab.clf()
#ssf2.plot()

#ssf4d1 = acls.SSF4D([ssf2])
#ssf4d1.add_ssf_in_ssf3d_list(ssf1)
#ssf4d1.ssf3d_list

#mlab.figure()
#ssf_out.plot()
#mlab.figure()
#ssf2.plot()

Zextract , Cextract = ssf_munk.make_a_SSP(obj=False)

z0=0
r0=0
theta0 = 80

theta_4_integ =  theta0
theta_4_clasik = - theta0
smax = 3456
smax = 4456
h=100

h1000 = 1000
iposi = [0,0,0]
iangle = [-theta0,90]
restylis = ['rkck']# , 'rk2'  ,  , 'rkck' ]
restylis = ["rk2",'rk4','rkf45','rkck','rkdp','euler']# , 'rk2'  ,  , 'rkck' ]

resultlis = []
result1000lis = []

#resty='euler'
#YYY,s,c,t,HH = acls.raytrace_ODE_2or3d(ssf_munk,iposi,iangle,smax,h,resotype=resty,return_mode='full')

for resty in restylis:
    print('=====' , resty , '=====')
    restuptmp = acls.raytrace_ODE_2or3d(ssf_munk,iposi,iangle,smax,h,resotype=resty,return_mode='full')
#    restuptmp1000 = acls.raytrace_ODE_2or3d(ssf_munk,iposi,iangle,smax,h1000,resotype=resty,return_mode='full')
    resultlis.append(restuptmp)
#    result1000lis.append(restuptmp1000)
    
resultlis[1][1]
       
t = 2.28687167063
t = 4.5925892334506768
t = 4.

zmax = 4999

#BIGTUP   = rt.raytrace_SnellDesc2(Zextract,Cextract,zstart=0,zmax=max(Zextract),tmax=t,theta=theta_4_clasik)
#X_clasik , Z_clasik , T_clasik  = rt.raytrace_SnellDesc_frontend(Zextract,Cextract,0,max(Zextract),tmax=t,theta=theta_4_clasik)

print('=====','classic','=====')

BIGTUP = rt.raytrace_SD3_frontend(Zextract,Cextract,0,zmax,tmax=t,
                                  theta=-theta_4_clasik,cumsum = 1,plotflag=0)
    
X_clasik , Z_clasik , T_clasik = BIGTUP

np.savetxt('/home/psakicki/Zmunk1',Zextract)
np.savetxt('/home/psakicki/Cmunk1',Cextract)

print("===== chargement du matlab =====")
Z_ml = np.loadtxt('/home/psakicki/Zmunk_matlab_1')
X_ml = np.loadtxt('/home/psakicki/Xmunk_matlab_1')

Z_py = BIGTUP[1]
X_py = BIGTUP[0]

print('=====','eiko direct','=====')
returnmode = 'full'
BIGTUP_EIKO = acls.raytrace_ODE_2or3d(ssf_munk,np.array([0,0]),
                              (theta_4_clasik,),  6.09480479e+03 ,
                              h=1,resotype='rk4',
                              return_mode=returnmode,
                              zmax_inp = zmax )

if returnmode == 'friendly':
    X_eiko = BIGTUP_EIKO[0][:,0]
    Z_eiko = BIGTUP_EIKO[0][:,1]
elif returnmode ==  'short':
    X_eiko = BIGTUP_EIKO[0][0]
    Z_eiko = BIGTUP_EIKO[0][1] 
elif returnmode ==  'full':
    X_eiko = BIGTUP_EIKO[0][:,0]
    Z_eiko = BIGTUP_EIKO[0][:,1]    
    
    
print('===== eiko indirect =====')

BIGTUP_EIKO_INDIRECT = acls.raytrace_ODE_stop(ssf_munk, np.array([0,0]), 
                                              (theta_4_clasik,) , 
                       stop_max= t , path_step = 1 ,
                       stop_mode = "t" , resotype = 'rk4',  
                       theta_orient_hor = True ,
                       return_mode = returnmode , zmax_inp = zmax , 
                       s_init = 10000 , tolerance = 10**-8) 


if returnmode == 'friendly':
    X_eiko_indi = BIGTUP_EIKO_INDIRECT[0][:,0]
    Z_eiko_indi = BIGTUP_EIKO_INDIRECT[0][:,1]
elif returnmode ==  'short':
    X_eiko_indi = BIGTUP_EIKO_INDIRECT[0][0]
    Z_eiko_indi = BIGTUP_EIKO_INDIRECT[0][1]  
elif returnmode ==  'full':
    X_eiko_indi = BIGTUP_EIKO_INDIRECT[0][:,0]
    Z_eiko_indi = BIGTUP_EIKO_INDIRECT[0][:,1] 
    S_eiko_indi = BIGTUP_EIKO_INDIRECT[1]


plt.clf()
plt.plot(X_ml,-Z_ml,'-g+')
plt.plot(X_py,-Z_py,'rx-')
plt.plot(X_eiko,-Z_eiko,'kv-')
plt.plot(X_eiko_indi,-Z_eiko_indi,'y^-')