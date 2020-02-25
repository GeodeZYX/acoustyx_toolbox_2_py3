# -*- coding: utf-8 -*-
"""
Created on Sun Mar  1 18:28:47 2015

@author: psakicki
"""

import acouclass as acls
import raytrace as rt
import numpy as np
import scipy.optimize as optimize
#from megalib import *
import time
import copy
import itertools
import timeit
import pickle

import multiprocessing as mp

from megalib import *

reload(acls)

import warnings
warnings.simplefilter('always', UserWarning)


# ==================== zone de test ====================
Zbazik  = np.arange(0,5000,1)
Cbazik1 = 1500 + 20  * np.cos(Zbazik  * (0.001))
Cbazik2 = 1500 + 20  * np.sin(Zbazik  * (0.001))
Cbazik3 = 1500 + 500 * np.sin(Zbazik  * (0.01))

SSPbazik1 = acls.SSP(Zbazik,Cbazik1,e= 0)
SSPbazik2 = acls.SSP(Zbazik,Cbazik2,e=10)
SSPbazik3 = acls.SSP(Zbazik,Cbazik3,e=10)

zmaxtop = 4000
zm,cm   = rt.munk(zmaxtop)
SSPmunk = acls.SSP(zm,cm,e=0)

#SSPbazik1.plot()
ssf_munk   = acls.make_SSF3D_from_SSP(SSPmunk)
ssf_bazik2 = acls.make_SSF3D_from_SSP(SSPbazik2)
ssf_bazik3 = acls.make_SSF3D_from_SSP(SSPbazik3)

ssf_opera = ssf_bazik2
ssf_opera = ssf_munk

Zextract , Cextract = ssf_opera.make_a_SSP(obj=False)

if 0:
    axiss0 =  0
    axiss1 =  1
    dmin  = -0.6
    dmax  =  0.6
    ssf_opera.add_a_gradient(axiss0,dmin,dmax)
    ssf_opera.add_a_gradient(axiss1,dmin,dmax)
# ====================== FIN DU SSP ======================
    
z0=0
r0=0
theta0 = 80

theta_4_integ  =   theta0
theta_4_clasik = - theta0
smax = 3456
smax = 4456
h=100

h1000 = 1000
iposi = [0,0,0]
iangle = [-theta0,90]
restylis = ['rkck']
restylis = ["rk2",'rk4','rkf45','rkck','rkdp','euler']

resultlis = []
result1000lis = []

zmax = 3000
smax = 3950
h    = 1

Xe = np.array([0,0,0])
Xr = np.array([2500,2500,zmax])
Xr = np.array([1000,1000,zmax])


Papri = acls.canonical_shooting_angle(Xe,Xr)

sol_dir_SD = rt.raytrace_SD1_frontend(Zextract,Cextract,
                                      Papri[0],zmax,cumsum=0,
                                      out_path_length=True,
                                      shift_90deg = True)


                              
sol_dir_eiko = acls.raytrace_ODE_2or3d(ssf_opera,Xe,
                                  (Papri[0],0),
                                  sol_dir_SD[-1],
                                  h=h,resotype='rkdp',
                                  return_mode='short')  

sol_dir_eiko = acls.raytrace_ODE_2or3d(ssf_opera,Xe,
                                  (Papri[0],0),
                                  sol_dir_SD[-1],
                                  h=h,resotype='rkdp',
                                  return_mode='short')  
                                  
equiv_param   = acls.equiv_get_param(Zextract,Cextract,Xr[-1])
sol_dir_equiv = acls.equiv_direct(equiv_param,Xe,Papri[0],sol_dir_SD[-1],'s')
sol_inv_eiko  = acls.raytrace_ODE_seek(Xe,Xr,ssf_opera,Papri,h,'rkdp',True)

#%%
sol_inv_SD = rt.raytrace_seek(Xe,Xr,Zextract,Cextract,
                              5,85,verbose=0,fulloutput=0)
                               
equiv_param = acls.equiv_get_param(Zextract,Cextract,Xr[-1])

sol_inv_equiv = acls.equiv_inverse(equiv_param,Xe,Xr)

hlis     = [1 , 10 , 100 , 500 , 1000 , 2000 , 3000  ]
restylis = ["rk2",'rk4','rkf45','rkck','rkdp','euler']

varargslis = list(itertools.product(hlis,restylis))

N = 20 

resultslis = []

for ii , (h , resotype) in enumerate(varargslis):
    if ii != 10 and False:
        continue
    print('exp',ii,'h = ', h , 'resotype = ' , resotype)
    results = acls.raytrace_ODE_seek(*(Xe,Xr,ssf_opera,Papri,h,resotype,1))
    TIMER = timeit.timeit("acls.raytrace_ODE_seek(*(Xe,Xr,ssf_opera,Papri,h,resotype))",
    "from __main__ import Xe,Xr,ssf_opera,Papri,h,resotype,acls",number=N)
    resultslis.append(results + (TIMER,N,h,resotype,TIMER/float(N)))

outdir = "/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/scripts/benchmark_rt/results/"
gf.pickle_saver(resultslis,outdir,gf.get_timestamp() + 'benchmark_rt_direct')
