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

zmaxtop = 1500
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

t = 2.28687167063
t = 6.
t = 4.5925892334506768

zmax = 1000
smax = 3950
h    = 1

Xe = [0,0,0]
Xr = [2500,2500,zmax]

Papri       = acls.canonical_shooting_angle(Xe,Xr)
Papri_compl = acls.canonical_shooting_angle(Xe,Xr,True)

sol_direct_sneldesc = rt.raytrace_SD1_frontend(Zextract,Cextract,
                                        Papri_compl[0],
                                        zmax,cumsum=0,out_path_length=1)

sol_invers_sneldesc = rt.raytrace_seek(np.array(Xe),np.array(Xr),
                                Zextract,Cextract,
                                5,85,verbose=0,fulloutput=0)                     
h        =  100
resotype = 'rkdp'

acls.raytrace_ODE_2or3d(ssf_opera,Xe,
                       Papri[0:2],
                       sol_direct_sneldesc[-1],
                       h=h,resotype=resotype,
                       return_mode='short')

#%%

reload(acls)
### SINGLE MODE ###
if 1:  
    h        =  100
    resotype = 'rkdp'

    sol , sol_direct = acls.raytrace_ODE_seek(Xe,Xr,ssf_opera,Papri,
                                 h,resotype,
                                 method='hybr')
                                 
    equiv_param = acls.equiv_get_param(Zextract,Cextract,Xr[-1])

    D,Ang, C, T = acls.equiv_inverse(equiv_param,Xe,Xr)    

    sol_opera = sol
    sol_invers_simil = acls.raytrace_ODE_2or3d(ssf_opera,Xe,
                                               sol_opera[0:2],
                                               sol_opera[2],
                                               h=h,resotype=resotype,
                                               return_mode='short')
                                               
    sol_invers_rkdp = acls.raytrace_ODE_2or3d(ssf_opera,Xe,
                                               sol_opera[0:2],
                                               sol_opera[2],
                                               h=h,resotype='rkdp',
                                               return_mode='short')                                    

    sol_invers_rkdp_friend = acls.raytrace_ODE_2or3d(ssf_opera,Xe,
                                               sol_opera[0:2],
                                               sol_opera[2],
                                               h=h,resotype='rkdp',
                                               return_mode='friendly') 
                                               
    print("coords du rec (simil mode)       " ,  sol_invers_simil[0])
    print("diff en tps du ping (simil mode) " ,  sol_invers_simil[-1] - sol_invers_sneldesc[-1])
    print("coords du rec (rkdp mode)        " ,  sol_invers_rkdp[0])
    print("diff en tps du ping (rkdp mode)  " ,  sol_invers_rkdp[-1]  - sol_invers_sneldesc[-1])

    ssf_opera.plot_a_ray(sol_invers_rkdp_friend[0],Xe,Xr)
    
#%%
### BATCH MODE ###
if 0:
    hlis     = [1 , 10 , 100 , 500 , 1000 , 2000 , 3000  ]
    restylis = ["rk2",'rk4','rkf45','rkck','rkdp','euler']

    argzlis = []
    argzlis_simple=[]
    for itertup in itertools.product(hlis,restylis):
        h , resotype  = itertup  
        h , resotype  = 1000 , 'rkf45'
        argzlis.append((Xe,Xr,ssf_opera,Papri,h,resotype))
        argzlis_simple.append((Xe,Xr,Papri,h,resotype))
    
    results_raw = []
    for arg in argzlis:
        results_raw.append(acls.raytrace_ODE_seek(*arg))
    
    res_invers_simil_lis = []
    res_invers_rkdp_lis  = []

    for res , arg in zip(results_raw,argzlis):
        h        = arg[4]
        resotype = arg[5]
        res_invers_simil = acls.raytrace_ODE_2or3d(ssf_opera,Xe,
                                               res[0:2],
                                               res[2]  ,
                                               h=h, resotype=resotype,
                                               return_mode='short')
                                               
        res_invers_rkdp  = acls.raytrace_ODE_2or3d(ssf_opera,Xe,
                                               res[0:2],
                                               res[2]  ,
                                               h=h, resotype='rkdp',
                                               return_mode='short')
                                               
                                          
        res_invers_simil_lis.append(res_invers_simil)
        res_invers_rkdp_lis.append(res_invers_rkdp)
        
        print("==============================================================")
        print("step :" , h , "resol. type : "     ,  resotype) 
        print("coords du rec (simil mode)       " , res_invers_simil[0])
        print("diff en tps du ping (simil mode) " , res_invers_simil[-1] - sol_invers_sneldesc[-1])
        print("coords du rec (rkdp mode)        " , res_invers_rkdp[0])
        print("diff en tps du ping (rkdp mode) "  , res_invers_rkdp[-1]  - sol_invers_sneldesc[-1])
        print("==============================================================")

if 0:
    outdir = "/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/working/COMPAR_RT_EIKO"
    outdir_work = gf.create_dir(os.path.join(gf.create_dir(outdir),
                                             gf.get_timestamp()))

    gf.pickle_saver(argzlis_simple,outdir_work,"argzlis_simple")
    gf.pickle_saver(results_raw,outdir_work,"results_raw")
    gf.pickle_saver(res_invers_rkdp_lis,outdir_work,"res_invers_rkdp_lis")    
    gf.pickle_saver(res_invers_rkdp_lis,outdir_work,"res_invers_simil_lis")
    gf.pickle_saver(res_invers_rkdp_lis,outdir_work,"res_invers_rkdp_lis")

#timer = timeit.timeit("acls.equiv_wrapper(Zextract,Cextract,Xe,Xr)","from __main__ import Zextract,Cextract,Xe,Xr,acls",number=1000)
