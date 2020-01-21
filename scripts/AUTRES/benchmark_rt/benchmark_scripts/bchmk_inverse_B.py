# -*- coding: utf-8 -*-
"""
Created on Sat Feb 20 14:15:01 2016

@author: psakicki

TENTATIVE DE //ISATION
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
zm,cm   = rt.munk(zmaxtop,.1)
SSPmunk = acls.SSP(zm,cm,e=0)

#SSPbazik1.plot()
ssf_munk   = acls.make_SSF3D_from_SSP(SSPmunk,1)
ssf_bazik2 = acls.make_SSF3D_from_SSP(SSPbazik2)
ssf_bazik3 = acls.make_SSF3D_from_SSP(SSPbazik3)

ssf_opera = ssf_bazik2
ssf_opera = ssf_munk

Zextract , Cextract = ssf_opera.make_a_SSP(obj=False)

if 0:
    axiss0 =  0
    axiss1 =  1
    dmin  = -0.5
    dmax  =  0.5
    ssf_opera.add_a_gradient(axiss0,dmin,dmax)
    ssf_opera.add_a_gradient(axiss1,dmin,dmax)
# ====================== FIN DU SSP ======================

Xe = np.array([0,0,0])

hlis     = [1 , 10 , 100 , 500 , 1000 , 2000 , 3000  ]
hlis     = [1 , 5, 10, 20, 30, 40, 50, 60, 70, 80, 90]
hlis     = [1 , 2, 3, 5, 10, 20, 30, 40, 50,  75, 100]
hlis     = [1 , 2, 3, 5, 10, 20, 30 , 50,  75, 100]

restylis = ["rk2",'rk4','rkf45','rkck','rkdp','euler']
adaplis  = (True , )
xrlis    = [500 , 1000 , 2000 , 3000] 
yrlis    = [500 , 1000 , 2000 , 3000] 
yrlis    = [0] 
zrlis    = [500 , 1000 , 2000 , 3000] 

varargslis_proto = list(itertools.product(hlis,restylis,adaplis,
                                    xrlis,yrlis,zrlis))
                                    
varargslis = []
for e in varargslis_proto:
    x,y = e[-3],e[-2]
    if y > x:
        continue
    else:
        varargslis.append((e,))
        
print(len(varargslis) , 'experiences')

resultdic = dict()

N  = 5
Nf = float(N)


def wrapperfct(intup):
    tup = intup
    h , resotype , adap , xr , yr , zr = tup
    
    Xr = np.array([xr,yr,zr])
    
    Papri = acls.canonical_shooting_angle(Xe,Xr)
    #consititution du sous dico , aject ajonction préalable
    #d'infos supplémentaire
    tup = tup + (Papri,)
    
    ### SNELL DESCARTES
    print("SNELL DESCARTES")
    strt = time.time()
    for i in range(N):
        results_sd = rt.raytrace_seek(Xe,Xr,Zextract,Cextract,verbose=0,fulloutput=0)

    TIMER_sd  = time.time() - strt
    
    _,_,_,Ssd = rt.raytrace_SD1_frontend(Zextract,Cextract,results_sd[0],Xr[2],
                                         out_path_length=True)
                                         
    results_sd = results_sd + (Ssd , TIMER_sd/Nf)
    ### EIKONAL
    print("EIKONAL")
    try:
        strt         = time.time()
        for i in range(N):
            results_eiko = acls.raytrace_ODE_seek(Xe,Xr,ssf_opera,Papri,h,resotype,1,
                                              adaptative_step=adap,exectime_out=True)
        TIMER_eiko   = time.time() - strt
        results_eiko = results_eiko + (TIMER_eiko/Nf,)
    except:
        print(e)
        print("ERR : EIKO CRASH")
        results_eiko = np.nan

    ### EQUIV
    print("EQUIVALENT")
    paramtuple    = acls.equiv_get_param(Zextract,Cextract,Xr[-1])
    strt         = time.time()
    for i in range(N):
        results_equiv = acls.equiv_inverse(paramtuple,Xe,Xr)
    TIMER_equiv   = time.time() - strt
    results_equiv = results_equiv + (TIMER_equiv/Nf,)

    return tup , results_sd , results_eiko , results_equiv
    
pool   = mp.Pool(processes=3)

varargslis = list(reversed(varargslis))
#varargslis = varargslis[:2]

R  = [pool.apply_async(wrapperfct, args=ar) for ar in varargslis]
R2 = [e.get() for e in R]

outdir  = "/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/scripts/benchmark_rt/results/"
outpath = gf.pickle_saver(R2,outdir,gf.get_timestamp() + 'benchmark_rt_inverse')

print("results in" , outpath)