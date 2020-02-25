# -*- coding: utf-8 -*-
"""
Created on Sat Feb 20 14:15:01 2016

@author: psakicki

DISCONTINUÉ CAR NON //ISÉ
ET C'EST LEEEENNNNNNT
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

Xe = np.array([0,0,0])

hlis     = [1 , 10 , 100 , 500 , 1000 , 2000 , 3000  ]
restylis = ["rk2",'rk4','rkf45','rkck','rkdp','euler']
adaplis  = (True , )
xrlis    = [500 , 1000 , 2000 , 3000] 
yrlis    = [500 , 1000 , 2000 , 3000] 
zrlis    = [500 , 1000 , 2000 , 3000] 

varargslis = list(itertools.product(hlis,restylis,adaplis,
                                    xrlis,yrlis,zrlis))

resultdic = dict()

N  = 1
Nf = float(N)


for ii , tup in enumerate(varargslis):
    
    h , resotype , adap , xr , yr , zr = tup
    print("exp",ii,tup)
    
    Xr = np.array([xr,yr,zr])

    if ii == 1001:
        break    
    if ii != 1000 and True:
        continue


    Papri = acls.canonical_shooting_angle(Xe,Xr)
    #consititution du sous dico , aject ajonction préalable
    #d'infos supplémentaire
    tup = tup + (Papri,)
    resultdic[tup] = dict()
    
    ### SNELL DESCARTES
    print("SNELL DESCARTES")
    results_sd = rt.raytrace_seek(Xe,Xr,Zextract,Cextract,verbose=0,fulloutput=0)
    TIMER = timeit.timeit("rt.raytrace_seek(Xe,Xr,Zextract,Cextract,verbose=0,fulloutput=0)",
    "from __main__ import Xe,Xr,Zextract,Cextract,rt",number=N)    
    resultdic[tup]['inv_eiko'] = results_sd + (TIMER/Nf,)
    
    ### EIKONAL
    print("EIKONAL")
    try:
        results_eiko = acls.raytrace_ODE_seek(Xe,Xr,ssf_opera,
                                              Papri,h,resotype,1,
                                              adaptative_step=adap,
                                              exectime_out=True)
        #TIMER = timeit.timeit("acls.raytrace_ODE_seek(Xe,Xr,ssf_opera,Papri,h,resotype,adaptative_step=adap)",
        #"from __main__ import Xe,Xr,ssf_opera,Papri,h,resotype,adap,acls",number=N)
        resultdic[tup]['inv_eiko'] = results_eiko # + (TIMER/Nf,)
    except:
        print(e)
        print("ERR : EIKO CRASH")
        resultdic[tup]['inv_eiko'] = np.nan

    ### EQUIV
    print("EQUIVALENT")
    paramtuple    = acls.equiv_get_param(Zextract,Cextract,Xr[-1])
    results_equiv = acls.equiv_inverse(paramtuple,Xe,Xr)
    TIMER = timeit.timeit("acls.equiv_inverse(paramtuple,Xe,Xr)",
    "from __main__ import paramtuple,Xe,Xr,acls",number=N)
    resultdic[tup]['inv_equiv'] = results_equiv + (TIMER/Nf,)    

outdir = "/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/scripts/benchmark_rt/results/"
gf.pickle_saver(resultslis,outdir,gf.get_timestamp() + 'benchmark_rt_inverse')

