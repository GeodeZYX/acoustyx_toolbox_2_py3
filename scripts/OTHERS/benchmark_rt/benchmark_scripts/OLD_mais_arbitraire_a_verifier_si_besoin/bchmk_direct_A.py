# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 11:01:32 2016

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

# ==================== DEBUT DU SSP ====================
Zbazik  = np.arange(0,5000,1)
Cbazik1 = 1500 + 20  * np.cos(Zbazik  * (0.001))
Cbazik2 = 1500 + 20  * np.sin(Zbazik  * (0.001))
Cbazik3 = 1500 + 500 * np.sin(Zbazik  * (0.01))

SSPbazik1 = acls.SSP(Zbazik,Cbazik1,e= 0)
SSPbazik2 = acls.SSP(Zbazik,Cbazik2,e=10)
SSPbazik3 = acls.SSP(Zbazik,Cbazik3,e=10)

zmaxtop = 4000
zm,cm   = rt.munk(zmaxtop,10**0)
SSPmunk = acls.SSP(zm,cm,e=0)

#SSPbazik1.plot()
ssf_munk   = acls.make_SSF3D_from_SSP(SSPmunk)
ssf_bazik2 = acls.make_SSF3D_from_SSP(SSPbazik2)
ssf_bazik3 = acls.make_SSF3D_from_SSP(SSPbazik3)

ssf_opera = ssf_bazik2
ssf_opera = ssf_munk

Zextract , Cextract = ssf_opera.make_a_SSP(obj=False)

if 0:
    axiss0       =  1
    axiss1       =  1
    dmin         = -1
    dmax         =  1

    z_max_smooth = 250
    Tensor_grad  = ssf_opera.add_a_gradient(axiss0,(dmin,dmax),
                                          z_max_smooth=z_max_smooth)
    Tensor_grad= ssf_opera.add_a_gradient(axiss1,(dmin,dmax),
                                          z_max_smooth=z_max_smooth)
  
#ssf_opera.plot(1)
# ====================== FIN DU SSP ======================

##### BLOC DES EXPERIENCES

hlis        = [10]
restylis    = ['rkdp']
zmaxlis     = [500, 1000, 2000, 3000]
shootanglis = np.arange(20,81,20)
h = hlis[0]
adaplis = (True,)

hlis        = [1 , 10 , 100 , 500 , 1000 , 2000 , 3000  ]
restylis    = ["rk2",'rk4','rkf45','rkck','rkdp','euler']
zmaxlis     = [500, 1000, 2000 , 3000] 
shootanglis = np.arange(20,81,20)
h = hlis[0]
adaplis = (True , False)

hlis        = [1 , 2, 3, 5, 10, 20, 30, 40, 50,  75, 100]
restylis    = ["rk2",'rk4','rkf45','rkck','rkdp','euler']
zmaxlis     = [500, 1000, 2000 , 3000] 
shootanglis = tuple(np.arange(20,81,20)) + (25,75)
h = hlis[0]
adaplis = (True,)



Xe = np.array([0,0,0])

N = 5
Nf = float(N)
#Papri = acls.canonical_shooting_angle(Xe,Xr)
resultdic = {}

varargslis = itertools.product(hlis,restylis,zmaxlis,shootanglis,adaplis)

with_eiko = 1

for iii , tup in enumerate(varargslis):
    h , resty , zmax , shootang , adap = tup
    print("exp",iii , tup)
    Papri = (-shootang,0)
    
    resultdic[tup] = dict()
    
    ### SNELL DESCARTES
    print("SNELL DESCARTES")
    sol_dir_SD = rt.raytrace_SD1_frontend(Zextract,Cextract,Papri[0],
                                          zmax,cumsum=0,out_path_length=True,
                                          shift_90deg = True)
                                          
    sol_dir_SD_TIMER = timeit.timeit("rt.raytrace_SD1_frontend(Zextract,Cextract,Papri[0], \
                                          zmax,cumsum=0,out_path_length=True, \
                                          shift_90deg = True)",
    "from __main__ import Zextract,Cextract,Papri,zmax,rt",number=N)
    
    resultdic[tup]['SD'] = (sol_dir_SD , sol_dir_SD_TIMER / N)
    
    ### EIKONAL
    print("EIKONAL")
    if with_eiko:
        sol_dir_eiko = acls.raytrace_ODE_2or3d(ssf_opera,Xe,Papri,sol_dir_SD[-1],
                                               h=h,resotype=resty,
                                               return_mode='short',
                                               adaptative_step=adap,
                                               verbose=0)  
    
        sol_dir_eiko_TIMER = timeit.timeit("acls.raytrace_ODE_2or3d(ssf_opera,Xe, \
                                               Papri,sol_dir_SD[-1], \
                                               h=h,resotype=resty,   \
                                               return_mode='short',  \
                                               adaptative_step=adap, \
                                               verbose=0)",
        "from __main__ import ssf_opera,Xe,Papri,sol_dir_SD,h,resty,adap,acls",number=N)
        
        resultdic[tup]['eiko'] = (sol_dir_eiko , sol_dir_eiko_TIMER / Nf)
        print("result SD  " , resultdic[tup]['SD'])
        print("result eiko" , resultdic[tup]['eiko'])
        
    ### EQUIV
    print("EQUIVALENT")
    s2 = geok.pythagore(sol_dir_SD[0],zmax)
    anggeom  = np.rad2deg(np.arctan2(zmax,sol_dir_SD[0]))

    a_4_dir          = Papri[0]                                                              
    equiv_param      = acls.equiv_get_param(Zextract,Cextract,zmax)
    sol_dir_equiv_z  = acls.equiv_direct(equiv_param,Xe,a_4_dir,zmax,'z')

    sol_dir_equiv_s  = acls.equiv_direct(equiv_param,Xe,a_4_dir,
                                         sol_dir_SD[-1],'s')

    sol_dir_equiv_s2 = acls.equiv_direct(equiv_param,Xe,a_4_dir,s2,'s')
    sol_dir_equiv_x  = acls.equiv_direct(equiv_param,Xe,a_4_dir,
                                         sol_dir_SD[0],'x')

    sol_dir_equiv_z_anggeom  = acls.equiv_direct(equiv_param,Xe,-anggeom,
                                                zmax,'z')
    sol_dir_equiv_s_anggeom  = acls.equiv_direct(equiv_param,Xe,-anggeom,
                                                sol_dir_SD[-1],'s')
    sol_dir_equiv_s2_anggeom = acls.equiv_direct(equiv_param,Xe,-anggeom,
                                                 s2,'s')
    sol_dir_equiv_x_anggeom  = acls.equiv_direct(equiv_param,Xe,-anggeom,
                                                sol_dir_SD[0],'x')
                                                
    sol_dir_equiv_z_TIMER  = timeit.timeit("acls.equiv_direct(equiv_param,Xe,a_4_dir,zmax,'z')","from __main__ import equiv_param,Xe,a_4_dir,zmax,acls",number=N)
    sol_dir_equiv_s_TIMER  = timeit.timeit("acls.equiv_direct(equiv_param,Xe,a_4_dir,sol_dir_SD[-1],'s')","from __main__ import equiv_param,Xe,a_4_dir,sol_dir_SD,acls",number=N)
    sol_dir_equiv_s2_TIMER = timeit.timeit("acls.equiv_direct(equiv_param,Xe,a_4_dir,s2,'s')","from __main__ import equiv_param,Xe,a_4_dir,sol_dir_SD,s2,acls",number=N)
    sol_dir_equiv_x_TIMER  = timeit.timeit("acls.equiv_direct(equiv_param,Xe,a_4_dir,sol_dir_SD[0],'x')" ,"from __main__ import equiv_param,Xe,a_4_dir,sol_dir_SD,acls",number=N)
    
    resultdic[tup]['equiv_z' ] = (sol_dir_equiv_z,sol_dir_equiv_z_TIMER / Nf)
    resultdic[tup]['equiv_s' ] = (sol_dir_equiv_s,sol_dir_equiv_s_TIMER / Nf)
    resultdic[tup]['equiv_s2'] = (sol_dir_equiv_s2,sol_dir_equiv_s_TIMER / Nf)
    resultdic[tup]['equiv_x' ] = (sol_dir_equiv_x,sol_dir_equiv_x_TIMER / Nf)

    resultdic[tup]['equiv_z_anggeom' ] = (sol_dir_equiv_z_anggeom,0)
    resultdic[tup]['equiv_s_anggeom' ] = (sol_dir_equiv_s_anggeom,0)
    resultdic[tup]['equiv_s2_anggeom'] = (sol_dir_equiv_s2_anggeom,0)
    resultdic[tup]['equiv_x_anggeom' ] = (sol_dir_equiv_x_anggeom,0)
    
gf.pickle_saver(resultdic,"/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/scripts/benchmark_rt/results/",
                gf.get_timestamp() + '_benmk_direct_monoproc')
