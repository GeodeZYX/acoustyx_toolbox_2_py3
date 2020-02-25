# -*- coding: utf-8 -*-
"""
Created on Sat Feb 20 14:15:01 2016

@author: psakicki

TENTATIVE DE //ISATION

Est consideré comme prioritaire pour les phases de developpement/upgrade
par rapport à son frère direct C gradient

Le pas du SSP de munk influe grandement sur le temps de calcul
(ce qui explique la diff de tps de calcul avec un SSP realistik)
"""

import acouclass as acls
import raytrace as rt
import numpy as np
import scipy.optimize as optimize
#from megalib import *
import time
import copy
import itertools
import genefun
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
ssf_munk   = acls.make_SSF3D_from_SSP(SSPmunk,
                xmin=-5000,
                xmax= 5000,
                ymin=-5000,
                ymax= 5000)
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

# ======= EXTRACTION DU GRADIENT pseudo munkizé D'UN Udic ==========


if genefun.get_computer_name() == 'calipso':
    dir_prefix = '/home/psakicki/THESE/CODES/CodePython/'
elif genefun.get_computer_name() == 'brisote':
    dir_prefix = '/home/psakicki/Documents/CODES/'

prm_Udic_path = dir_prefix + "/acoustyx_toolbox_2/exemple/Udic/UdicFINAL_meteor.pik"
prm_Udic_path = dir_prefix + "/acoustyx_toolbox_2/exemple/Udic/UdicFINAL_atal.pik"

Udic = gf.pickle_loader(prm_Udic_path)
X4I , Z4I , CgridXZ4I = Udic['X'] , Udic['Z_X'] ,  Udic['CgridXZ']

SSF_Grad_on =  acls.make_SSF3D_from_SSP_n_distance((X4I,Z4I,CgridXZ4I),'X',
                               xmin=-5000,
                               xmax= 5000,
                               ymin=-5000,
                               ymax= 5000)

C_grad_realistik = SSF_Grad_on.Cgrad[0][0,0,:]
Z_grad_realistik = SSF_Grad_on.Z[0,0,:]
Grad_realistik   = (Z_grad_realistik , C_grad_realistik)

SSPctr = SSF_Grad_on.make_a_SSP(*SSF_Grad_on.center(median=True),obj=1)
SSF_Grad_off    = acls.make_SSF3D_from_SSP(SSPctr,
                xmin=-5000,
                xmax= 5000,
                ymin=-5000,
                ymax= 5000)

SSF_Grad_on_added_over_Grad_off = copy.deepcopy(SSF_Grad_off)

Tensor_grad = SSF_Grad_on_added_over_Grad_off.add_a_gradient_from_vectors(0,Grad_realistik,
                                                                          add_direct=True)
SSF_Grad_on_added_over_Grad_off.Cgrad

SSF_Grad_on.C[0] - SSF_Grad_on_added_over_Grad_off.C[0]


D1 = SSF_Grad_on_added_over_Grad_off.C - SSF_Grad_off.C
D2 = SSF_Grad_on.C                     - SSF_Grad_off.C

D1[:,:,1]  
D2[:,:,1]  

# ======= FIN DE L'EXTRACTION DU GRADIENT pseudo munkizé D'UN Udic ==========

# ======= Préparation du Mode munk_n_add ==========

ssf_munk_non_grad = copy.deepcopy(ssf_opera)
ssf_munk_grad     = copy.deepcopy(ssf_opera)

Tensor_grad       = ssf_munk_grad.add_a_gradient_from_vectors(0,Grad_realistik,
                                          add_direct=True)

# ======= Fin de Préparation du Mode munk_n_add ==========


Xe = np.array([0,0,0])

#hlis     = [1 , 10 , 100 , 500 , 1000 , 2000 , 3000  ]
#hlis     = [1 , 5, 10, 20, 30, 40, 50, 60, 70, 80, 90]
#hlis     = [1 , 2, 3, 5, 10, 20, 30, 40, 50,  75, 100]
#hlis     = [1 , 2, 3, 5, 10, 20, 30 , 50,  75, 100]
#
#restylis = ["rk2",'rk4','rkf45','rkck','rkdp','euler']

#### mode

# synthetik        : gradients artificiels comme dans le papier v2016
# realistik_simple : gradients issus d'un Udic, avec la methode 
#                    tel que dans les scripts de fabrikation GNSS/A
# realistik_add    : on prend les grands issus des Udic, MAIS
#                    on fabrique les 2 SFF avec la methode add
# munk_n_add       : on ajoute avec la methode add les gradients des Udic
#                    MAIS à un un profil de munk (au profil opera + exactement)

mode = 'munk_n_add'
mode = 'realistik_simple'
mode = 'realistik_add'

suffix = mode

# SHORT VERSION
hlis        = [10]
restylis    = ['rkck']
adaplis  = (True , )
xrlis    = [  3000 ] 
yrlis    = [  3000 ] 
yrlis    = [ 0 ] 
zrlis    = [ 3000 ] 
z_max_smooth_lis = [ 500 ]
smoothstyle_lis  = ('lin','log')
x_grad_lis       = [0]
y_grad_lis       = [0]


# NOMINAL
hlis        = [10]
restylis    = ['rkck']
adaplis  = (True , )
xrlis    = [ 500 , 1000 , 1500 , 2000 , 2500 , 3000 ] 
yrlis    = [ 500 , 1000 , 2000 , 3000 ] 
yrlis    = [ 0 ] 
zrlis    = [ 500 , 1000 , 1500 , 2000 , 2500 , 3000 ] 
z_max_smooth_lis = [ 50 ,100 ,150, 200 ,250,300,400,500 ]
smoothstyle_lis  = ('lin','log')
x_grad_lis       = [0,10**-6,5*10**-6,10**-5,5 * 10**-5,10**-4]
y_grad_lis       = [0]

#%%


if not mode == 'synthetik':
    print("not synthetik gradient, so z_max_smooth_lis, smoothstyle_lis, x/y_grad_lis useless")
    smoothstyle_lis  = (None,)
    x_grad_lis       = [None]
    y_grad_lis       = [None]
    z_max_smooth_lis = [None]
    suffix = suffix + '_' + prm_Udic_path.split('_')[-1].split('.')[0]

    

varargslis_proto = list(itertools.product(hlis,restylis,adaplis,
                                    xrlis,yrlis,zrlis,z_max_smooth_lis,
                                    smoothstyle_lis,x_grad_lis,y_grad_lis))

varargslis = []
for e in varargslis_proto:
    x,y = e[-3],e[-2]
    if y > x:
        continue
    else:
        varargslis.append((e,))
        
print(len(varargslis) , 'experiences')

resultdic = dict()

N  = 1
Nf = float(N)
with_equiv = False

ssf_opera_proto = ssf_opera

def wrapperfct(intup):
    tup = intup
    h , resotype , adap , xr , yr , zr , z_max_smooth , smoothstyle , x_grad , y_grad = tup
    
    print(tup)
    
    if mode == 'synthetik':
        # Application du gradient 
        ssf_opera_grad     = copy.deepcopy(ssf_opera_proto)
        ssf_opera_non_grad = copy.deepcopy(ssf_opera_proto)
    
        if smoothstyle == 'lin':
            smooth_style_fct = np.linspace
        elif smoothstyle == 'log':
            smooth_style_fct = np.logspace
            Tensor_grad_x = ssf_opera_grad.add_a_gradient(0,x_grad,z_max_smooth=z_max_smooth,
                                                  smooth_style_fct = smooth_style_fct)
            Tensor_grad_y = ssf_opera_grad.add_a_gradient(1,y_grad,z_max_smooth=z_max_smooth,
                                                  smooth_style_fct = smooth_style_fct)
    
    elif mode == 'realistik_simple':
        #Tensor_grad = ssf_opera.add_a_gradient_from_vectors(0,Grad_realistik,
        #                                          add_direct=True)
        ssf_opera_grad       = SSF_Grad_on
        ssf_opera_non_grad   = SSF_Grad_off
    elif mode == 'realistik_add':
        ssf_opera_grad       = SSF_Grad_on_added_over_Grad_off    
        ssf_opera_non_grad   = SSF_Grad_off
    elif mode == 'munk_n_add':
        ssf_opera_grad       = ssf_munk_grad    
        ssf_opera_non_grad   = ssf_munk_non_grad   

    Xr = np.array([xr,yr,zr])
    
    Papri = acls.canonical_shooting_angle(Xe,Xr)
    #consititution du sous dico , aject ajonction préalable
    #d'infos supplémentaire
    tup = tup + (Papri,)
    
    ### SNELL DESCARTES
    print("SNELL DESCARTES")
    strt = time.time()
    for i in range(N):
        results_sd = rt.raytrace_seek(Xe,Xr,Zextract,Cextract,
                                      verbose=0,fulloutput=0)

    TIMER_sd  = time.time() - strt
    
    _,_,_,Ssd = rt.raytrace_SD1_frontend(Zextract,Cextract,results_sd[0],Xr[2],
                                         out_path_length=True)
                                         
    results_sd = results_sd + (Ssd , TIMER_sd/Nf)
    ### EIKONAL
    print("EIKONAL")
    try:
        strt         = time.time()
        for i in range(N):
            results_eiko = acls.raytrace_ODE_seek(Xe,Xr,ssf_opera_grad,Papri,h,resotype,1,
                                              adaptative_step=adap,exectime_out=True)

        results_eiko_proto = acls.raytrace_ODE_seek(Xe,Xr,ssf_opera_non_grad,Papri,h,resotype,1,
                                          adaptative_step=adap,exectime_out=True)
        TIMER_eiko   = time.time() - strt
        results_eiko = results_eiko + results_eiko_proto + (TIMER_eiko/Nf,)
    except:
        print(e)
        print("ERR : EIKO CRASH")
        results_eiko = np.nan

    ### EQUIV
    if 0:
        print("EQUIVALENT")
        paramtuple    = acls.equiv_get_param(Zextract,Cextract,Xr[-1])
        strt         = time.time()
        for i in range(N):
            results_equiv = acls.equiv_inverse(paramtuple,Xe,Xr)
        TIMER_equiv   = time.time() - strt
        results_equiv = results_equiv + (TIMER_equiv/Nf,)
    else:
        print("EQUIV SKIPED")
        results_equiv = np.nan

    return tup , results_sd , results_eiko , results_equiv
  

def wrapperfct2(inptup):
    return gf.timeout(wrapperfct,(inptup,) , timeout_duration=100)
        

pool   = mp.Pool(processes=3)

varargslis = list(reversed(varargslis))
#varargslis = varargslis[:2]

R  = [pool.apply_async(wrapperfct2, args=ar) for ar in varargslis]
R2 = [e.get() for e in R]

outdir  = dir_prefix + "/acoustyx_toolbox_2/scripts/benchmark_rt/results_REBOOT1703/"
outpath = gf.pickle_saver(R2,outdir,gf.get_timestamp() + 'benmk_inverse_grd' + suffix)

print("results in" , outdir)