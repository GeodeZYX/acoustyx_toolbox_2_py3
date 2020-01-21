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
import genefun
import timeit
import pickle

import multiprocessing as mp

from megalib import *

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

# ====================== FIN DU SSP ======================

# ======= EXTRACTION DU GRADIENT pseudo munkizé D'UN Udic ==========


if genefun.get_computer_name() == 'calipso':
    dir_prefix = '/home/psakicki/THESE/CODES/CodePython/'
elif genefun.get_computer_name() == 'brisote':
    dir_prefix = '/home/psakicki/Documents/CODES/'
elif genefun.get_computer_name() == 'psakicki-MS-16Y1':
    dir_prefix = '/home/adminuser/Documents/CODES/python2/'

prm_Udic_path = dir_prefix + "//acoustyx_toolbox_2/exemple/input_data_for_simulation/Udic/UdicFINAL_atal.pik"
prm_Udic_path = dir_prefix + "//acoustyx_toolbox_2/exemple/input_data_for_simulation/Udic/UdicFINAL_meteor.pik"

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

if 1:     
    export_ssf = SSF_Grad_on_added_over_Grad_off
    export_path = '/home/adminuser/aaa_FOURBI/'
    Gpath = export_path + 'Gradtensor' 
    Zpath = export_path + 'Ztensor' 
    Xpath = export_path + 'Xtensor' 
    Ypath = export_path + 'Ytensor' 
    Cpath = export_path + 'Ctensor' 

    np.save(Cpath,export_ssf.C)
    np.save(Xpath,export_ssf.X)
    np.save(Ypath,export_ssf.Y)
    np.save(Zpath,export_ssf.Z)
    np.save(Gpath,Tensor_grad)






# ======= FIN DE L'EXTRACTION DU GRADIENT pseudo munkizé D'UN Udic ==========

# ======= Préparation du Mode munk_n_add ==========

ssf_munk_non_grad = copy.deepcopy(ssf_opera)
ssf_munk_grad     = copy.deepcopy(ssf_opera)

Tensor_grad       = ssf_munk_grad.add_a_gradient_from_vectors(0,Grad_realistik,
                                          add_direct=True)

# ======= Fin de Préparation du Mode munk_n_add ==========




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


hlis        = [10]
restylis    = ['rkck']
zmaxlis     = [500, 1000, 2000 , 3000] 
shootanglis = tuple(np.arange(20,81,20)) + (25,75)
adaplis          = (True,)
z_max_smooth_lis = reversed([ 50 ,100 ,150, 200 ,250,300,400,500 ])
smoothstyle_lis  = ('lin','log')
x_grad_lis       = reversed([0,10**-6,5*10**-6,10**-5,5 * 10**-5,10**-4])
y_grad_lis       = [0]

if not mode == 'synthetik':
    print("not synthetik gradient, so z_max_smooth_lis, smoothstyle_lis, x/y_grad_lis useless")
    smoothstyle_lis  = (None,)
    x_grad_lis       = [None]
    y_grad_lis       = [None]
    z_max_smooth_lis = [None]
    suffix = suffix + '_' + prm_Udic_path.split('_')[-1].split('.')[0]

Xe = np.array([0,0,0])

N = 5
Nf = float(N)
#Papri = acls.canonical_shooting_angle(Xe,Xr)
resultdic = {}

varargslis = list(itertools.product(hlis,restylis,zmaxlis,shootanglis,adaplis,
                               z_max_smooth_lis,smoothstyle_lis,
                               x_grad_lis,y_grad_lis))

with_eiko = 1

#for iii , tup in enumerate(varargslis):    

def wrapperfct(tup):
    h , resty   , zmax   , shootang , adap , z_max_smooth ,  \
    smoothstyle , x_grad , y_grad = tup
    
    print("exp",  tup)
    Papri = (-shootang,0)
    
    resultdic[tup] = dict()
    
    
    
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
  
    
    if 0:     
        export_path = '/home/adminuser/aaa_FOURBI/'
        Cpath = export_path + 'Ctensor' 
        Xpath = export_path + 'Xtensor' 
        Ypath = export_path + 'Ytensor' 
        Zpath = export_path + 'Ztensor' 
        
        np.save(Cpath,ssf_opera_grad.C)
        np.save(Xpath,ssf_opera_grad.X)
        np.save(Ypath,ssf_opera_grad.Y)
        np.save(Zpath,ssf_opera_grad.Z)


    
    
    
    
    ### SNELL DESCARTES
    print("SNELL DESCARTES")
    sol_dir_SD = rt.raytrace_SD1_frontend(Zextract,Cextract,Papri[0],
                                          zmax,cumsum=0,out_path_length=True,
                                          shift_90deg = True)
    
    resultdic[tup]['SD'] = (sol_dir_SD , np.nan)
    
    ### EIKONAL
    print("EIKONAL")
    if with_eiko:
        sol_dir_eiko = acls.raytrace_ODE_2or3d(ssf_opera_grad,Xe,Papri,sol_dir_SD[-1],
                                               h=h,resotype=resty,
                                               return_mode='short',
                                               adaptative_step=adap,
                                               verbose=0)

        sol_dir_eiko_proto = acls.raytrace_ODE_2or3d(ssf_opera_non_grad,Xe,Papri,sol_dir_SD[-1],
                                               h=h,resotype=resty,
                                               return_mode='short',
                                               adaptative_step=adap,
                                               verbose=0)

        
        resultdic[tup]['eiko'] = (sol_dir_eiko , sol_dir_eiko_proto)
        print("result SD  " , resultdic[tup]['SD'])
        print("result eiko" , resultdic[tup]['eiko'])
        
    ### EQUIV
    print("EQUIVALENT")
    s2       = geok.pythagore(sol_dir_SD[0],zmax)
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
                                                

    resultdic[tup]['equiv_z' ] = (sol_dir_equiv_z,np.nan )
    resultdic[tup]['equiv_s' ] = (sol_dir_equiv_s,np.nan)
    resultdic[tup]['equiv_s2'] = (sol_dir_equiv_s2,np.nan)
    resultdic[tup]['equiv_x' ] = (sol_dir_equiv_x,np.nan)

    resultdic[tup]['equiv_z_anggeom' ] = (sol_dir_equiv_z_anggeom,0)
    resultdic[tup]['equiv_s_anggeom' ] = (sol_dir_equiv_s_anggeom,0)
    resultdic[tup]['equiv_s2_anggeom'] = (sol_dir_equiv_s2_anggeom,0)
    resultdic[tup]['equiv_x_anggeom' ] = (sol_dir_equiv_x_anggeom,0)
    
    return resultdic

pool = mp.Pool(processes=2)
#varargslis = varargslis[:1]

R  = [pool.apply_async(wrapperfct, args=(ar,)) for ar in varargslis]
R2 = [e.get() for e in R]  

finresultdic = dict()

for D in R2:
    for k,it in D.items():
        finresultdic[k] = it



if genefun.get_computer_name() == 'calipso':
    outdir  = "/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/scripts/benchmark_rt/results_REBOOT1703/"
elif genefun.get_computer_name() == 'brisote':
    outdir  = '/home/psakicki/Documents/CODES/acoustyx_toolbox_2/scripts/benchmark_rt/results_REBOOT1703'

gf.create_dir(outdir)
outpath = gf.pickle_saver(finresultdic,outdir,
                gf.get_timestamp() + '_benmk_direct_grd' + suffix)
                
print("results in" , outpath)