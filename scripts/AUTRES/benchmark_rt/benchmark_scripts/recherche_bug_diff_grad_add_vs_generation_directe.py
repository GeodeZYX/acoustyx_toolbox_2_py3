#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 15:04:15 2017

@author: psakicki
"""

# -*- coding: utf-8 -*-
"""
Created on Sat Feb 20 14:15:01 2016

@author: psakicki

Le dit bug vient du choix du milieu
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
#prm_Udic_path = dir_prefix + "/acoustyx_toolbox_2/exemple/Udic/UdicFINAL_atal.pik"

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

SSF_Grad_off.make_a_SSP(0,0,0)


SSF_Grad_on_added_over_Grad_off = copy.deepcopy(SSF_Grad_off)

Tensor_grad = SSF_Grad_on_added_over_Grad_off.add_a_gradient_from_vectors(0,Grad_realistik,
                                                                          add_direct=True)
SSF_Grad_on_added_over_Grad_off.Cgrad

SSF_Grad_on.C[0] - SSF_Grad_on_added_over_Grad_off.C[0]


D1 = SSF_Grad_on_added_over_Grad_off.C - SSF_Grad_off.C
D2 = SSF_Grad_on.C                     - SSF_Grad_off.C


A = SSF_Grad_on.C[:,:,1]

np.diff(A,axis=0)

B = SSF_Grad_on_added_over_Grad_off.C[:,:,1]

np.diff(B,axis=0)
