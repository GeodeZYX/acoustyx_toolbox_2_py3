# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 18:22:29 2015

Falbication des shoots acoustiques + la trajectoire
 du bateau dans un ENVIRONNEMENT SIMPLE

@author: psakicki
"""

import glob
import acouclass as acls
import numpy as np
import raytrace as rt

from megalib import *

path = "/home/pierre/Documents/CODES/acoustyx_toolbox/working/data_2/"
outpath = "/home/pierre/Documents/CODES/acoustyx_toolbox/working/data_4/"

#path = '/home/psakicki/THESE/bd_CTD/NEW_DATA_2015/MOVE/MOVE3_epoch/data_2'
outpath='/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox/working/data_20/'


# === IMPORT DU SSP ===
#Z = np.loadtxt(path + 'Zout2000.txt')
#C = np.loadtxt(path + 'Cout2000.txt')

# Import d'un SSP réaliste, matlab-like trouvé dans un coin 
# ne doit plus être utilisé en prod 
import scipy.io as ioo
mat = ioo.loadmat('/home/psakicki/THESE/CODES/CodePython/GPShApy/SSPs/SSPGwada1A.mat')
Z = np.squeeze(mat['Prof'].T)
C = np.squeeze(mat['Vit'].T)

# FABRICATION DE LA TRAJECTOIRE
XYZ , E , _ = acls.fabriq_traject_droite(2000,2000,0,0,0,3,5000,45,vit=800,epoch_init=0,plot=0)

theta_stk = []
R_stk = []
T_stk = []

PXP = np.array([-0,-0,4000])
PXP = np.array([-2500,-2500,4000])


args_lis = []
pool = mp.Pool(processes=4)

xyzlis,PXPlis,Zlis,Clis = [],[],[],[]

for i in range(XYZ.shape[0]):
    xyz = XYZ[i,:]
    args_lis.append((xyz,PXP,Z,C,0,88,False,False))

results = [pool.apply(rt.raytrace_seek, args=x) for x in args_lis]  

theta_stk = [e[0] for e in results]
R_stk     = [e[1] for e in results]
T_stk     = [e[2] for e in results]

Tarr = np.expand_dims(T_stk,1)
Thetaarr = np.expand_dims(theta_stk,1)
Rarr = np.expand_dims(R_stk,1)

# noising
# the ships coordinate
sigma_x = 0.1
sigma_y = 0.1
sigma_z = 0.1
sigma_t = 10**-6 

XYZnoise = XYZ + np.random.randn(*XYZ.shape) * np.repeat([[sigma_x,sigma_y,sigma_z]],XYZ.shape[0],axis=0)
Tnoise = Tarr + np.random.randn(*Tarr.shape) * sigma_t

Mout = np.hstack((XYZnoise,Tnoise,XYZ,Tarr,Rarr,Thetaarr))

np.savetxt(outpath + 'Mout2000.txt',Mout)
np.savetxt(outpath + 'Zout2000.txt',Z)
np.savetxt(outpath + 'Cout2000.txt',C)
