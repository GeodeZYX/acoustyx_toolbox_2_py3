# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 18:22:29 2015

@author: psakicki
"""

import glob
import acouclass as acls
import numpy as np
import raytrace as rt

from megalib import *

path='/home/psakicki/THESE/bd_CTD/NEW_DATA_2015/MOVE/MOVE3_epoch/SSP*'
outpath='/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox/working/data_2_4_PXP/'

filis = glob.glob(path)

SSPlis = []
for f in filis:
    SSPlis.append(acls.read_ssp_file(f))

SSPlis[0].C

SSPlis.sort(key=lambda x: x.t)
# fit d'une double expo pour avoir un ssp un peu plus sioux
SSPlis = [acls.fit_dbl_exp_SSP_objt(s) for s in SSPlis]
# ajout d'une valeur de s quelcoque à la surface
SSPlis = [acls.add_z0_mes_in_a_ssp(s,1537.,10.) for s in SSPlis]

tref = SSPlis[0].t

SSPlis[0].C

ssf3d_lis = []
for ssp in SSPlis:
    ssp.e = (ssp.t - tref).total_seconds()/3600.
    
    ssf3d = acls.make_SSF3D_from_SSP(ssp)
    ssf3d_lis.append(ssf3d)
    
    
ssf4d = acls.SSF4D(ssf3d_lis)
    
    
I = ssf4d.Interp_C

I((0,0,0,0.1))

ssf4d.get_SSF3D_fct_epoch(2)
A  = acls.get_SSF3D_from_SSF4D_fct_epoch(ssf4d,8)
A.make_a_SSP(0,0).plot()

plt.clf()
Z1 = ssf4d.make_a_SSP(e=20).Z
C1 = ssf4d.make_a_SSP(e=20).C

nbshoot = 100
XYZ , E , _ = acls.fabriq_traject_droite(2000,2000,0,0,0,2,nbshoot,80,vit=800,epoch_init=0)

np.cumsum(np.linalg.norm(np.diff(XYZ,axis=0),axis=1)) / 2000

PXP = np.array([-2500,-2500,4000])
PXP_coord_lis = [np.array([-2500,-2500,4000]),
                 np.array([2500,-2500,4000]) ,
                 np.array([-2500,2500,4000]) ,
                 np.array([2500,2500,4000])]

PXP = np.array([-2500,-2500,4000])

args_lis = []

# noising
# the ships coordinate
sigma_x = 0.03
sigma_y = 0.03
sigma_z = 0.06

XYZnoise = XYZ + np.random.rand(*XYZ.shape) * np.repeat([[sigma_x,sigma_y,sigma_z]],XYZ.shape[0],axis=0)

xyzlis,PXPlis,Zlis,Clis = [],[],[],[]



for i,e in enumerate(E):
    print(i)
    Zlis.append(ssf4d.make_a_SSP(e=e).Z)
    Clis.append(ssf4d.make_a_SSP(e=e).C)


for ipxp,PXP in enumerate(PXP_coord_lis):
    pool = mp.Pool(processes=4)
    print('PXP',ipxp)
    R_stk, theta_stk , T_stk = [],[],[]   
    
    args_lis = []
    for i in range(XYZ.shape[0]):
        
        print(i)
        xyz = XYZ[i] 
        e = E[i]
        Z = Zlis[i]
        C = Clis[i]

        args_lis.append((xyz,PXP,Z,C,1,89,False,False))
    
        #    theta , R , T = rt.raytrace_seek(xyz,PXP,Z,C,1,89,True,False)
        #    theta_stk.append(theta)
        #    R_stk.append(R)
        #    T_stk.append(T)
    
    results = [pool.apply_async(rt.raytrace_seek, args=x) for x in args_lis]  
    
    theta_stk = [e.get()[0] for e in results]
    R_stk     = [e.get()[1] for e in results]
    T_stk     = [e.get()[2] for e in results]

    Thetaarr = np.expand_dims(theta_stk,1)
    Tarr = np.expand_dims(T_stk,1)
    Rarr = np.expand_dims(R_stk,1)

    # noising the time
    sigma_t = 10**-3
    Tnoise = Tarr + np.random.rand(*Tarr.shape) * sigma_t
    # Noising the PXP
    sigma_xy_PXP = 10.
    sigma_z_PXP  = 15.
    PXP_noise = np.array([PXP[0] + np.random.rand() * sigma_xy_PXP,
                          PXP[1] + np.random.rand() * sigma_xy_PXP,
                          PXP[2] + np.random.rand() * sigma_z_PXP ])
                          
    PXP_clean_4_print = np.expand_dims(PXP,0).repeat(len(Tnoise),0)
    PXP_noise_4_print = np.expand_dims(PXP_noise,0).repeat(len(Tnoise),0)

    Mout = np.hstack((XYZnoise,
                      Tnoise,
                      XYZ,
                      Tarr,
                      Rarr,
                      Thetaarr,
                      PXP_noise_4_print,
                      PXP_clean_4_print))

    np.savetxt(outpath + 'Mout_'+str(nbshoot)+'_PXP'+str(ipxp),Mout)
        

Zout = ssf4d.make_a_SSP(e=0).Z
Cout = ssf4d.make_a_SSP(e=0).C
np.savetxt(outpath + 'Zout'+str(nbshoot),Zout)
np.savetxt(outpath + 'Cout'+str(nbshoot),Cout)