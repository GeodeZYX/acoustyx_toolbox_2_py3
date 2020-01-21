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
outpath='/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox/working/data_1/'

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

XYZ , E , _ = acls.fabriq_traject_droite(2000,2000,0,0,0,2,1000,80,vit=800,epoch_init=0)

np.cumsum(np.linalg.norm(np.diff(XYZ,axis=0),axis=1)) / 2000

theta_stk = []
R_stk = []
T_stk = []


PXP = np.array([-2500,-2500,4000])

args_lis = []
pool = mp.Pool(processes=4)

print("aaa")

xyzlis,PXPlis,Zlis,Clis = [],[],[],[]

##préparation des Z et C
#arg_e_lis = []
#for i in range(XYZ.shape[0]):
#    arg_e_lis.append(E[i])
 
#results_Z = [pool.apply(acls.ZfromSSP_wrapper, args=(ssf4d,e)) for e in arg_e_lis]
#results_C = [pool.apply(acls.CfromSSP_wrapper, args=(ssf4d,e)) for e in arg_e_lis]
#Zlis = [z.get() for z in results_Z]
#Clis = [c.get() for c in results_C]

for i in range(XYZ.shape[0]):
    print(i)
    xyz = XYZ[i] 
    e = E[i]
#    Zlis = results_Z[i]
#    Clis = results_C[i]
    Z = ssf4d.make_a_SSP(e=e).Z
    C = ssf4d.make_a_SSP(e=e).C
    args_lis.append((xyz,PXP,Z,C,1,89,False,False))
    
#    theta , R , T = rt.raytrace_seek(xyz,PXP,Z,C,1,89,True,False)
#    theta_stk.append(theta)
#    R_stk.append(R)
#    T_stk.append(T)

results = [pool.apply_async(rt.raytrace_seek, args=x) for x in args_lis]  
#results = pool.map_async(raytrace_seek_wrap,args_lis)

theta_stk = [e.get()[0] for e in results]
R_stk     = [e.get()[1] for e in results]
T_stk     = [e.get()[2] for e in results]

Tarr = np.expand_dims(T_stk,1)
Thetaarr = np.expand_dims(theta_stk,1)
Rarr = np.expand_dims(R_stk,1)

# noising
# the ships coordinate
sigma_x = 0.03
sigma_y = 0.03
sigma_z = 0.06
sigma_t = 10**-3

XYZnoise = XYZ + np.random.rand(*XYZ.shape) * np.repeat([[sigma_x,sigma_y,sigma_z]],XYZ.shape[0],axis=0)
Tnoise = Tarr + np.random.rand(*Tarr.shape) * sigma_t

Mout = np.hstack((XYZnoise,Tnoise,XYZ,Tarr,Rarr,Thetaarr))

np.savetxt(outpath + 'Mout2000.txt',Mout)
Zout = ssf4d.make_a_SSP(e=0).Z
Cout = ssf4d.make_a_SSP(e=0).C
np.savetxt(outpath + 'Zout2000.txt',Zout)
np.savetxt(outpath + 'Cout2000.txt',Cout)

#print 'fin'