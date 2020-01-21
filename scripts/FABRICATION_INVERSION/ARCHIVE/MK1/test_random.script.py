# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 11:07:31 2015

@author: psakicki
"""

import acouclass as acls
import raytrace as rt
import numpy as np
import scipy.optimize as optimize


Zbazik = np.arange(0,5000,1)
Cbazik1 = 1500 + 20 * np.cos(Zbazik  * (0.001))
Cbazik2 = 1500 + 20 * np.sin(Zbazik  * (0.001))

SSPbazik1 = acls.SSP(Zbazik,Cbazik1,e=0)
SSPbazik2 = acls.SSP(Zbazik,Cbazik2,e=10)
zm,cm = rt.munk(6000)

SSPmunk=acls.SSP(zm,cm,e=0)

SSPbazik1.plot()
SSPbazik2.plot()

xmin,xmax,ymin,ymax,xstep,ystep = -2500 , 2500 , -2500 , 2500 , 500 , 500

ssf1 = acls.make_SSF3D_from_SSP(SSPbazik1,xmin,xmax,ymin,ymax,xstep,ystep)
ssf2 = acls.make_SSF3D_from_SSP(SSPbazik2,xmin,xmax,ymin,ymax,xstep,ystep)
ssf_munk = acls.make_SSF3D_from_SSP(SSPmunk)

#fig = mlab.figure()
#ssf1.plot()
#time.sleep(2)
##mlab.figure()
#mlab.clf()
#ssf2.plot()

#ssf4d1 = acls.SSF4D([ssf2])
#ssf4d1.add_ssf_in_ssf3d_list(ssf1)
#ssf4d1.ssf3d_list

#mlab.figure()
#ssf_out.plot()
#mlab.figure()
#ssf2.plot()

Zextract , Cextract = ssf_munk.make_a_SSP()

Z1 , C1  = ssf1.get_flatten_uniq('ZC')

def random_noise_Z_C(Z,C,zup,zdown,sigma=0.1,randgen=0):
    if randgen == 0:
        randgen = np.random.randint(1,999999)
    R = np.random.RandomState(randgen)
    Bools = (zup <= Z1) * (Z1 < zdown)
    dC = R.randn(*C1.shape) * sigma * Bools
    return dC , randgen
    

    
def random_noise_SSF(ssfin,z_interval_lis,sigma_lis,randgen_lis=[]):
    if randgen_lis == []:
        randgen_lis = list(np.zeros(len(z_interval_lis)))
    if not len(z_interval_lis) == len(sigma_lis) == len(randgen_lis):
        raise Exception('not len(zup_lis) == len(zdown_lis) == len(sigma_lis)')

    Z , C  = ssfin.get_flatten_uniq('ZC')
    dC_final = np.zeros(len(Z))
    randgen_lis_final = []
    for zinter, sig , randgen in zip(z_interval_lis,sigma_lis,randgen_lis):
        zup , zdown  = zinter
        dC , randgen = random_noise_Z_C(Z,C,zup,zdown,sig,randgen)
        dC_final = dC + dC_final
        
        randgen_lis_final.append(randgen)
        
    return dC_final , randgen_lis_final

z_inter = [(0,800) , (1600,9999)]
sigma_lis = [0.1 , 0.5]

dc ,rangen = random_noise_SSF(ssf1,z_inter,sigma_lis)

plt.figure()
plt.plot(dc)

