# -*- coding: utf-8 -*-
"""
Created on Sat Jun 25 16:12:40 2016

@author: psakicki

cette fonction est l'implémentation de Franchi1973a
"""

from netCDF4 import Dataset
import scipy.interpolate as interpp
from megalib import *

path = '/home/psakicki/Téléchargements/OS_MOVE3_11_D_CURRENTMETER.nc'
D = Dataset(path)

np.squeeze(np.array(D.variables['CSPD'][:]))
np.squeeze(np.array(D.variables['PRES'][:]))



V =  [ 0.241,  0.167]
P =  [1373. , 3123.]
I = interpp.interp1d(P,V,fill_value='extrapolate')

I(0)


d   = 1.64 * 10**-3
rho = 1.02
b  = -7.75 * 10**-2

vs = 1.45

a = -(b - 1.5848 * (10 ** -2) * rho ) * d * vs


d   = 1.64 * 10**-3 #fixe
rho = 1.023
b   =   -0.05559274  #-7.75 * 10**-2
vs  = 0.24
a = -(b - 1.5848 * (10 ** -2) * rho ) * d * vs

V =  [0.09,  0]
P =  [0 , 5000]
I = interpp.interp1d(P,V,fill_value='extrapolate')
I(5000)

ZC = np.loadtxt("/home/psakicki/Documents/CODES/acoustyx_toolbox_2/exemple/SSP/SSP_NOAA_dep5645_20030608000000")

ssp.munk_pseudo_coef_calc(ZC[:,0],ZC[:,1])

import seawater

seawater.eos80.dens0(36.0409,27.6246)

#acls.SSP_2_bilin_grads(ZC[:,0],ZC[:,1],4000)





