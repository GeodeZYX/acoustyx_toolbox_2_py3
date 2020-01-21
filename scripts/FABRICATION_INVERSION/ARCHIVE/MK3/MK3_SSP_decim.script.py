# -*- coding: utf-8 -*-
"""
Created on Thu Oct  1 17:54:41 2015

@author: psakicki
"""

from megalib import *

filis = ["/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/exemple/SSP/SSP_NOAA_dep5645_20030608000000",
"/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/exemple/SSP/SSP_NOAA_dep5781_20050827000000"]


for f in filis:    
    ZC = np.loadtxt(f)
    comm = [l[1:] for l in open(f) if l[0] == '#' and not 'step' in l]
    Z = ZC[:,0]
    C = ZC[:,1]
    Zlight , Clight = ssp.SSP_light(Z,C)
    ZClight = np.column_stack((Zlight,Clight))
    comm.append('decimated ' + str(genefun.get_timestamp(0)) + '\n')
    commfin = ''.join(comm)
    np.savetxt(f+'_decim_'+genefun.get_timestamp(), ZClight,header=commfin)
    