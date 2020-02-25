#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 24 17:43:45 2018

@author: psakicki
"""

from megalib import *

p="/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/1709_GNSSA_InSitu_EXPs/1808_reboot_GEODESEA_for_VB-HDR/REBOOT1808/classic/IUEM_LAPTOP-3Beacontracking/IUEM_LAPTOP-3Beacontracking__20180830_144020.smartV"
p="/home/psakicki/GFZ_WORK/RENDU/1709_CANOPUS_Cloture/report_canopus/1505-1/Detection15mai__20170727_115538__iter4.smartV"


A = np.loadtxt(p)

A[:,0]

DIVACOU  = [1,3,4]
GEODESEA = [3,5,7]

for i in range(3):
    AA = A[A[:,0] == i]
    T  = AA[:,1]
    Tdt = geok.posix2dt(T)
    D = AA[:,2]
    
    string = GEODESEA[i]
    string = DIVACOU[i]


        
    plt.plot(Tdt , D * 1500 * .5,".",label="AMT " + str(string),markersize=1)
    
plt.legend( markerscale=10)
plt.xlabel("Time")
plt.ylabel("Residual (m)")
plt.tight_layout()