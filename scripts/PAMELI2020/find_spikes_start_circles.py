#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 23:41:27 2020

@author: psakicki
"""

#### Import star style
from geodezyx import *                   # Import the GeodeZYX modules
from geodezyx.externlib import *         # Import the external modules
from geodezyx.megalib.megalib import *   # Import the legacy modules names


p1 = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/07_/PAM_BST_v211_m2507-1158_m3_d04_bea1_RGF93_RTKLIB/PAM_BST_v211_m2507-1158_m3_d04_bea1_RGF93_RTKLIB.csv"
p3 = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/07_/PAM_BST_v221_m2507-1109_m5_d02_bea2_RGF93_RTKLIB/PAM_BST_v221_m2507-1109_m5_d02_bea2_RGF93_RTKLIB.csv"
p2 = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/07_/PAM_BST_v231_m2507-1245_m4_d06_bea3_RGF93_RTKLIB/PAM_BST_v231_m2507-1245_m4_d06_bea3_RGF93_RTKLIB.csv"

P = [p1,p2,p3]

OUTDICT = dict()

for p in P:
    DF = pd.read_csv(p)
    
    bea = int(p.split("bea")[1][0])
    
    fig,ax = plt.subplots(2,1)
    
    T = conv.strdate2dt(DF["date_rec"])
    
    ax[0].plot(T,DF["E_AHD_rec"],"x")
    ax[1].plot(T,DF["N_AHD_rec"],"x")
    
    N = DF["N_AHD_rec"].values
    E = DF["E_AHD_rec"].values
    
    A = np.rad2deg(np.arctan2(E - np.nanmean(E),N - np.nanmean(N)))
    
    Tp = conv.dt2posix(T)
    
    plt.figure()
    
    plt.plot(Tp,A,".")
    
    I,_ = scipy.signal.find_peaks(A,height=None,distance=30,prominence=10)
    plt.plot(Tp[I],A[I],"+")
    
    START_raw = list(I / len(A))
    
    LENGTH = np.array(np.diff([0] + START_raw + [1]))
    START = np.array([0] + START_raw)
    END = START + LENGTH
    
    A = np.column_stack((START,LENGTH,END))
    
    OUTDICT[bea] = A
    
    print(len(A))
    
    utils.chunkIt(A,5)
    
utils.pickle_saver(OUTDICT,"/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING","spike_dict")
    
    
