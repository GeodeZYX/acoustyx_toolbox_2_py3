#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 14:35:20 2020

@author: psakicki
"""

#### Import star style
from geodezyx import *                   # Import the GeodeZYX modules
from geodezyx.externlib import *         # Import the external modules
from geodezyx.megalib.megalib import *   # Import the legacy modules names

import raytrace as rt
import SSP as ssp



## Immersion: statique mediane
depthdic = dict()

depthdic[1] = 37.57
depthdic[2] = 40.60
depthdic[3] = 37.84
depthdic[4] = 39.38

p_obs="/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/01_/PAMELI_BREST_vJupyter_1/PAMELI_BREST_vJupyter_1.csv"

DFbig = pd.read_csv(p_obs,index_col=0)

cols_posi = ['N_AHD_emi', 'E_AHD_emi', 'D_AHD_emi',
             'N_AHD_rec', 'E_AHD_rec', 'D_AHD_rec']

cols_posi_emi = cols_posi[0:3]
cols_posi_rec = cols_posi[3:6]

Z,C = ssp.munk(6000)

OPTIMstk = []

Xrec_apri_all = np.array([[24.31494786,  -7.20555817, 40.07],
                          [-26.20205757, -7.53913662, 40.66],
                          [2.23923942,   22.42808402, 39.68]])


for idbea in np.sort(DFbig['ID_BEA'].unique()):

    print("BEACON",idbea)
    
    Xrec_apri = Xrec_apri_all[idbea-1]
    
    
    def fct_inter_lsq(Xrec_apri_inp):
    
        DFbea = DFbig[DFbig['ID_BEA'] == idbea]
        
        Xemi = DFbea[cols_posi_emi].values
        Xrec = DFbea[cols_posi_rec].values
        TWTT = DFbea['TWTT'].values * 10**-6
       
        ZIP = list(zip(Xemi,Xrec,TWTT))
     
        ADTstk = []
        RESstk = []    
     
        for xemi,xrec,twtt in ZIP:
            
            if np.sum(np.isnan(xemi)) > 0:
                continue
            adt = rt.raytrace_seek((xemi,xrec),
                             Xrec_apri_inp,Z,C,
                             verbose=False,
                             fulloutput=False,
                             severe=True,
                             legacy=True)
                    
            ADTstk.append(adt)
            
            RESstk.append(adt[-1] - twtt)
            
        RES = np.array(RESstk)
        
        print("SUM RES",np.sum(np.square(RES)))
        
        return RES
    
    OPTIM = scipy.optimize.least_squares(fct_inter_lsq,Xrec_apri)
    OPTIMstk.append(OPTIM)
    
    
    
            
            
    
    
        


