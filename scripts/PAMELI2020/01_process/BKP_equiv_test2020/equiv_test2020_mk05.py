#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 18:14:03 2020

@author: psakicki

Test differences between theoretical and observed ranges
"""

#### Import star style
from geodezyx import *                   # Import the GeodeZYX modules
from geodezyx.externlib import *         # Import the external modules
from geodezyx.megalib.megalib import *   # Import the legacy modules names

p_obs1 = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/01_/PAMELI_BREST_vJupyter_6_4transducers_newRF_newDTime/PAMELI_BREST_vJupyter_6_4transducers_newRF_newDTime.csv"
p_obs2 = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/01_/PAMELI_BREST_vJupyter_5_4transducers_new_lever_arms_GINS_PPP_coords/PAMELI_BREST_vJupyter_5_4transducers_new_lever_arms.csv"

DFbig1 = pd.read_csv(p_obs1,index_col=0)
DFbig_orig1 = DFbig1.copy()
DFbig1 = DFbig1.dropna()
DFbig1['VALID'] = np.ones(len(DFbig1['TWTT'])).astype(bool)
DFbig1['date_rec'].apply(conv.string_date2dt)
DFbig1['TWTT_obs'] = DFbig1['TWTT'] * 10**-6


DFbig2 = pd.read_csv(p_obs2,index_col=0)
DFbig_orig2 = DFbig2.copy()
DFbig2 = DFbig2.dropna()
DFbig2['VALID'] = np.ones(len(DFbig2['TWTT'])).astype(bool)
DFbig2['date_rec'].apply(conv.string_date2dt)
DFbig2['TWTT_obs'] = DFbig2['TWTT'] * 10**-6


cols_posi = ['N_AHD_emi',
             'E_AHD_emi', 
             'D_AHD_emi',
             'N_AHD_rec',
             'E_AHD_rec',
             'D_AHD_rec']

'N_AHD_rec','E_AHD_rec','D_AHD_rec'

cols_posi_emi = ["date_emi"] +  cols_posi[0:3]
cols_posi_rec = ["date_emi"] +  cols_posi[3:6]


cols_posi_emi = cols_posi[0:3]
cols_posi_rec = cols_posi[3:6]

DFbig1[cols_posi_emi] - DFbig2[cols_posi_emi]
DFbig1[cols_posi_rec] - DFbig2[cols_posi_rec]



DFbig = DFbig1

Xbea_apri_all_init = np.array([[   -7.20555817, 24.31494786, 40.07],
                               [   -7.53913662,-26.20205757, 40.66],
                               [  22.42808402,  2.23923942,  39.68]])

Xbea_apri_all_init = np.array([[24.31494786,   -7.20555817, 40.07],
                               [-26.20205757,  -7.53913662, 40.66],
                               [2.23923942,    22.42808402, 39.68]])

Xbea_apri_all = Xbea_apri_all_init

ID_bea_list_orig = np.sort(DFbig['ID_BEA'].unique())
ID_bea_list = ID_bea_list_orig

c = 1520.4208032999682

def fct_obs_equiv(xbea_x_apri_in,
                  xbea_y_apri_in,
                  xbea_z_apri_in,
                  xemi_in,
                  xrec_in,
                  c_bea_in):
    
    xbea_apri_in =  np.array([xbea_x_apri_in,
                              xbea_y_apri_in,
                              xbea_z_apri_in])
    
    #d_emi = conv.dist(xbea_apri_in,xemi_in)
    #d_rec = conv.dist(xbea_apri_in,xrec_in)

    d_emi = np.linalg.norm(xbea_apri_in-xemi_in,axis=1)
    d_rec = np.linalg.norm(xbea_apri_in-xrec_in,axis=1)
    
    print("d_emi",d_emi)
    
    #t_emi = d_emi / c_bea_in
    #t_rec = d_rec / c_bea_in
    
    t_sum = (d_emi + d_rec) / c_bea_in
    
    return t_sum



for i_idbea, idbea in enumerate(ID_bea_list):
    xbea_apri = Xbea_apri_all[idbea-1]

    Dtheo_emi = np.linalg.norm(DFbig[cols_posi_emi].values - xbea_apri,axis=1)
    Dtheo_emi = Dtheo_emi[:,np.newaxis]
    Dtheo_rec = np.linalg.norm(DFbig[cols_posi_rec].values - xbea_apri,axis=1)
    Dtheo_rec = Dtheo_rec[:,np.newaxis]
    
    Dtheo = Dtheo_rec + Dtheo_emi
    
    Dobs = DFbig['TWTT_obs'] * c
    dD   = Dobs - np.squeeze(Dtheo)
    
    Ttheo = Dtheo / c
    Tobs = DFbig['TWTT_obs']
    
    np.squeeze(Ttheo) - Tobs
        
    fct_obs_equiv(xbea_apri[0],
                  xbea_apri[1],
                  xbea_apri[2],
                  DFbig[cols_posi_emi].values,
                  DFbig[cols_posi_rec].values,
                  c)
    
    print(dD / c)