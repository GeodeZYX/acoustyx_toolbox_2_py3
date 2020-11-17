#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 18:03:16 2020

@author: psakicki
"""

#### Import star style
from geodezyx import *                   # Import the GeodeZYX modules
from geodezyx.externlib import *         # Import the external modules
from geodezyx.megalib.megalib import *   # Import the legacy modules names





#F = "/home/psakicki/CODES/acoustyx_toolbox_2_py3/working/SIMU_SUMMER2020_1000_R50_nois1-0.0____/SIMU_SUMMER2020_1000_R50_nois1-0.0_____.PXP1.P.dat"

p = "/home/psakicki/CODES/acoustyx_toolbox_2_py3/working/SIMU_SUMMER2020_1000_R50_nois1-0.0____/"
p = "/home/psakicki/CODES/acoustyx_toolbox_2_py3/working/SIMU_SUMMER2020_3x333_x500_y500_ang0_nois1-0.0____"

L = utils.find_recursive(p,"*P.dat")

DFstk = []

for f in L:
    M  = np.loadtxt(f)
    DF = pd.DataFrame(M)
    
    rename_dict = dict()
    
    # E_emi_noise   date_emi
    # X_emi_noise   E_AHD_emi
    # Y_emi_noise   N_AHD_emi
    # Z_emi_noise   D_AHD_emi
    # T_emi_noise   NA
    # E_rec_noise   date_rec
    # X_rec_noise   E_AHD_rec
    # Y_rec_noise   N_AHD_rec
    # Z_rec_noise   D_AHD_rec
    # T_rec_noise   NA
    # TAT           NA
    # TWTT_noise    TWTT_obs
    # E_emi_clean
    # X_emi_clean
    # Y_emi_clean
    # Z_emi_clean
    # T_emi_clean
    # E_rec_clean
    # X_rec_clean
    # Y_rec_clean
    # Z_rec_clean
    # T_rec_clean
    # TWTT_clean
    # Noise_emi
    # Noise_rec
    
    rename_dict[0]='E_emi_noise'
    rename_dict[1]='X_emi_noise'
    rename_dict[2]='Y_emi_noise'
    rename_dict[3]='Z_emi_noise'
    rename_dict[4]='T_emi_noise'
    rename_dict[5]='E_rec_noise'
    rename_dict[6]='X_rec_noise'
    rename_dict[7]='Y_rec_noise'
    rename_dict[8]='Z_rec_noise'
    rename_dict[9]='T_rec_noise'
    rename_dict[10]='TAT'
    rename_dict[11]='TWTT_noise'
    rename_dict[12]='E_emi_clean'
    rename_dict[13]='X_emi_clean'
    rename_dict[14]='Y_emi_clean'
    rename_dict[15]='Z_emi_clean'
    rename_dict[16]='T_emi_clean'
    rename_dict[17]='E_rec_clean'
    rename_dict[18]='X_rec_clean'
    rename_dict[19]='Y_rec_clean'
    rename_dict[20]='Z_rec_clean'
    rename_dict[21]='T_rec_clean'
    rename_dict[22]='TWTT_clean'
    rename_dict[23]='Noise_emi'
    rename_dict[24]='Noise_rec'
    
    DF2 = DF.rename(columns=rename_dict)
    
    rename_dict = dict()
    
    if False:
        rename_dict['E_emi_noise'] = 'date_emi'
        rename_dict['X_emi_noise'] = 'E_AHD_emi'
        rename_dict['Y_emi_noise'] = 'N_AHD_emi'
        rename_dict['Z_emi_noise'] = 'D_AHD_emi'
        rename_dict['T_emi_noise'] = 'NA'
        rename_dict['E_rec_noise'] = 'date_rec'
        rename_dict['X_rec_noise'] = 'E_AHD_rec'
        rename_dict['Y_rec_noise'] = 'N_AHD_rec'
        rename_dict['Z_rec_noise'] = 'D_AHD_rec'
        rename_dict['T_rec_noise'] = 'NA'
        rename_dict['TAT'        ] = 'NA'
        rename_dict['TWTT_noise' ] = 'TWTT_obs'
    else:
        rename_dict['E_emi_clean'] = 'date_emi'
        rename_dict['X_emi_clean'] = 'E_AHD_emi'
        rename_dict['Y_emi_clean'] = 'N_AHD_emi'
        rename_dict['Z_emi_clean'] = 'D_AHD_emi'
        rename_dict['T_emi_clean'] = 'NA'
        rename_dict['E_rec_clean'] = 'date_rec'
        rename_dict['X_rec_clean'] = 'E_AHD_rec'
        rename_dict['Y_rec_clean'] = 'N_AHD_rec'
        rename_dict['Z_rec_clean'] = 'D_AHD_rec'
        rename_dict['T_rec_clean'] = 'NA'
        rename_dict['TAT'        ] = 'NA'
        rename_dict['TWTT_clean' ] = 'TWTT_obs'
    
    DF3 = DF2.rename(columns=rename_dict)
    DF3["TWTT"] = DF3["TWTT_obs"]
    
    DF3["ID_BEA"] = int(f[-7])
    
    print(DF3["ID_BEA"])
    
    DF3["VALID"]  = True
    
    DFstk.append(DF3)
    

DF4 = pd.concat(DFstk)

DF4 = DF4.reset_index()



DF4.to_csv(p + "/DF_new_format.csv")
    
######### BASELINE
    
pp = "/home/psakicki/CODES/acoustyx_toolbox_2_py3/working/SIMU_SUMMER2020_1000_R50_nois1-0.0____/SIMU_SUMMER2020_1000_R50_nois1-0.0_____.B.dat"

MM = np.loadtxt(pp)

BLstk = []
for i in range(len(MM)):
    for j in range(len(MM)):
        BLstk.append([i+1,j+1,MM[i,j]])

DFbl = pd.DataFrame(BLstk)

DFbl.rename(columns={0:"ID_BEA_A"},inplace=True)
DFbl.rename(columns={1:"ID_BEA_B"},inplace=True)
DFbl.rename(columns={2:"dist"},inplace=True)

DFbl.to_csv("/home/psakicki/CODES/acoustyx_toolbox_2_py3/working/SIMU_SUMMER2020_1000_R50_nois1-0.0____/SIMU_SUMMER2020_1000_R50_nois1-0.0_____.B.csv_RAW")