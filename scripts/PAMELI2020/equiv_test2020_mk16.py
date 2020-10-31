#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 15:18:25 2020

@author: psakicki

mk04 : retour a des Arrays, trop lent sinon
mk06 : improvement of mk04, now the BL can ben switched off and then only one Beacon
mk07 : fct are transfered in a lib
mk08 : we can work in microsec now
mk09 : add the angle
mk10 : restructure thee angle option
mk11-12 : re implementation of the dir vectors (previouly angle)
mk14 : window 
mk15 : loop
"""
#### Import star style
from geodezyx import *                   # Import the GeodeZYX modules
from geodezyx.externlib import *         # Import the external modules
from geodezyx.megalib.megalib import *   # Import the legacy modules names

from geodezyx import conv
from geodezyx import stats

import SSP as ssp
import raytrace as rt
import acouclass as acls

import numpy as np
import pandas as pd

import time

from scipy.spatial.transform import Rotation

### from pygoat.seafloorpos import ixblue_process_fcts as ibpf
import ixblue_process_fcts as ibpf

np.set_printoptions(threshold=10)

########################### LOOP PARAMETER ###################################

param_mode = 5

if param_mode == 1: 
    ###### Range initial
    Tab_estim_c                = [0,0,0,0,1,1,1,1]
    Tab_with_direction_vectors = [0,1,1,1,0,1,1,1]
    Tab_p_twtt                 = [10**-0,10**-0,10**-6,10**0]*2
    Tab_p_dir_vect             = [10**-0,10**-0,10**-0,10**-6]*2
elif param_mode == 2: #weight wide range
    ###### Influence of the DirVect weights
    Tab_p_dir_vect             = [10**float(n) for n in np.arange(-6,7)]
    Tab_p_twtt                 = np.ones(len(Tab_p_dir_vect))
    Tab_estim_c                = np.ones(len(Tab_p_dir_vect))
    Tab_with_direction_vectors = np.ones(len(Tab_p_dir_vect))
elif param_mode  == 3: #weight range better
    Tab_p_twtt_proto = np.arange(1*10**-5,10*10**-5,2*10**-5)
    ### Tab_p_dir_vect_proto = np.arange(1*10**-3,12*10**-2,2*10**-2) too loose
    Tab_p_dir_vect_proto = np.arange(1*10**-3,1.2*10**-2,2*10**-3)
    
    Tab_p_proto = np.vstack(list(itertools.product(Tab_p_twtt_proto,Tab_p_dir_vect_proto)))
    
    Tab_p_twtt                 = Tab_p_proto[:,0]
    Tab_p_dir_vect             = Tab_p_proto[:,1]

    Tab_estim_c                = np.ones(len(Tab_p_twtt))
    Tab_with_direction_vectors = np.ones(len(Tab_p_twtt))
elif param_mode  == 4: #weight range mono
    #INFO: weight TWTT       : 5e-05
    #INFO: weight Dir Vector : 0.02
    #Tab_p_twtt     = [8e-05]
    #Tab_p_dir_vect = [0.0042]    
    
    #Tab_p_twtt     = [1]
    #Tab_p_dir_vect = [1]    

    #Tab_p_twtt     = [10e-05]
    #Tab_p_dir_vect = [0.005]   

    #Tab_p_twtt     = [2.8e-05]
    #Tab_p_dir_vect = [0.0115]
    
    #Tab_p_twtt     = [3e-05]
    #Tab_p_dir_vect = [0.00015]

    #Tab_p_twtt     = [5e-05]
    #Tab_p_dir_vect = [0.0075]

    Tab_p_twtt     = [2e-05]
    Tab_p_dir_vect = [0.1]    
    Tab_p_dir_vect = [0.01]    
    Tab_p_dir_vect = [0.001]    

    Tab_estim_c                = np.ones(len(Tab_p_twtt))
    Tab_with_direction_vectors = np.ones(len(Tab_p_twtt))
    
elif param_mode == 5: #range good

    Tab_p_twtt     = [2e-05]*8
    Tab_p_dir_vect = [0,0.1,0.01,0.001]*2
    
    Tab_estim_c                = [0,0,0,0,1,1,1,1]
    Tab_with_direction_vectors = [0,1,1,1,0,1,1,1]
   
elif param_mode == 6: #range good windows
    Tab_p_twtt     = [2e-05]*6
    Tab_p_dir_vect = [0.1,0.01,0.001]*6
    Tab_estim_c                = [1,1,1,1,1,1]
    Tab_with_direction_vectors = [1,1,1,1,1,1]
    

Range_iTab_param = np.arange(0,8)  ### range initial
Range_iTab_param = [5]             ### range initial favorite 

Range_iTab_param = np.arange(0,13) ### weight wide range 
Range_iTab_param = [7]             ### weight wide range favorite 

Range_iTab_param = [0]                          ### weight range mono => worst one
Range_iTab_param = [len(Tab_p_twtt)-1]          ### weight range mono => best one
Range_iTab_param = np.arange(0,len(Tab_p_twtt)) ### weight range better

Range_iTab_path  = [9] ## 9
Range_iTab_path  = np.arange(0,3)
Range_iTab_path  = np.arange(0,12)

Range_iTab_win   = np.arange(0,5)
Range_iTab_win   = [-1] 

Range_iTab_bea = np.arange(1,4)
Range_iTab_bea = [-1]

for iTab_path,iTab_param,iTab_win,iTab_bea  in itertools.product(Range_iTab_path,
                                                                 Range_iTab_param,
                                                                 Range_iTab_win,
                                                                 Range_iTab_bea):
    
    print("################################################")
    print("iTab_param #",iTab_param)
    print("iTab_path  #",iTab_path)
    print("iTab_win   #",iTab_win)
    print("iTab_bea   #",iTab_bea)
    print("################################################")
    
    exp_suffix = "".join([str(e) for e in [iTab_path,iTab_param,iTab_win,iTab_bea]])

    ############################ LOADING DATA ####################################
    check_res = 1
    
    ############## REAL
    
    if utils.get_computer_name() == "TPX1-GFZ":
        path_computer = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/"
    else:
        path_computer = "/wrk/sakic_xchg/ggsp_pf/PLAYGROUND/psakicki/0000_PROJECTS_OTHERS/2009_PAMELi_GNSSA/"       
        path_computer = "/wrk/sakic_xchg/2009_PAMELi_GNSSA/"       
    path_bl  = path_computer + "/02_PREPROCESSING/01_/PAMELI_BREST_vJupyter_2_4transducers/PAMELI_BREST_vJupyter_2_4transducers_BASELINE.csv"
    pZ = path_computer + "/02_PREPROCESSING/01_/PAMELI_BREST_vJupyter_2_4transducers/PAMELI_BREST_vJupyter_2_4transducers.Z.dat"
    pC = path_computer + "/02_PREPROCESSING/01_/PAMELI_BREST_vJupyter_2_4transducers/PAMELI_BREST_vJupyter_2_4transducers.C.dat"    
  
    
    ##################################### TESTS 03 & 05
    path_obs = path_computer + "/02_PREPROCESSING/05_/"
    P_obs_stk = []
    ########## REAL TIME IXBLUE - Up from iXBlue
    p_obs = path_obs + "PAM_BST_v213_m2507-1158_m3_d04_bea1_RTiXBlue/PAM_BST_v213_m2507-1158_m3_d04_bea1_RTiXBlue.csv"
    P_obs_stk.append(p_obs)
    p_obs = path_obs + "PAM_BST_v223_m2507-1109_m5_d02_bea2_RTiXBlue/PAM_BST_v223_m2507-1109_m5_d02_bea2_RTiXBlue.csv"
    P_obs_stk.append(p_obs)
    p_obs = path_obs + "PAM_BST_v233_m2507-1245_m4_d06_bea3_RTiXBlue/PAM_BST_v233_m2507-1245_m4_d06_bea3_RTiXBlue.csv"
    P_obs_stk.append(p_obs)
    down_correction = 0 ### No more !!
    ########## REAL TIME IXBLUE - Up from RTKLIB directly in RGF93 (GPS only)
    p_obs = path_obs + "PAM_BST_v214_m2507-1158_m3_d04_bea1_UpPostRTiXBlue/PAM_BST_v214_m2507-1158_m3_d04_bea1_UpPostRTiXBlue.csv"
    P_obs_stk.append(p_obs)
    p_obs = path_obs + "PAM_BST_v224_m2507-1109_m5_d02_bea2_UpPostRTiXBlue/PAM_BST_v224_m2507-1109_m5_d02_bea2_UpPostRTiXBlue.csv"
    P_obs_stk.append(p_obs)
    p_obs = path_obs + "PAM_BST_v234_m2507-1245_m4_d06_bea3_UpPostRTiXBlue/PAM_BST_v234_m2507-1245_m4_d06_bea3_UpPostRTiXBlue.csv"
    P_obs_stk.append(p_obs)
    down_correction = 0 ### No more !!    
    ########## RTKLIB ITRF14 ***CORRECTED*** to RGF93
    p_obs = path_obs + "PAM_BST_v211_m2507-1158_m3_d04_bea1_ITRF14_RTKLIB/PAM_BST_v211_m2507-1158_m3_d04_bea1_ITRF14_RTKLIB.csv"
    P_obs_stk.append(p_obs)
    p_obs = path_obs + "PAM_BST_v221_m2507-1109_m5_d02_bea2_ITRF14_RTKLIB/PAM_BST_v221_m2507-1109_m5_d02_bea2_ITRF14_RTKLIB.csv"
    P_obs_stk.append(p_obs)
    p_obs = path_obs + "PAM_BST_v231_m2507-1245_m4_d06_bea3_ITRF14_RTKLIB/PAM_BST_v231_m2507-1245_m4_d06_bea3_ITRF14_RTKLIB.csv"
    P_obs_stk.append(p_obs)
    down_correction = 0
    ########## RTKLIB directly in RGF93 (multiGNSS)
    p_obs = path_obs + "PAM_BST_v211_m2507-1158_m3_d04_bea1_RGF93_RTKLIB/PAM_BST_v211_m2507-1158_m3_d04_bea1_RGF93_RTKLIB.csv"
    P_obs_stk.append(p_obs)
    p_obs = path_obs + "PAM_BST_v221_m2507-1109_m5_d02_bea2_RGF93_RTKLIB/PAM_BST_v221_m2507-1109_m5_d02_bea2_RGF93_RTKLIB.csv"
    P_obs_stk.append(p_obs)
    p_obs = path_obs + "PAM_BST_v231_m2507-1245_m4_d06_bea3_RGF93_RTKLIB/PAM_BST_v231_m2507-1245_m4_d06_bea3_RGF93_RTKLIB.csv"
    P_obs_stk.append(p_obs)
    down_correction = 0
    ########## RTKLIB directly in RGF93 (GPS only)
    p_obs = path_obs + "PAM_BST_v211_m2507-1158_m3_d04_bea1_RGF93_RTKLIB_GPSonly/PAM_BST_v211_m2507-1158_m3_d04_bea1_RGF93_RTKLIB_GPSonly.csv"
    P_obs_stk.append(p_obs)
    p_obs = path_obs + "PAM_BST_v221_m2507-1109_m5_d02_bea2_RGF93_RTKLIB_GPSonly/PAM_BST_v221_m2507-1109_m5_d02_bea2_RGF93_RTKLIB_GPSonly.csv"
    P_obs_stk.append(p_obs)
    p_obs = path_obs + "PAM_BST_v231_m2507-1245_m4_d06_bea3_RGF93_RTKLIB_GPSonly/PAM_BST_v231_m2507-1245_m4_d06_bea3_RGF93_RTKLIB_GPSonly.csv"
    P_obs_stk.append(p_obs)
    down_correction = 0

    ##################################### TESTS 04 & 06
    # path_obs = path_computer + "/02_PREPROCESSING/06_/"
    # P_obs_stk = []
    
    # # p_obs = path_obs + "PAM_BST_v401a_m21ma_d01_ITRF14_RTKLIB_Cntr/PAM_BST_v401a_m21ma_d01_ITRF14_RTKLIB_Cntr.csv"
    # # P_obs_stk.append(p_obs)
    # # p_obs = path_obs + "PAM_BST_v402a_m21so_d10_ITRF14_RTKLIB_Cntr/PAM_BST_v402a_m21so_d10_ITRF14_RTKLIB_Cntr.csv"
    # # P_obs_stk.append(p_obs)
    # # p_obs = path_obs + "PAM_BST_v403a_m22so_d11_ITRF14_RTKLIB_Cntr/PAM_BST_v403a_m22so_d11_ITRF14_RTKLIB_Cntr.csv"
    # # P_obs_stk.append(p_obs)
    # # down_correction = 0

    # p_obs = path_obs + "PAM_BST_v511_m3_d05_bea1_RGF93_RTKLIB_GPSonly_Abov/PAM_BST_v511_m3_d05_bea1_RGF93_RTKLIB_GPSonly_Abov.csv"
    # P_obs_stk.append(p_obs)
    # p_obs = path_obs + "PAM_BST_v521_m5_d03_bea2_RGF93_RTKLIB_GPSonly_Abov/PAM_BST_v521_m5_d03_bea2_RGF93_RTKLIB_GPSonly_Abov.csv"
    # P_obs_stk.append(p_obs)
    # p_obs = path_obs + "PAM_BST_v531_m4_d07_bea3_RGF93_RTKLIB_GPSonly_Abov/PAM_BST_v531_m4_d07_bea3_RGF93_RTKLIB_GPSonly_Abov.csv"
    # P_obs_stk.append(p_obs)
    # down_correction = 0

    
    ########## used p_obs ##########
    p_obs = P_obs_stk[iTab_path]
    ########## used p_obs ##########

    if "RTiXBlue" in p_obs and False:
        down_correction = -51.92176511738174
    else:
        down_correction = 0
        
    Zinp = np.loadtxt(pZ) 
    Cinp = np.loadtxt(pC)
    RS = np.random
    Cfake = Cinp + RS.rand(len(Cinp)) * 10**-5
    
    import acouclass as acls
    
    equiv_param = acls.equiv_get_param(Zinp,
                                       Cfake,
                                       zref=30,
                                       polydeg=12,
                                       zranging=(-5,5,1))
    
    # SSPpath = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/01_DATA/07-2019_IXBLUE_BREST/0_RAW_DATA/data-gaps/data/CEL_2019-07-24_08h00_utc.txt"
    # SSPpath = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/01_DATA/07-2019_IXBLUE_BREST/0_RAW_DATA/data-gaps/data/CEL_2019-07-24_11h30_utc.txt"
    # SSP = np.loadtxt(SSPpath)
    # Zinp = SSP[:,0]
    # Cinp = SSP[:,1]
    
        
    ############################ LOADING DATA ####################################
    
    
    ############################ SET ENVIRONEMENT VARIABLE #######################
    
    
    #####p_c    = 10**0
    DictIteraStore = dict()

    log_path = os.path.dirname(p_obs) + "/log"
    utils.create_dir(log_path)
    log_name = os.path.basename(p_obs)[:-4] + "_" + exp_suffix
    F = utils.Tee_frontend(log_path,log_name,print_timestamp=True)

    n_itera = 5
    p_bl   = 10**0
    
    log_name = "_".join((log_name,str(iTab_param)))
    
    estim_c                = Tab_estim_c[iTab_param]               
    with_direction_vectors = Tab_with_direction_vectors[iTab_param] 
    p_twtt                 = Tab_p_twtt[iTab_param]                
    p_dir_vect             = Tab_p_dir_vect[iTab_param]
    
    with_period_select  = False
    apply_twtt_correction_factor = True
    with_numerical_diff = False
    with_bl             = False
    fast_mode           = True ############ NOT PROPERLY IMPLEMENTED
    snelldesc_compute   = False
    snelldesc_use       = False
    with_apri_Cdict     = True
    
    if with_period_select:
        Tab_win_srt_twtt = [  0,0.2,0.4,0.6]
        Tab_win_len_twtt = [0.4,0.4,0.4,0.4]
        
        Tab_win_srt_twtt = np.arange(0,0.9,.1)
        Tab_win_len_twtt = [0.2] * len(Tab_win_srt_twtt)
        
        Tab_win_srt_twtt = np.arange(0,1.,.2)
        Tab_win_len_twtt = [0.2] * len(Tab_win_srt_twtt)
        
        
        win_srt_twtt=Tab_win_srt_twtt[iTab_win]
        win_len_twtt=Tab_win_len_twtt[iTab_win]
        win_srt_dirvect=win_srt_twtt
        win_len_dirvect=win_len_twtt
    
    ############################ SET ENVIRONEMENT VARIABLE ########################
    
    ########################### LOAD AND GENERATE DATAFRAME #######################
    
    RS = np.random.RandomState(42)
    Cfake = Cinp + RS.rand(len(Cinp)) * 10**-5
    Cfake = Cinp
        
    ### bugs in lever arms
    # p_obs = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/01_/PAMELI_BREST_vJupyter_2_4transducers/PAMELI_BREST_vJupyter_2_4transducers.csv"
    
    DFbig = pd.read_csv(p_obs,index_col=0)
    DFbig_orig = DFbig.copy()
    DFbig = DFbig.dropna()
    DFbig['VALID'] = np.ones(len(DFbig['TWTT'])).astype(bool)
    
    if with_period_select:
        DFbig['VALID'] = DFbig['VALID'].apply(np.logical_not) 
        n_win   = len(DFbig['VALID'])
        idx_srt = int(np.floor(win_srt_twtt*n_win))
        idx_len = int(np.floor(win_len_twtt*n_win))
        print("INFO:idx_str,idx_len",idx_srt,idx_len) 
        Idx = pd.Index(np.arange(idx_srt,idx_srt+idx_len))
        Idx = DFbig.iloc[Idx].index
        DFbig.loc[Idx,"VALID"] = True
    
    DFbig['date_rec'] = DFbig['date_rec'].apply(conv.string_date2dt)
    DFbig['date_emi'] = DFbig['date_emi'].apply(conv.string_date2dt)
    
    DFbl = pd.read_csv(path_bl,index_col=0)
    DFbl = DFbl.drop_duplicates('dist')
    DFbl = DFbl[DFbl.ID_BEA_A != DFbl.ID_BEA_B]
    DFbl.reset_index(inplace=True)
    
    cols_posi = ['N_AHD_emi', 'E_AHD_emi', 'D_AHD_emi',
                 'N_TDC_rec', 'E_TDC_rec', 'D_TDC_rec',
                 'N_AHD_rec', 'E_AHD_rec', 'D_AHD_rec']
    
    cols_posi_emi = cols_posi[0:3]
    cols_posi_rec = cols_posi[3:6]
    cols_posi_rec_AHD = cols_posi[3:6]
    
    ID_bea_list_orig = np.sort(DFbig['ID_BEA'].unique())
    
    if iTab_bea == -1:
        ID_bea_list = [3,4]
        ID_bea_list = [2]
        ID_bea_list = ID_bea_list_orig
    else:
        ID_bea_list = [iTab_bea]
        
    #### Creat the Dir Vector DF
    DFdirvect_big = DFbig[["date_emi","date_rec","ID_BEA",
                           "vN","vE","vD",
                           "vN_nrm","vE_nrm","vD_nrm"] + cols_posi_rec_AHD].copy()
    
    DFdirvect_big["date_rec_POSIX"] = DFdirvect_big["date_rec"].apply(conv.dt2posix)
    
    DFdirvect_big_grp = DFdirvect_big.groupby("date_emi")
    DFdirvect_big = DFdirvect_big_grp.mean()
    
    DFdirvect_big["date_rec_mean"] = DFdirvect_big_grp.mean()["date_rec_POSIX"].apply(conv.posix2dt)
    DFdirvect_big.reset_index(inplace=True)
    
    DFdirvect_big['VALID'] = np.ones(len(DFdirvect_big['ID_BEA'])).astype(bool)
    
    if with_period_select:
        n_win   = len(DFdirvect_big['VALID'])
        idx_srt = int(np.floor(win_srt_dirvect*n_win))
        idx_len = int(np.floor(win_len_dirvect*n_win))
        Idx = pd.Index(np.arange(idx_srt,idx_srt+idx_len))
        Idx = DFdirvect_big.iloc[Idx].index
    
        DFdirvect_big["VALID"]         = False    
        DFdirvect_big.loc[Idx,"VALID"] = True
        
    ########################### LOAD AND GENERATE DATAFRAME ###########################
    
    
    ############################ SET INITIAL APRI ################################
    c_ref_init = 1510.
    
    c_ref_init = ibpf.equiv_straght_ray(Zinp, Cinp, 40)
    
    if not apply_twtt_correction_factor:
        c_ref_init = c_ref_init * 10**-6
        
    Xbea_apri_all_init = np.array([[7500  , 7500  , 4010],
                                   [12500 , 12500 , 3980],
                                   [ 7500 , 12500 , 4030],
                                   [ 12500 , 7500 , 3960]])
    
    Xbea_apri_all_init = np.array([[7500  , 7500  , 4010],
                                   [12500 , 12500 , 3980],
                                   [ 12500 , 7500 , 4030],
                                   [ 7500  , 12500, 3960]])
    
    ################################ OLD COORDINATES FROM AN UNIDENTIFED RUN ##########################'
    #### END
    Xbea_apri_all_init = np.array([[24.31494786,   -7.20555817, 40.07],
                                   [-26.20205757,  -7.53913662, 40.66],
                                   [2.23923942,    22.42808402, 39.68]])
    
    #### NED (This one is theoretically the good one)
    Xbea_apri_all_init_array = np.array([[1,   -7.20555817, 24.31494786, 40.07],
                                         [2,   -7.53913662,-26.20205757, 40.66],
                                         [3,  22.42808402,  2.23923942,  39.68]])
    
    #### NED from a 1st LSQ
    Xbea_apri_all_init_array = np.array([[1,   -7.20555817, 24.31494786, 40.07],
                                         [2,   -7.53913662,-26.20205757, 40.66],
                                         [3,   22.57260714, 2.24218961, 39.011]])
    
    
    Xbea_apri_all_init = pd.DataFrame(Xbea_apri_all_init_array,
                                      columns=("ID_BEA","N","E","D"))
    
    Xbea_apri_all_init.set_index("ID_BEA",inplace=True)
    
    
    ################################ OLD COORDINATES FROM AN UNIDENTIFED RUN ##########################'
    
    ############### COORDINATES WRT [4236447.16018761 -329431.33460403 4740657.9228155 ] ###############
    ###END
    Xbea_apri_all_init_dict = {1: np.array([24.29822317,  -7.1673863 , 40.07]),
                               2: np.array([-26.21878226, -7.5009646 , 40.66]),
                               3: np.array([ 2.22251482,  22.46625595, 39.68])}
    
    ###NED
    Xbea_apri_all_init_dict = {1: np.array([-7.1673863 , 24.29822317, 40.07]),
                               2: np.array([-7.5009646 ,-26.21878226, 40.66]),
                               3: np.array([22.46625595, 2.22251482 , 39.68])}
    
    ### NED 
    ### /home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/03_g_weight_range_good_2/PAM_BST_v231_m2507-1245_m4_d06_bea3_RGF93_RTKLIB_GPSonly/log
    ### idem for bea1 n 2 
    Xbea_apri_all_init_dict = {1: np.array([  -7.09260362,   24.39713696,   37.80365583 ]),
                               2: np.array([  -7.36952229,  -26.22316284,   38.19201817]),
                               3: np.array([  22.61862272,    2.203012,     37.69790289])}

    
    Xbea_apri_all_init = pd.DataFrame(Xbea_apri_all_init_dict)
    Xbea_apri_all_init = Xbea_apri_all_init.transpose()
    Xbea_apri_all_init.columns = ["N","E","D"]
    Xbea_apri_all_init.index.rename("ID_BEA",inplace=True)
    
    C_apri_all_init_dict = {1: 1536.32220929,
                            2: 1534.79615206,
                            3: 1533.59767494}
                             
    
     
    ############### COORDINATES WRT [4236447.16018761 -329431.33460403 4740657.9228155 ] ###############
    
    lever_arm_RPY_dic = dict()
    # !!!!! CHECK LEVER ARMS and DIRECTIONs !!!!!
    # -0.039 m / +0.003 m / +1.481 m
    lever_arm_RPY_dic["GPS"]  = np.array([-0.039,+0.003,-1.481])
    lever_arm_RPY_dic["AHD"]  = np.array([0.,0.,0.3592])
    
    lever_arm_RPY_dic["TDC1"] = np.array([ 0.10735,0.,0.5095])
    lever_arm_RPY_dic["TDC2"] = np.array([-0.10735,0.,0.5095])
    lever_arm_RPY_dic["TDC3"] = np.array([0.,-0.10735,0.5733])
    lever_arm_RPY_dic["TDC4"] = np.array([0., 0.10735,0.5733])
    
    
    for i in range(2,5):
        A  = lever_arm_RPY_dic["TDC" + str(i)] - lever_arm_RPY_dic["TDC1"]
        print(np.linalg.norm(A))
    
    ############################ SET INITIAL APRI ################################
    
    
    ############################ SET APRI FOR THE 1st ITERA ######################
    
    c_ref = c_ref_init    
    
    Xbea_apri_all = Xbea_apri_all_init.copy()  
    Xbea_apri_used_init = Xbea_apri_all_init.loc[ID_bea_list].values
    
    
    ############################ SET APRI FOR THE 1st ITERA ######################
    
    for i_itera in range(n_itera):
        print("###############################################")
        print("INFO:","iteration",i_itera)
    
        if apply_twtt_correction_factor:
            DFbig['TWTT_obs'] = DFbig['TWTT'] * 10**-6
        else:
            DFbig['TWTT_obs'] = DFbig['TWTT']
            
        DFbig['TWTT_mod']          = np.zeros(len(DFbig['TWTT']))
        DFbig['TWTT_mod_equiv']    = np.zeros(len(DFbig['TWTT']))
        DFbig['TWTT_mod_sneldesc'] = np.zeros(len(DFbig['TWTT']))
    
        DFbig['B_TWTT']   = np.zeros(len(DFbig['TWTT']))
        
        DFbig['A_xbea_x'] = np.zeros(len(DFbig['TWTT']))
        DFbig['A_xbea_y'] = np.zeros(len(DFbig['TWTT']))
        DFbig['A_xbea_z'] = np.zeros(len(DFbig['TWTT']))
        DFbig['A_c']      = np.zeros(len(DFbig['TWTT']))
        
        if with_direction_vectors:
            DFdirvect_big['vN_mod'] = np.zeros(len(DFdirvect_big['vN_nrm']))
            DFdirvect_big['vE_mod'] = np.zeros(len(DFdirvect_big['vN_nrm']))
            DFdirvect_big['vD_mod'] = np.zeros(len(DFdirvect_big['vN_nrm']))
    
            DFdirvect_big['A_vN_x'] = np.zeros(len(DFdirvect_big['vN_nrm']))
            DFdirvect_big['A_vN_y'] = np.zeros(len(DFdirvect_big['vN_nrm']))
            DFdirvect_big['A_vN_z'] = np.zeros(len(DFdirvect_big['vN_nrm']))
    
            DFdirvect_big['A_vE_x'] = np.zeros(len(DFdirvect_big['vN_nrm']))
            DFdirvect_big['A_vE_y'] = np.zeros(len(DFdirvect_big['vN_nrm']))
            DFdirvect_big['A_vE_z'] = np.zeros(len(DFdirvect_big['vN_nrm']))
            
            DFdirvect_big['A_vD_x'] = np.zeros(len(DFdirvect_big['vN_nrm']))
            DFdirvect_big['A_vD_y'] = np.zeros(len(DFdirvect_big['vN_nrm']))
            DFdirvect_big['A_vD_z'] = np.zeros(len(DFdirvect_big['vN_nrm']))
    
            DFdirvect_big['B_vN'] = np.zeros(len(DFdirvect_big['vN_nrm']))
            DFdirvect_big['B_vE'] = np.zeros(len(DFdirvect_big['vN_nrm']))
            DFdirvect_big['B_vD'] = np.zeros(len(DFdirvect_big['vN_nrm']))    
        
        Xbea_apri_use     = []   
    
    
        for i_idbea, idbea in enumerate(ID_bea_list):
        
            print("INFO:","partial derivative for beacon",idbea)
            

            
            DFbea_orig = DFbig[DFbig['ID_BEA'] == idbea]
            DFbea      = DFbea_orig.copy()
            
            DFdirvect_bea_orig = DFdirvect_big[DFdirvect_big['ID_BEA'] == idbea]
            DFdirvect_bea      = DFdirvect_bea_orig.copy()
            
            if i_itera == 0:
                DFbea_good        = DFbea[DFbea['VALID']]
                _ , Bool_outliers = stats.outlier_mad(DFbea_good['TWTT_obs'])
                Bool_outliers     = pd.DataFrame(Bool_outliers)
                Bool_outliers.set_index(DFbea_good.index,inplace=True)
                
                # Set the new bools
                DFbea.loc[Bool_outliers.index,'VALID'] = Bool_outliers.values  
                
            
            print("INFO: ping valids",np.sum(DFbea['VALID']),"/",len(DFbea))
            
            xbea_apri = Xbea_apri_all.loc[idbea].values
            
            Xbea_apri_use.append(xbea_apri)
            
            c_bea  = ibpf.equiv_straght_ray(Zinp, Cinp, xbea_apri[2])
            c_used = c_ref
            
            if with_apri_Cdict and i_itera == 0:
                c_used = C_apri_all_init_dict[idbea]
             
            Jac_xbea_stk = []
            Jac_c_stk    = []
            
            iiiiii=0
            Afast_stk = []
            Afast_twtt_diff_only_stk = []
            Bfast_stk = []
            Bfast_twtt_only_stk = []
            
            n_rowbea = len(DFbea)
            
            for i_rowbea , rowbea in DFbea.iterrows():
                #print("INFO:",i_rowbea,"/",n_rowbea)
                
                if not rowbea['VALID']:
                    continue
                
                iiiiii+=1
                
                # #cols_posi = ['N_AHD_emi', 'E_AHD_emi', 'D_AHD_emi',
                # #             'N_AHD_rec', 'E_AHD_rec', 'D_AHD_rec']
    
                # cols_posi_emi = ['N_AHD_emi', 'E_AHD_emi', 'D_AHD_emi']
                # tdcstr = "TDC" + str(rowbea['ID_TDC'])
                # cols_posi_rec = ['N_' +tdcstr+ '_rec',
                #                  'E_' +tdcstr+ '_rec',
                #                  'D_' +tdcstr+ '_rec']
                
                xemi = np.array(rowbea[cols_posi_emi].values,dtype=float)
                xrec = np.array(rowbea[cols_posi_rec].values,dtype=float)
                
                
                #print("equiv",acls.equiv_inverse(equiv_param,xemi,xbea_apri))

                            
                #### !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                xemi[2] += down_correction
                xrec[2] += down_correction
                #### !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                
                twtt_obs = rowbea['TWTT_obs']
                                  
                ###################### START MODELING & DERIVATION ###############
                
                args_for_partial_dev =  (xbea_apri[0],
                                         xbea_apri[1],
                                         xbea_apri[2],
                                         xemi,
                                         xrec,
                                         c_used)
                
                ########## Modeled observations
                twtt_mod_equiv = ibpf.fct_obs_equiv(*args_for_partial_dev)
                
                # twtt_mod_sneldesc_leg = rt.raytrace_seek((xemi,xrec),
                #                                          xbea_apri,
                #                                          Zinp,
                #                                          Cfake,
                #                                          legacy=True,
                #                                          verbose=False,
                #                                          fulloutput=False)
                
                if snelldesc_compute:
                    twtt_mod_sneldesc_new = rt.raytrace_seek((xemi,xrec),
                                                             xbea_apri,
                                                             Zinp,
                                                             Cfake,
                                                             legacy=False,
                                                             verbose=False,
                                                             fulloutput=False)
                else:
                    twtt_mod_sneldesc_new = (0,0)
                    
                if snelldesc_use:
                    twtt_mod = twtt_mod_sneldesc_new[-1]
                else:
                    twtt_mod = twtt_mod_equiv
                
                if not fast_mode:
                    rowbea['TWTT_mod_equiv']    = twtt_mod_equiv
                    rowbea['TWTT_mod_sneldesc'] = twtt_mod_sneldesc_new[-1]
                    
                    if snelldesc_use:
                        rowbea['TWTT_mod']          = rowbea['TWTT_mod_sneldesc'] 
                    else:
                        rowbea['TWTT_mod']          = rowbea['TWTT_mod_equiv'] 
        
                    rowbea['B_TWTT'] = rowbea['TWTT_obs'] - rowbea['TWTT_mod']  
                    
                    #if with_direction_vectors:
                    #    rowbea[['vN_mod','vE_mod','vD_mod']] = dir_vect_mod
                    #    rowbea[['B_vN','B_vE','B_vD']]       = B_dir_vect
                        
                else:
                    Bfast_stk.append(twtt_obs - twtt_mod)
                    Bfast_twtt_only_stk.append(twtt_obs - twtt_mod)
    
                ########## derivative of position
                if with_numerical_diff:
                    ######### Numerical derivation
                    jac_xbea_line = []
                    for i in range(3):
                        jac_xbea = stats.partial_derive(fct_obs_equiv, i,
                                                        args_f = args_for_partial_dev)
                        jac_xbea_line.append(jac_xbea)
    
                    ########## derivative of c
                    jac_c = stats.partial_derive(fct_obs_equiv, 'c_bea_in',
                                                 args_f = args_for_partial_dev)  
    
                    if not fast_mode:
                        rowbea[['A_xbea_x','A_xbea_y','A_xbea_z','A_c']] = jac_xbea_line[0],jac_xbea_line[1],jac_xbea_line[2],jac_c
    
                    else:
                        Afast_row = [jac_xbea_line[0],
                                     jac_xbea_line[1],
                                     jac_xbea_line[2],
                                     jac_c]
                        
                        Afast_stk.append(Afast_row)
    
                else:
                    ######### Analytic derivation
                    jac_ana_line = np.array(ibpf.fct_obs_equiv_deriv_ana(*args_for_partial_dev))
                    
                    if not fast_mode:
                        rowbea[['A_xbea_x','A_xbea_y','A_xbea_z','A_c']] = jac_ana_line[0],jac_ana_line[1],jac_ana_line[2],jac_ana_line[3]
                    else:
                        Afast_row = [jac_ana_line[0],
                                     jac_ana_line[1],
                                     jac_ana_line[2],
                                     jac_ana_line[3]]
                        
                        Afast_stk.append(Afast_row)
                        Afast_twtt_diff_only_stk.append(Afast_row)
                        
                ######### Symbolic derivation
                #jac_sym = fct_obs_equiv_deriv_symbo(*args_for_partial_dev)
                            
                ####### Update the row
                if not fast_mode:
                    DFbea.loc[i_rowbea] = rowbea
                    
            if fast_mode:
                Afast = np.array(Afast_stk)
                Bfast = np.array(Bfast_stk)
    
                Bool4Afast = (DFbea["ID_BEA"] == idbea) & (DFbea["VALID"] == True)
                DFbea.loc[Bool4Afast,['A_xbea_x','A_xbea_y','A_xbea_z','A_c']] = Afast_twtt_diff_only_stk
                DFbea.loc[Bool4Afast,"B_TWTT"] = Bfast_twtt_only_stk
    
            
            ###################### END MODELING & DERIVATION #######################               
            
            #### Outlier cleaning
            # Detect the outliers
            DFbea_good        = DFbea[DFbea['VALID']]
            _ , Bool_outliers = stats.outlier_mad(DFbea_good['B_TWTT'])
            Bool_outliers     = pd.DataFrame(Bool_outliers)
            Bool_outliers.set_index(DFbea_good.index,inplace=True)
            
            # Set the new bools
            DFbea.loc[Bool_outliers.index,'VALID'] = Bool_outliers.values  
            
            ####### Update the big DataFrame
            DFbig.loc[DFbea.index] = DFbea
            
        
        ######################### BASELINES #########################
    
        DFbl['dist_mod'] = np.zeros(len(DFbl))
        DFbl['B_bl']     = np.zeros(len(DFbl))
        
        if with_bl:
            DFbl['A_bl_A_x'] = np.zeros(len(DFbl))
            DFbl['A_bl_A_y'] = np.zeros(len(DFbl))
            DFbl['A_bl_A_z'] = np.zeros(len(DFbl))
            
            DFbl['A_bl_B_x'] = np.zeros(len(DFbl))
            DFbl['A_bl_B_y'] = np.zeros(len(DFbl))
            DFbl['A_bl_B_z'] = np.zeros(len(DFbl))
        
            for i_rowbl , rowbl in DFbl.iterrows():
            
                xbea_A = Xbea_apri_all.loc[rowbl.ID_BEA_A].values
                xbea_B = Xbea_apri_all.loc[rowbl.ID_BEA_B].values
                
                dist_mod = conv.dist(xbea_A, xbea_B)
                dist_obs = rowbl['dist'] 
                
                xbea_A_diff , xbea_B_diff = conv.dist_diff(xbea_A, xbea_B)
                
                rowbl['dist_mod'] = dist_mod
                rowbl['B_bl']     = dist_obs - dist_mod
                
                rowbl['A_bl_A_x'] = xbea_A_diff[0]
                rowbl['A_bl_A_y'] = xbea_A_diff[1]
                rowbl['A_bl_A_z'] = xbea_A_diff[2]
                
                rowbl['A_bl_B_x'] = xbea_B_diff[0]
                rowbl['A_bl_B_y'] = xbea_B_diff[1]
                rowbl['A_bl_B_z'] = xbea_B_diff[2]
                        
                ####### Update the row
                DFbl.loc[i_rowbl] = rowbl
                
    ######################### DIRECTION VECTORS #########################
    
        Bdir_vect_stk = []    
        Adir_vect_stk = []
        ### O for Obseervations
        Odir_vect_stk = []
        
        Irow_stk = []
        
        print("INFO: dir vect valids",np.sum(DFdirvect_bea['VALID']),"/",len(DFdirvect_bea))
        
        if with_direction_vectors: 
            Valid_indexes_dirvect = []
            prev_rowbea_date_emi = None
            for i_rowdirvectbea , rowdirvectbea in DFdirvect_bea.iterrows():
                
                if rowdirvectbea.date_emi != prev_rowbea_date_emi:
                    prev_rowbea_date_emi = rowdirvectbea.date_emi
                    add_dir_vect_to_AnB = True
                else:
                    add_dir_vect_to_AnB = False
                    
                
                if not rowdirvectbea['VALID']:
                    continue
                else:
                    Valid_indexes_dirvect.append(i_rowdirvectbea)
                    
                
                xrec_AHD = rowdirvectbea[cols_posi_rec_AHD].values.astype(np.float64)
                xrec_AHD[2] += down_correction
                dir_vect_obs = rowdirvectbea[["vN_nrm","vE_nrm","vD_nrm"]].values.astype(np.float64)
            
                args_for_partial_dev_dir_vect = (xbea_apri[0],
                                                 xbea_apri[1],
                                                 xbea_apri[2],
                                                 xrec_AHD)      
                
                dir_vect_mod = ibpf.fct_obs_dir_vect(*args_for_partial_dev_dir_vect)
                dir_vect_mod = np.array(dir_vect_mod)
                
                
                dir_vect_siteazi_mod = np.rad2deg(ibpf.ned2site_azim(dir_vect_mod))
                dir_vect_siteazi_obs = np.rad2deg(ibpf.ned2site_azim(dir_vect_obs))

                #print(dir_vect_siteazi_mod , dir_vect_siteazi_obs)                
                #print(dir_vect_siteazi_mod , dir_vect_siteazi_obs)
                        
                b_dir_vect = dir_vect_obs - dir_vect_mod
                #from pygoat.seafloorpos import ixblue_process_fcts as ibpf
                    
                DirVectDIFF1 = ibpf.fct_obs_dir_vect_deriv_ana(*args_for_partial_dev_dir_vect)
                if False:
                    DirVectDIFF2 = ibpf.fct_obs_dir_vect_deriv_symbo(*args_for_partial_dev_dir_vect) 
                    boll_equiv   = np.all(np.isclose(np.array(DirVectDIFF1,np.float64) , np.array(DirVectDIFF2)))
                    print(boll_equiv)
                
                DirVectDIFF  = DirVectDIFF1
                
                if add_dir_vect_to_AnB:
                    
                    Bdir_vect_stk = Bdir_vect_stk + list(b_dir_vect)
                    Odir_vect_stk = Odir_vect_stk + list(dir_vect_obs)
    
                    Irow_stk.append(i_rowdirvectbea)
                    
                    if not fast_mode:
                        rowdirvectbea[['A_vN_x','A_vN_y','A_vN_z']] = DirVectDIFF[0:3]
                        rowdirvectbea[['A_vE_x','A_vE_y','A_vE_z']] = DirVectDIFF[3:6]
                        rowdirvectbea[['A_vD_x','A_vD_y','A_vD_z']] = DirVectDIFF[6:9]
                        
                    else:
                        if estim_c:
                            lc = [0]
                        else:
                            lc = []
                        
                        Afast_row = list(DirVectDIFF[0:3]) + lc
                        Adir_vect_stk.append(Afast_row)
                        Afast_row = list(DirVectDIFF[3:6]) + lc
                        Adir_vect_stk.append(Afast_row)
                        Afast_row = list(DirVectDIFF[6:9]) + lc
                        Adir_vect_stk.append(Afast_row)
        
        
            Bdir_vect = np.array(Bdir_vect_stk)
            Adir_vect = np.array(Adir_vect_stk)
            Odir_vect = np.array(Odir_vect_stk)
        
        
        
        ################# BUILD FINAL DESIGN MATRIX AND VECTOR ###################
        A_c_stk    = []
        A_xbea_stk = []
        B_twtt_stk = []
        
        for i_idbea, idbea in enumerate(ID_bea_list):
            
            Bool_idbea = ((DFbig['ID_BEA'] == idbea) & DFbig['VALID']) 
            
            B_twtt_stk.append(DFbig.loc[Bool_idbea,'B_TWTT'].values)
            A_c_stk.append(DFbig.loc[Bool_idbea,'A_c'].values)
            A_xbea_stk.append(DFbig.loc[Bool_idbea,['A_xbea_x',
                                                    'A_xbea_y',
                                                    'A_xbea_z']].values)  
            
        if with_bl:
            A_bl_stk = []
            B_bl_stk = []
            for inat_rowbl , (i_rowbl , rowbl) in enumerate(DFbl.iterrows()):
                
                if estim_c:
                    sc = 1
                else:
                    sc = 0
                
                A_bl_line = np.zeros(3*len(ID_bea_list) + sc)
                iA = int(rowbl.ID_BEA_A-1)
                iB = int(rowbl.ID_BEA_B-1)
                
                A_bl_line[3*iA:3*iA+3] = DFbl.loc[i_rowbl,['A_bl_A_x',
                                                           'A_bl_A_y',
                                                           'A_bl_A_z']]
                
                A_bl_line[3*iB:3*iB+3] = DFbl.loc[i_rowbl,['A_bl_B_x',
                                                           'A_bl_B_y',
                                                           'A_bl_B_z']]
                A_bl_stk.append(A_bl_line)
                B_bl_stk.append(rowbl['B_bl'])
            
    
        #### FINAL STACK
        A1 = scipy.linalg.block_diag(*A_xbea_stk)
        A2 = np.hstack(A_c_stk)
        
        if estim_c and with_bl:
            A3 = np.column_stack((A1,A2))
            A  = np.vstack((A3,*A_bl_stk))
        elif estim_c and not with_bl:
            A  = np.column_stack((A1,A2))
        elif not estim_c and with_bl:
            A = np.vstack((A1,np.vstack(A_bl_stk)))
        elif not estim_c and not with_bl:
            A = A1
        
        B_twtt = np.hstack(B_twtt_stk)
        n_B_twtt = len(B_twtt)
        
        if with_bl:
            B = np.hstack((B_twtt,B_bl_stk))
        else:
            B = B_twtt
            
        if with_direction_vectors:
            Atwtt_bkp = A.copy()
            Btwtt_bkp = B.copy()
            
            A = np.vstack((A,Adir_vect))
            B = np.hstack((B,Bdir_vect))
            n_B_dir_vect = len(Bdir_vect)
            
        if i_itera == 0:
            fuvi=1
        
        if with_bl:
            K , Q , P = stats.weight_mat([p_twtt,p_bl],
                                         [len(B_twtt),len(B_bl_stk)],
                                         fuvinp=fuvi)
        elif with_direction_vectors:
            K , Q , P = stats.weight_mat([p_twtt,p_dir_vect],
                                         [len(B_twtt),len(Bdir_vect)],
                                         fuvinp=fuvi)
        else:
            K , Q , P = stats.weight_mat([p_twtt],
                                         [len(B_twtt)],
                                         fuvinp=fuvi)  
            
        #P = P/fuv
            

        ################# BUILD FINAL DESIGN MATRIX AND VECTOR ###################
    
    
        ############################# INVERSION ##################################    
        At   = A.T
        
        N    = At.dot(P).dot(A)
        Ninv = scipy.linalg.inv(N)
        
        dX   = Ninv.dot(At).dot(P).dot(B)
        ############################# INVERSION ##################################    
    
    
        ############################# GET NEW VALUES #############################        
        X = np.hstack(Xbea_apri_use)
        Xprev = X
        
        if estim_c:
            X = np.append(X,c_used)
        
        Xnew = X + dX
        
        if estim_c:
            # Xnew2 is Xnew without the c element
            Xnew2     = Xnew[:-1]
            c_ref_new = Xnew[-1]
        else:
            Xnew2     = Xnew
            c_ref_new = c_ref
    
        Xbea_apri_all_new = Xnew2.reshape(len(ID_bea_list),3)
            
        V = B - A.dot(dX)
        
        dof = (A.shape[0] - A.shape[1])
        
        fuv = V.dot(P).dot(V) / (dof)
        Varcovar = fuv * Ninv
        Sigma    = np.sqrt(np.diag(Varcovar))
            
        
        Xnew_reshape = Xnew2.reshape(len(ID_bea_list),3)
        BLnew = reffram.BL_from_points(Xnew_reshape)
        
        
        CHI2 = stats.chi2_test_lsq(V,A,P)
        
        B_curiter       = DFbig["B_TWTT"]
        B_curiter_bool  = DFbig["VALID"]
        B_curiter_valid = B_curiter[B_curiter_bool] 
        
                
        ############################# GET NEW VALUES #############################        
        
        
        ############################ STORE NEW VALUES ############################ 
        DictIteraStore_iter = dict()
        DictIteraStore[i_itera] = DictIteraStore_iter
    
        DictIteraStore_iter['globals']           = utils.globals_filtered().copy()
        DictIteraStore_iter['beacon']            = idbea
        DictIteraStore_iter['Xbea_apri_all']     = Xbea_apri_all
        DictIteraStore_iter['c_ref']             = c_ref
        DictIteraStore_iter['Xbea_apri_all_new'] = Xbea_apri_all_new
        DictIteraStore_iter['c_ref_new']         = c_ref_new
        DictIteraStore_iter['V']                 = V
        DictIteraStore_iter['B']                 = B
        DictIteraStore_iter['dX']                = dX
        DictIteraStore_iter['Xapri']             = Xbea_apri_used_init
        DictIteraStore_iter['Xprev']             = Xprev
        DictIteraStore_iter['Xnew']              = Xnew
        DictIteraStore_iter['Xnew_reshape']      = Xnew_reshape
        DictIteraStore_iter['BLnew']             = BLnew
        DictIteraStore_iter['BLobs']             = DFbl["dist"]
        DictIteraStore_iter['Sigma']             = Sigma
        DictIteraStore_iter['fuv']               = fuv
        DictIteraStore_iter['chi2']              = CHI2
        DictIteraStore_iter['p_twtt']            = p_twtt
        DictIteraStore_iter['p_dir_vect']        = p_dir_vect        
        DictIteraStore_iter['DFbig']             = copy.deepcopy(DFbig)
        

    
        ############################ STORE NEW VALUES ############################     
        print("INFO: estim_c           :",estim_c)
        print("INFO: direction_vectors :",with_direction_vectors)
        print("INFO: weight TWTT       :",p_twtt)
        print("INFO: weight Dir Vector :",p_dir_vect)
        if with_period_select:
            print("INFO: period selected TWTT       :",win_srt_twtt,win_len_twtt)
            print("INFO: period selected Dir Vector :",win_srt_dirvect,win_len_dirvect)
            DictIteraStore_iter['win_srt_twtt'] = win_srt_twtt
            DictIteraStore_iter['win_len_twtt'] = win_len_twtt
            DictIteraStore_iter['win_srt_dirvect'] = win_srt_dirvect
            DictIteraStore_iter['win_len_dirvect'] = win_len_dirvect
    
        ####################### COMPUTE NATURAL RESIDUALS ########################
        B_TWTT_Stk = []
        for i_idbea, idbea in enumerate(ID_bea_list):
        
            print("INFO:","Natural residuals for beacon",idbea)
            
            DFbea = DFbig[(DFbig['ID_BEA'] == idbea) & (DFbig['VALID'])]
            
            xbea_apri = Xbea_apri_all_new[i_idbea-1]
            
            Xemi     = DFbea[cols_posi_emi].values
            Xrec     = DFbea[cols_posi_rec].values
            if apply_twtt_correction_factor:
                TWTT_obs = DFbea['TWTT'].values * 10**-6
            else:
                TWTT_obs = DFbea['TWTT'].values            
             
            TWTT_mod_stk = []
            RESstk = []   
            
            c_bea = ibpf.equiv_straght_ray(Zinp, Cinp, xbea_apri[2])
            
            c_used = c_ref_new
             
            Jac_xbea_stk = []
            Jac_c_stk    = []
        
            ZIP = list(zip(Xemi,Xrec,TWTT_obs))    
            for xemi,xrec,twtt_obs in ZIP:
            
                ###################### DERIVATION #######################
                args_for_partial_dev =  (xbea_apri[0],
                                         xbea_apri[1],
                                         xbea_apri[2],
                                         xemi,
                                         xrec,
                                         c_used)
                
                ########## Modeled observations
                twtt_mod = ibpf.fct_obs_equiv(*args_for_partial_dev)
                TWTT_mod_stk.append(twtt_mod)
            
            ### Stacking the values in Matrices/Vectors
            TWTT_mod     = np.array(TWTT_mod_stk)
            B_TWTT       = (TWTT_obs - TWTT_mod)
            B_TWTT_Stk.append(B_TWTT)
            
        B_TWTT_Stk = np.hstack(B_TWTT_Stk)
        DictIteraStore_iter['Vnatural'] = B_TWTT_Stk
        Vnatural = DictIteraStore_iter['Vnatural']
        
        
        ####################### COMPUTE NATURAL RESIDUALS ########################
    
        
        ##################### REPLACE APRIORI WITH NEW VALUES ####################  
              
        Xbea_apri_all = Xbea_apri_all_new
        
        Xbea_apri_all = pd.DataFrame(np.column_stack((ID_bea_list,Xbea_apri_all_new)),
                                      columns=("ID_BEA","N","E","D"))
        Xbea_apri_all.set_index("ID_BEA",inplace=True)
        
        
        c_ref = c_ref_new
        ##################### REPLACE APRIORI WITH NEW VALUES ####################   
    
    
        ######### COMPUTE ANGLE FRIENDLY DIRECVTION VECTOR AND OUTLIERS #########
            
        if with_direction_vectors:
            O1 = Odir_vect.reshape((int(len(Odir_vect)/3),3))
            B1 = B[n_B_twtt:].reshape((int(len(Odir_vect)/3),3))
            V1 = V[n_B_twtt:].reshape((int(len(Odir_vect)/3),3))
            
            O2 = O1 + V1
            
            Angle_residuals = []
            
            for o1,o2 in zip(O1,O2):
                scalprod = np.dot(o1,o2)
                scalprod_norm = scalprod / (np.linalg.norm(o1) * np.linalg.norm(o2))
                angle_residual = np.rad2deg(np.arccos(scalprod_norm))
                Angle_residuals.append(angle_residual)
                
                Angle_good , Angle_bool = stats.outlier_mad(Angle_residuals)
                
            DFdirvect_big.loc[Valid_indexes_dirvect,"VALID"] = Angle_bool

            DictIteraStore_iter['DFdirvect_big'] = DFdirvect_big
            DictIteraStore_iter['Rdv_angle'] = Angle_residuals
            DictIteraStore_iter['Odv1'] = O1
            DictIteraStore_iter['Odv2'] = O2
            DictIteraStore_iter['Vdv']  = V1
            DictIteraStore_iter['Bdv']  = B1

        ###### Last Saves
        DictIteraStore_iter['globals_end'] = utils.globals_filtered().copy()
    
            
        ######### COMPUTE ANGLE FRIENDLY DIRECVTION VECTOR AND OUTLIERS #########
        print("Xapri    ",i_itera,DictIteraStore[i_itera]["Xapri"])
        print("Xprev    ",i_itera,DictIteraStore[i_itera]["Xprev"])
        print("Xnew     ",i_itera,DictIteraStore[i_itera]["Xnew"])
        print("dX       ",i_itera,DictIteraStore[i_itera]["dX"])
        print("Sigma    ",i_itera,DictIteraStore[i_itera]["Sigma"])
        print("B        ",i_itera,DictIteraStore[i_itera]["B"])
        print("V        ",i_itera,DictIteraStore[i_itera]["V"])
        print("Vnat     ",i_itera,DictIteraStore[i_itera]["Vnatural"])
        print("std BTWTT",i_itera,np.std(B_curiter))
        print("S B**2   ",i_itera,np.sum(DictIteraStore[i_itera]["B"]**2))
        print("S V**2   ",i_itera,np.sum(DictIteraStore[i_itera]["V"]**2))
        print("RMS B    ",i_itera,stats.rms_mean(DictIteraStore[i_itera]["B"]**2))
        print("RMS V    ",i_itera,stats.rms_mean(DictIteraStore[i_itera]["V"]**2))
        print("S Vnat**2",i_itera,np.sum(DictIteraStore[i_itera]["Vnatural"]**2))
        print("BL new   ",i_itera,acls.BL_from_PXPlist(Xbea_apri_all_new))
        print("BL obs   ",i_itera,DictIteraStore_iter['BLobs'])
        print("Barycntr ",i_itera,acls.barycenter_calc(Xbea_apri_all_new))
        print("c_ref_old",i_itera,DictIteraStore_iter['c_ref'] )
        print("c_ref_new",i_itera,DictIteraStore_iter['c_ref_new'] )
        print("fuv      ",i_itera,DictIteraStore_iter['fuv'] )
        print("chi2     ",i_itera,DictIteraStore_iter['chi2'] )
        
        B_lastiter = DictIteraStore[i_itera]['DFbig']["B_TWTT"]

    time.sleep(1)
    utils.pickle_saver(DictIteraStore, 
                       outdir = log_path, 
                       outname = log_name,
                       ext='.pik',
                       timestamp = True)
    
    F.stop()
        
    print("***********************************************************************")
        
    #B_1stiter        = DictIteraStore[0]['DFbig']["B_TWTT"]
    #B_1stiter_bool   = DictIteraStore[0]['DFbig']["VALID"]
    #B_1stiter_valid  = B_1stiter[B_1stiter_bool] 
    #B_lastiter       = DictIteraStore[i_itera]['DFbig']["B_TWTT"]
    #B_lastiter_bool  = DictIteraStore[i_itera]['DFbig']["VALID"]
    #B_lastiter_valid = B_lastiter[B_lastiter_bool] 
    
    #nbin = 50
    
    #fighist,axhist = plt.subplots()
    #fighist.suptitle("Histogram of the TWTT residuals")
    #axhist.hist(B_1stiter_valid,nbin,density=False,label="1st iter",color='C3')
    #xpdf,ypdf = utils.gaussian_for_plot(B_1stiter_valid,False,nbin)
    #axhist.plot(xpdf, ypdf,c='C1')
    
    #axhist.hist(B_lastiter,nbin,density=False,label="last iter",color='C0')
    #xpdf,ypdf = utils.gaussian_for_plot(B_lastiter_valid,False,nbin)
    #axhist.plot(xpdf,ypdf,c='C1')
        
    
    # for i in range(n_itera):
    #     print("B   ",i,DictIteraStore[i]["B"])
    #     print("V   ",i,DictIteraStore[i]["V"])
    #     print("Vnat",i,DictIteraStore[i]["Vnatural"])
    #     print("S B**2   ",i,np.sum(DictIteraStore[i]["B"]**2))
    #     print("S V**2   ",i,np.sum(DictIteraStore[i]["V"]**2))
    #     print("S Vnat**2",i,np.sum(DictIteraStore[i]["Vnatural"]**2))
    

