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

from scipy.spatial.transform import Rotation

from pygoat.seafloorpos import ixblue_process_fcts as ibpf


############################ LOADING DATA ####################################
check_res = 1

############## SIMU
pC = "/home/psakicki/CODES/acoustyx_toolbox_2_py3/working/SIMU_SUMMER2020_1000_R50_nois1-0.0____/SIMU_SUMMER2020_1000_R50_nois1-0.0_____.C.dat"
pZ = "/home/psakicki/CODES/acoustyx_toolbox_2_py3/working/SIMU_SUMMER2020_1000_R50_nois1-0.0____/SIMU_SUMMER2020_1000_R50_nois1-0.0_____.Z.dat"
p_obs = "/home/psakicki/CODES/acoustyx_toolbox_2_py3/working/SIMU_SUMMER2020_1000_R50_nois1-0.0____/SIMU_SUMMER2020_1000_R50_nois1-0.0_____.PXP1.P.csv"
p_obs = "/home/psakicki/CODES/acoustyx_toolbox_2_py3/working/SIMU_SUMMER2020_3x333_x500_y500_ang0_nois1-0.0____/DF_new_format.csv"
p_bl  = "/home/psakicki/CODES/acoustyx_toolbox_2_py3/working/SIMU_SUMMER2020_1000_R50_nois1-0.0____/SIMU_SUMMER2020_1000_R50_nois1-0.0_____.B.csv"

############## REAL
pZ = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/01_/PAMELI_BREST_vJupyter_2_4transducers/PAMELI_BREST_vJupyter_2_4transducers.Z.dat"
pC = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/01_/PAMELI_BREST_vJupyter_2_4transducers/PAMELI_BREST_vJupyter_2_4transducers.C.dat"
### Corrections in lever arms
p_obs = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/01_/PAMELI_BREST_vJupyter_5_4transducers_new_lever_arms/PAMELI_BREST_vJupyter_5_4transducers_new_lever_arms.csv"
p_bl  = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/01_/PAMELI_BREST_vJupyter_2_4transducers/PAMELI_BREST_vJupyter_2_4transducers_BASELINE.csv"
p_obs = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/01_/PAMELI_BREST_vJupyter_5_4transducers_new_lever_arms_GINS_PPP_coords/PAMELI_BREST_vJupyter_5_4transducers_new_lever_arms.csv"
p_obs = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/01_/PAMELI_BREST_vJupyter_7_4transducers_newRF_newDTime_BAD_BEFORE_TIMING_CORR/PAMELI_BREST_vJupyter_7_4transducers_newRF_newDTime_BAD_BEFORE_TIMING_CORR.csv"

p_obs = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/01_/PAMELI_BREST_vJupyter_6_4transducers_newRF_newDTime/PAMELI_BREST_vJupyter_6_4transducers_newRF_newDTime.csv"
p_obs = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/01_/PAMELI_BREST_vJupyter_8_4transducers_AmbigCorrect01/PAMELI_BREST_vJupyter_8_4transducers_AmbigCorrect01.csv"

p_att = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/01_/PAMELI_BREST_vJupyter_5_4transducers_new_lever_arms/PAMELI_BREST_vJupyter_5_4transducers_new_lever_arms_AttitudeInterpolator.pik"

######## Version Ambig Corriged
p_obs = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/02_/PAM_BREST_AmbigCorr02/PAM_BREST_AmbigCorr02_bea1_mk01.csv"
p_obs = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/02_/PAM_BREST_AmbigCorr02/PAM_BREST_AmbigCorr02_bea2_mk01.csv"

############### Version 02_ sept 2020
## p_obs = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/02_/PAM_BST_v221_m2507-1109_n02_bea2_ITRF14_RTKLIB/PAM_BST_v221_m2507-1109_n02_bea2_ITRF14_RTKLIB.csv"

path_obs = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/02_/"

p_obs = path_obs + "PAM_BST_v201_m2507-1446_m2_d11_4bea_ITRF14_RTKLIB/PAM_BST_v201_m2507-1446_m2_d11_4bea_ITRF14_RTKLIB.csv"

p_obs = path_obs + "PAM_BST_v232_m2507-1158_m4_d06_bea3_RTiXBlue/PAM_BST_v232_m2507-1158_m4_d06_bea3_RTiXBlue.csv"

########## RTKLIB ITRF14 ***CORRECTED*** to RGF93
p_obs = path_obs + "PAM_BST_v231_m2507-1245_m4_d06_bea3_ITRF14_RTKLIB/PAM_BST_v231_m2507-1245_m4_d06_bea3_ITRF14_RTKLIB.csv"
p_obs = path_obs + "PAM_BST_v221_m2507-1109_m5_d02_bea2_ITRF14_RTKLIB/PAM_BST_v221_m2507-1109_m5_d02_bea2_ITRF14_RTKLIB.csv"
p_obs = path_obs + "PAM_BST_v211_m2507-1158_m3_d04_bea1_ITRF14_RTKLIB/PAM_BST_v211_m2507-1158_m3_d04_bea1_ITRF14_RTKLIB.csv"
########## REAL TIME IXBLUE
p_obs = path_obs + "PAM_BST_v232_m2507-1245_m4_d06_bea3_RTiXBlue/PAM_BST_v231_m2507-1245_m4_d06_bea3_RTiXBlue.csv"
p_obs = path_obs + "PAM_BST_v221_m2507-1109_m5_d02_bea2_RTiXBlue/PAM_BST_v221_m2507-1109_m5_d02_bea2_RTiXBlue.csv"
p_obs = path_obs + "PAM_BST_v211_m2507-1158_m3_d04_bea1_RTiXBlue/PAM_BST_v211_m2507-1158_m3_d04_bea1_RTiXBlue.csv"


Zinp = np.loadtxt(pZ) 
Cinp = np.loadtxt(pC)


log_path = os.path.dirname(p_obs) + "/log"
utils.create_dir(log_path)
log_name = os.path.basename(p_obs)[:-4]


# SSPpath = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/01_DATA/07-2019_IXBLUE_BREST/0_RAW_DATA/data-gaps/data/CEL_2019-07-24_08h00_utc.txt"
# SSPpath = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/01_DATA/07-2019_IXBLUE_BREST/0_RAW_DATA/data-gaps/data/CEL_2019-07-24_11h30_utc.txt"
# SSP = np.loadtxt(SSPpath)
# Zinp = SSP[:,0]
# Cinp = SSP[:,1]

RS = np.random.RandomState(42)
Cfake = Cinp + RS.rand(len(Cinp)) * 10**-5
Cfake = Cinp

IAtt = utils.pickle_loader(p_att)

### bugs in lever arms
# p_obs = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/01_/PAMELI_BREST_vJupyter_2_4transducers/PAMELI_BREST_vJupyter_2_4transducers.csv"

DFbig = pd.read_csv(p_obs,index_col=0)
DFbig_orig = DFbig.copy()
DFbig = DFbig.dropna()
DFbig['VALID'] = np.ones(len(DFbig['TWTT'])).astype(bool)

DFbig['date_rec'] = DFbig['date_rec'].apply(conv.string_date2dt)
DFbig['date_emi'] = DFbig['date_emi'].apply(conv.string_date2dt)

DFbl = pd.read_csv(p_bl,index_col=0)
DFbl = DFbl.drop_duplicates('dist')
DFbl = DFbl[DFbl.ID_BEA_A != DFbl.ID_BEA_B]
DFbl.reset_index(inplace=True)

############################ LOADING DATA ####################################

F = utils.Tee_frontend(log_path,log_name,print_timestamp=True)

############################ SET ENVIRONEMENT VARIABLE #######################
cols_posi = ['N_AHD_emi', 'E_AHD_emi', 'D_AHD_emi',
             'N_TDC_rec', 'E_TDC_rec', 'D_TDC_rec',
             'N_AHD_rec', 'E_AHD_rec', 'D_AHD_rec']

cols_posi_emi = cols_posi[0:3]
cols_posi_rec = cols_posi[3:6]
cols_posi_rec_AHD = cols_posi[3:6]


ID_bea_list_orig = np.sort(DFbig['ID_BEA'].unique())

ID_bea_list = [3,4]
ID_bea_list = [2]
ID_bea_list = ID_bea_list_orig

n_itera = 5

DictIteraStore = dict()

p_twtt     = 10**-0
p_dir_vect = 10**-0

#####p_c    = 10**0
p_bl   = 10**0


with_numerical_diff = False
with_bl             = False


### if only one beacon
### estim_c = False
### with_bl = False

fast_mode         = True ############ NOT PROPERLY IMPLEMENTED
snelldesc_compute = True
snelldesc_use     = False
estim_c           = 1
apply_twtt_correction_factor = True
with_direction_vectors = 0
############################ SET ENVIRONEMENT VARIABLE #######################

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


Xbea_apri_all_init = pd.DataFrame(Xbea_apri_all_init_dict)
Xbea_apri_all_init = Xbea_apri_all_init.transpose()
Xbea_apri_all_init.columns = ["N","E","D"]
Xbea_apri_all_init.index.rename("ID_BEA",inplace=True)
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



########################### DEFINE THE DIRECTION VECTOR ################################

#%%
if 0:
    if apply_twtt_correction_factor:
        DFbig['TWTT_obs']      = DFbig['TWTT'] * 10**-6
        
    else:
        DFbig['TWTT_obs']      = DFbig['TWTT']

    DFbig['TWTT_mod']          = np.zeros(len(DFbig['TWTT']))
    DFbig['TWTT_mod_equiv']    = np.zeros(len(DFbig['TWTT']))
    DFbig['TWTT_mod_sneldesc'] = np.zeros(len(DFbig['TWTT']))

    DFbig['B_TWTT']   = np.zeros(len(DFbig['TWTT']))
    
    DFbig['A_xbea_x'] = np.zeros(len(DFbig['TWTT']))
    DFbig['A_xbea_y'] = np.zeros(len(DFbig['TWTT']))
    DFbig['A_xbea_z'] = np.zeros(len(DFbig['TWTT']))
    DFbig['A_c']      = np.zeros(len(DFbig['TWTT']))
    
    DFbig['Dir_RPY_R'] = np.zeros(len(DFbig['TWTT'])) * np.nan
    DFbig['Dir_RPY_P'] = np.zeros(len(DFbig['TWTT'])) * np.nan
    DFbig['Dir_RPY_Y'] = np.zeros(len(DFbig['TWTT'])) * np.nan
    
    DFbig['Dir_NED_N'] = np.zeros(len(DFbig['TWTT'])) * np.nan
    DFbig['Dir_NED_E'] = np.zeros(len(DFbig['TWTT'])) * np.nan
    DFbig['Dir_NED_D'] = np.zeros(len(DFbig['TWTT'])) * np.nan
    
    for dateemi in DFbig['date_emi'].unique():
        DFepoc  = DFbig[DFbig['date_emi'] == dateemi]   
        
        ############### VLBI Direction ###############
        Dir_RPY = ibpf.direction_vector_finder(DFepoc,
                                               lever_arm_RPY_dic,
                                               cin=Cinp[0])
        
        DFbig.loc[DFepoc.index,'Dir_RPY_R'] = Dir_RPY[0]
        DFbig.loc[DFepoc.index,'Dir_RPY_P'] = Dir_RPY[1]
        DFbig.loc[DFepoc.index,'Dir_RPY_Y'] = Dir_RPY[2]
        
        #Head,Pitch,Roll
        Euler = DFepoc[["head_rec","pitc_rec","roll_rec"]].drop_duplicates().values
        Rot = Rotation.from_euler("zyx",Euler,degrees=False)
        Dir_NED = Rot.apply(Dir_RPY).squeeze()
        Dir_NED = Dir_NED / np.linalg.norm(Dir_NED)
    
        DFbig.loc[DFepoc.index,'Dir_NED_N'] = Dir_NED[0]
        DFbig.loc[DFepoc.index,'Dir_NED_E'] = Dir_NED[1]
        DFbig.loc[DFepoc.index,'Dir_NED_D'] = Dir_NED[2]
        
        dist = np.mean(DFepoc["TWTT_obs"] * 0.5 * c_ref_init)
        DirExt_NED = Dir_NED * dist

        DFbig.loc[DFepoc.index,'DirExt_NED_N'] = DirExt_NED[0]
        DFbig.loc[DFepoc.index,'DirExt_NED_E'] = DirExt_NED[1]
        DFbig.loc[DFepoc.index,'DirExt_NED_D'] = DirExt_NED[2]        

        DFbig.loc[DFepoc.index,'N_BEA_DirExt'] = DirExt_NED[0] + DFbig.loc[DFepoc.index,'N_AHD_rec']
        DFbig.loc[DFepoc.index,'E_BEA_DirExt'] = DirExt_NED[1] + DFbig.loc[DFepoc.index,'E_AHD_rec']
        DFbig.loc[DFepoc.index,'D_BEA_DirExt'] = DirExt_NED[2] + DFbig.loc[DFepoc.index,'D_AHD_rec']       
        
        TWTTsort , IDsort = utils.sort_binom_list(DFepoc["TWTT_raw"], DFepoc["ID_TDC"],True)
        IDsort_str = [ str(e) for e in IDsort ] 
        
        if 0:
            print("*****************")
            for i0,i in enumerate(IDsort[1:]):
                D = np.linalg.norm(lever_arm_RPY_dic["TDC" + str(i)] - lever_arm_RPY_dic["TDC" + str(IDsort[0])])
                T = (TWTTsort[i0+1] - TWTTsort[0]) * 10**-6
                #print("D",D)
                #print("T",T)
                C = D/T
                print("C",C)
            print("*****************")
     
    if False:   
        DFbig[DFbig.ID_BEA == 1].plot('date_rec','Dir_NED_E')
        DFbig[DFbig.ID_BEA == 2].plot('date_rec','Dir_NED_E')
        DFbig[DFbig.ID_BEA == 3].plot('date_rec','Dir_NED_E')
            
        DFbig[DFbig.ID_BEA == 1].plot('date_rec','Dir_NED_N')
        DFbig[DFbig.ID_BEA == 2].plot('date_rec','Dir_NED_N')
        DFbig[DFbig.ID_BEA == 3].plot('date_rec','Dir_NED_N')
        
    DFbig['Azi_NED'] = np.rad2deg(np.arctan2(DFbig['Dir_NED_E'],DFbig['Dir_NED_N']))
    
    
    
    
    
    ################# PLOT AND TEST ZONE
    if False:
        DFbig[DFbig.ID_BEA == 1].plot("date_rec","Azi_NED",style="x") 
        DFbig[DFbig.ID_BEA == 2].plot("date_rec","Azi_NED",style="x") 
        DFbig[DFbig.ID_BEA == 3].plot("date_rec","Azi_NED",style="x") 
        
        DFbig[DFbig.ID_BEA == 1].plot("E_BEA_DirExt","N_BEA_DirExt",style="x") 
        plt.axis("equal")
        DFbig[DFbig.ID_BEA == 2].plot("E_BEA_DirExt","N_BEA_DirExt",style="x") 
        plt.axis("equal")
        DFbig[DFbig.ID_BEA == 3].plot("E_BEA_DirExt","N_BEA_DirExt",style="x") 
        plt.axis("equal")
        
        
        DFbig[DFbig.ID_BEA == 1].plot("date_rec","Azi_NED",style="x") 
        DFbig[DFbig.ID_BEA == 2].plot("date_rec","Azi_NED",style="x") 
        DFbig[DFbig.ID_BEA == 3].plot("date_rec","Azi_NED",style="x")
    
    
    DFbig[["date_rec","E_AHD_emi"]].drop_duplicates().plot("date_rec","E_AHD_emi",style="x")
    
#%% 

########################### DEFINE THE DIRECTION VECTOR ################################




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
        DFdirvect_bea      = DFdirvect_big.copy()
        
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
            
            #### !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            up_offset = 0
            xemi[2] += up_offset
            xrec[2] += up_offset
            #### !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            twtt_obs = rowbea['TWTT_obs']
                              
            ###################### START MODELING & DERIVATION #######################
            
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
            
            dir_vect_obs = rowdirvectbea[["vN_nrm","vE_nrm","vD_nrm"]].values.astype(np.float64)
        
            args_for_partial_dev_dir_vect = (xbea_apri[0],
                                             xbea_apri[1],
                                             xbea_apri[2],
                                             xrec_AHD)      
            
            dir_vect_mod = ibpf.fct_obs_dir_vect(*args_for_partial_dev_dir_vect)
            dir_vect_mod = np.array(dir_vect_mod)
                    
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

    if with_bl:
        K , Q , P = stats.weight_mat([p_twtt,p_bl],[len(B_twtt),len(B_bl_stk)])
    elif with_direction_vectors:
        K , Q , P = stats.weight_mat([p_twtt,p_dir_vect],[len(B_twtt),len(Bdir_vect)])
    else:
        K , Q , P = stats.weight_mat([p_twtt],[len(B_twtt)])        
        
        
            
    
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
    
    Xnew_reshape = Xnew2.reshape(len(ID_bea_list),3)
    BLnew = reffram.BL_from_points(Xnew_reshape)
    ############################# GET NEW VALUES #############################        
    
    
    ############################ STORE NEW VALUES ############################ 
    DictIteraStore_iter = dict()
    DictIteraStore[i_itera] = DictIteraStore_iter

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

    ############################ STORE NEW VALUES ############################ 
    
    
    print("INFO: estim_c           :",estim_c)
    print("INFO: direction_vectors :",with_direction_vectors)
    print("INFO: weight TWTT       :",p_twtt)
    print("INFO: weight Dir Vector :",p_dir_vect)

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

        
    ######### COMPUTE ANGLE FRIENDLY DIRECVTION VECTOR AND OUTLIERS #########
    



             

    print("Xapri    ",i_itera,DictIteraStore[i_itera]["Xapri"])
    print("Xprev    ",i_itera,DictIteraStore[i_itera]["Xprev"])
    print("Xnew     ",i_itera,DictIteraStore[i_itera]["Xnew"])
    print("dX       ",i_itera,DictIteraStore[i_itera]["dX"])
    print("B        ",i_itera,DictIteraStore[i_itera]["B"])
    print("V        ",i_itera,DictIteraStore[i_itera]["V"])
    print("Vnat     ",i_itera,DictIteraStore[i_itera]["Vnatural"])
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
    
    
print("***********************************************************************")
    

# for i in range(n_itera):
#     print("B   ",i,DictIteraStore[i]["B"])
#     print("V   ",i,DictIteraStore[i]["V"])
#     print("Vnat",i,DictIteraStore[i]["Vnatural"])
#     print("S B**2   ",i,np.sum(DictIteraStore[i]["B"]**2))
#     print("S V**2   ",i,np.sum(DictIteraStore[i]["V"]**2))
#     print("S Vnat**2",i,np.sum(DictIteraStore[i]["Vnatural"]**2))



#%%
###################### PLOT FOR TRAJECTORY
if False:
    ############ Trajectory EN
    fig,ax = plt.subplots()
    ax.plot(DFbig["E_AHD_emi"],DFbig["N_AHD_emi"],"+")
    ax.plot(DFbig["E_TDC_rec"],DFbig["N_TDC_rec"],"x")
    DFtdc1 = DFbig[DFbig.ID_TDC == 1]
    ax.plot(DFtdc1["E_TDC_rec"],DFtdc1["N_TDC_rec"],".")
    ax.axis("equal")
    
    ############ Trajectory Time Series        
    figt,(axte,axtn) = plt.subplots(2,1)
    figt.suptitle("Trajectory TimeSeries")
    
    axte.plot(DFbig["date_emi"],DFbig["E_AHD_emi"],"+")
    axtn.plot(DFbig["date_emi"],DFbig["N_AHD_emi"],"+")
    axte.plot(DFbig["date_rec"],DFbig["E_TDC_emi"],"x")
    axtn.plot(DFbig["date_rec"],DFbig["N_TDC_emi"],"x")

    figttdc,(axttdce,axttdcn) = plt.subplots(2,1)
    figttdc.suptitle("Trajectory TimeSeries per transducer")

    for tdc in DFbig.ID_TDC.unique():
        DFtsc = DFbig[DFbig.ID_TDC == tdc]
        tdcstr = "TDC" + str(tdc)
        cols_posi_rec = ['N_' +tdcstr+ '_rec',
                         'E_' +tdcstr+ '_rec',
                         'D_' +tdcstr+ '_rec']
        axttdce.plot(DFtsc["date_rec"],DFtsc[cols_posi_rec[1]],"x")
        axttdcn.plot(DFtsc["date_rec"],DFtsc[cols_posi_rec[0]],"x")
    
    figup,axup = plt.subplots()
    axup.plot(DFbig["date_rec"],DFbig["D_AHD_emi"],"x")
    figup.suptitle("Up + TWTT")
    
    axtwtt = axup.twinx()
    axtwtt.plot(DFbig["date_rec"],DFbig["TWTT_obs"],"+")
    
    fighistangl, axhistangl = plt.subplots()
    axhistangl.hist(Angle_residuals,100)
    fighistangl.suptitle("Histogram of Direction cosine angle")
        
    F.stop()
    
###################### PLOT FOR TRAJECTORY

    
###################### PLOT FOR RESIDUALS
figres,axres   = plt.subplots()
fighist,axhist = plt.subplots()

fighist.suptitle("Histogram of the residuals")

figres.suptitle("TWTT residuals")
for idbea in DFbig.ID_BEA.unique():    
    DFbeaplt = DFbig[DFbig.ID_BEA == idbea] 
    axres.plot(DFbeaplt["date_rec"],DFbeaplt["B_TWTT"],"x")
    DFbeaplt_valid = DFbeaplt[DFbeaplt["VALID"]]
    axres.plot(DFbeaplt_valid["date_rec"],DFbeaplt_valid["B_TWTT"],"+")
    
    

axhist.hist(DFbeaplt_valid["B_TWTT"],50,density=False,label="last iter")
axhist.hist(DictIteraStore[0]["B"],50,density=False,label="1st iter")
axhist.legend()
fighist.savefig(os.path.dirname(p_obs) + "/hist_res_TWTT_angle.png")





#%%





