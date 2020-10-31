#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 15:18:25 2020

@author: psakicki

mk04 : retour a des Arrays, trop lent sinon
mk06 : improvement of mk04, now the BL can ben switched off and then only one Beacon
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

##############################################################################

def equiv_straght_ray(Zinp,Cinp,zrec,
                      extrapolate_first = True,
                      harmonic_mean = True):
    
    """
    harmonic_mean = False => Aritmetic Mean
    """
    
    
    if extrapolate_first:
        Z,C = ssp.SSP_extrapolate(Zinp, Cinp, 10000, 1)
    else:
        Z,C = Zinp,Cinp
    
    Z , C = ssp.SSP_cut(Z, C, zrec)
    
    if harmonic_mean:
        #### Harmonic Mean
        m = (Z[-1] - Z[0]) * ((1/np.trapz(1/C,Z)))
    else:
        #### Aritmetic Mean
        m = (1/(Z[-1] - Z[0])) * np.trapz(C,Z)
        
    return m

def fct_obs_equiv(xbea_x_apri_in,
                  xbea_y_apri_in,
                  xbea_z_apri_in,
                  xemi_in,
                  xrec_in,
                  c_bea_in):
    
    xbea_apri_in =  np.array([xbea_x_apri_in,
                              xbea_y_apri_in,
                              xbea_z_apri_in])
    
    d_emi = conv.dist(xbea_apri_in,xemi_in)
    d_rec = conv.dist(xbea_apri_in,xrec_in)
    
    #t_emi = d_emi / c_bea_in
    #t_rec = d_rec / c_bea_in
    
    t_sum = (d_emi + d_rec) / c_bea_in
    
    return t_sum


def fct_obs_equiv_deriv_ana(xbea_x_apri_in,
                            xbea_y_apri_in,
                            xbea_z_apri_in,
                            xemi_in,
                            xrec_in,
                            c_bea_in):
    
    # xbea_apri_in =  np.array([xbea_x_apri_in,
    #                           xbea_y_apri_in,
    #                           xbea_z_apri_in])
    
    # d_emi = conv.dist(xbea_apri_in,xemi_in)
    # d_rec = conv.dist(xbea_apri_in,xrec_in)
    
    # keof_denom = 1
    # keof_numer = 1
    
    # denom = (keof_denom*c_bea_in*(d_emi + d_rec))
    
    # diff_xbea_x = (keof_numer*xemi_in[0] - keof_numer*xbea_x_apri_in + keof_numer*xrec_in[0] - keof_numer*xbea_x_apri_in)/denom
    # diff_xbea_y = (keof_numer*xemi_in[1] - keof_numer*xbea_y_apri_in + keof_numer*xrec_in[1] - keof_numer*xbea_y_apri_in)/denom
    # diff_xbea_z = (keof_numer*xemi_in[2] - keof_numer*xbea_z_apri_in + keof_numer*xrec_in[2] - keof_numer*xbea_z_apri_in)/denom
    
    # diff_xbea_x = (keof_numer*xbea_x_apri_in - keof_numer*xemi_in[0] + keof_numer*xbea_x_apri_in - keof_numer*xrec_in[0])/denom
    # diff_xbea_y = (keof_numer*xbea_y_apri_in - keof_numer*xemi_in[1] + keof_numer*xbea_y_apri_in - keof_numer*xrec_in[1])/denom
    # diff_xbea_z = (keof_numer*xbea_z_apri_in - keof_numer*xemi_in[2] + keof_numer*xbea_z_apri_in - keof_numer*xrec_in[2])/denom
    
    # diff_c = -1 * ((d_emi + d_rec)/(c_bea_in**2))
    
    xb,yb,zb =  xbea_x_apri_in,xbea_y_apri_in,xbea_z_apri_in
    xe,ye,ze = xemi_in
    xr,yr,zr = xrec_in
    c = c_bea_in
    
    ##### Dirty copy/paste from Maple        
    diff_xbea_x = (((xb - xe) ** 2 + (yb - ye) ** 2 + (zb - ze) ** 2) ** (-0.1e1 / 0.2e1) * (2 * xb - 2 * xe) / 2 + ((xb - xr) ** 2 + (yb - yr) ** 2 + (zb - zr) ** 2) ** (-0.1e1 / 0.2e1) * (2 * xb - 2 * xr) / 2) / c
    diff_xbea_y = (((xb - xe) ** 2 + (yb - ye) ** 2 + (zb - ze) ** 2) ** (-0.1e1 / 0.2e1) * (2 * yb - 2 * ye) / 2 + ((xb - xr) ** 2 + (yb - yr) ** 2 + (zb - zr) ** 2) ** (-0.1e1 / 0.2e1) * (2 * yb - 2 * yr) / 2) / c
    diff_xbea_z = (((xb - xe) ** 2 + (yb - ye) ** 2 + (zb - ze) ** 2) ** (-0.1e1 / 0.2e1) * (2 * zb - 2 * ze) / 2 + ((xb - xr) ** 2 + (yb - yr) ** 2 + (zb - zr) ** 2) ** (-0.1e1 / 0.2e1) * (2 * zb - 2 * zr) / 2) / c
    
    diff_c = -(math.sqrt((xb - xe) ** 2 + (yb - ye) ** 2 + (zb - ze) ** 2) + math.sqrt((xb - xr) ** 2 + (yb - yr) ** 2 + (zb - zr) ** 2)) / c ** 2

    return diff_xbea_x , diff_xbea_y , diff_xbea_z , diff_c

import sympy

def fct_obs_equiv_deriv_symbo(xbea_x_apri_in,
                              xbea_y_apri_in,
                              xbea_z_apri_in,
                              xemi_in,
                              xrec_in,
                              c_bea_in):
    
    xbea,ybea,zbea = sympy.symbols('xbea ybea zbea')
    xemi,yemi,zemi = sympy.symbols('xemi yemi zemi')
    xrec,yrec,zrec = sympy.symbols('xrec yrec zrec')
    c = sympy.symbols('c')
    
    demi = sympy.sqrt((xbea-xemi)**2 + (ybea-yemi)**2 + (zbea-zemi)**2)
    drec = sympy.sqrt((xbea-xrec)**2 + (ybea-yrec)**2 + (zbea-zrec)**2)
    
    t = (demi + drec) / c
    
    tdiff_c = sympy.diff(t,c)
    
    tdiff_xbea = sympy.diff(t,xbea)
    tdiff_ybea = sympy.diff(t,ybea)
    tdiff_zbea = sympy.diff(t,zbea)
    
    dict_for_eval = dict()    

    dict_for_eval['xbea']=xbea_x_apri_in
    dict_for_eval['ybea']=xbea_y_apri_in
    dict_for_eval['zbea']=xbea_z_apri_in
    dict_for_eval['xemi']=xemi_in[0]
    dict_for_eval['yemi']=xemi_in[1]
    dict_for_eval['zemi']=xemi_in[2]
    dict_for_eval['xrec']=xrec_in[0]
    dict_for_eval['yrec']=xrec_in[1]
    dict_for_eval['zrec']=xrec_in[2]
    dict_for_eval['c']=c_bea_in
    
    tdiff_xbea_eva = tdiff_xbea.evalf(subs=dict_for_eval)
    tdiff_ybea_eva = tdiff_ybea.evalf(subs=dict_for_eval)
    tdiff_zbea_eva = tdiff_zbea.evalf(subs=dict_for_eval)
    tdiff_c_eva    = tdiff_c.evalf(subs=dict_for_eval)
    
    return tdiff_xbea_eva,tdiff_ybea_eva,tdiff_zbea_eva,tdiff_c_eva


def fct_obs_vlbi(B_in,S_in,c_in):
    dt = -1 * (np.dot(B_in,S_in)) / c_in
    return dt

def fct_obs_vlbi_deriv_symbo(B_in,S_in,c_in):

    bx,by,bz = sympy.symbols('bx by bz')    
    sx,sy,sz = sympy.symbols('sx sy sz')
    c = sympy.symbols('c')
    
    dt = -1 * (bx*sx + by*sy + bz*sz) / c
    
    dt_diff_sx = sympy.diff(dt,sx)
    dt_diff_sy = sympy.diff(dt,sy)
    dt_diff_sz = sympy.diff(dt,sz)
    
    dict_for_eval = dict()    

    dict_for_eval['bx']=B_in[0]
    dict_for_eval['by']=B_in[1]
    dict_for_eval['bz']=B_in[2]
    dict_for_eval['sx']=S_in[0]
    dict_for_eval['sy']=S_in[1]
    dict_for_eval['sz']=S_in[2]
    dict_for_eval['c']=c_in
    
    dt_diff_sx_eva = np.float64(dt_diff_sx.evalf(subs=dict_for_eval))
    dt_diff_sy_eva = np.float64(dt_diff_sy.evalf(subs=dict_for_eval))
    dt_diff_sz_eva = np.float64(dt_diff_sz.evalf(subs=dict_for_eval))
            
    return dt_diff_sx_eva,dt_diff_sy_eva,dt_diff_sz_eva


def fct_obs_vlbi_deriv_ana(B_in,S_in,c_in):
    
    dt_diff_sx_eva = - B_in[0] / c_in
    dt_diff_sy_eva = - B_in[1] / c_in
    dt_diff_sz_eva = - B_in[2] / c_in
        
    return dt_diff_sx_eva,dt_diff_sy_eva,dt_diff_sz_eva



def direction_vector_finder(DFepoc,
                            lever_arm_RPY_dic,
                            id_tdc_ref = 1,
                            Sin = np.array([0,0,0]),
                            cin=1500):
    
        if len(DFepoc) != 4:
            return np.ones(3) * np.nan
        
        B = []
        Aline_stk = []        
        
        DFepoc_ref = DFepoc[DFepoc['ID_TDC'] == id_tdc_ref]
        t_ref  = DFepoc_ref['TWTT_obs'].values[0]    
        
        for ilin , lin in DFepoc.iterrows():
            
            id_tdc = lin['ID_TDC']
            
            if id_tdc == id_tdc_ref:
                continue

            t_tdc  = lin['TWTT_obs']

            dt_obs = t_tdc - t_ref
            Baseline_in = lever_arm_RPY_dic["TDC" + str(id_tdc)] - lever_arm_RPY_dic["TDC" + str(id_tdc_ref)] 
        
            dt_mod = fct_obs_vlbi(Baseline_in,Sin,cin)
            
            B.append(dt_obs - dt_mod)
            
            #Aline  = fct_obs_vlbi_deriv_symbo(Baseline_in,Sin,cin)
            Aline = fct_obs_vlbi_deriv_ana(Baseline_in,Sin,cin)

            #Aline_stk.append(Aline)
            Aline_stk.append(Aline)
                    
        A  = np.vstack(Aline_stk)
        At = np.transpose(A)
        B  = np.array(B)
        
        N = At.dot(A)
        Ninv = scipy.linalg.inv(N)
        
        dX = Ninv.dot(At).dot(B)
        
        B - A.dot(dX)
        
        Snew = Sin + dX
        
        Snew = Snew / np.linalg.norm(Snew)
        
        return Snew
    
##############################################################################




############################ LOADING DATA ####################################
check_res = 1

############## SIMU
pC = "/home/psakicki/CODES/acoustyx_toolbox_2_py3/working/SIMU_SUMMER2020_1000_R50_nois1-0.0____/SIMU_SUMMER2020_1000_R50_nois1-0.0_____.C.dat"
pZ = "/home/psakicki/CODES/acoustyx_toolbox_2_py3/working/SIMU_SUMMER2020_1000_R50_nois1-0.0____/SIMU_SUMMER2020_1000_R50_nois1-0.0_____.Z.dat"
p_obs = "/home/psakicki/CODES/acoustyx_toolbox_2_py3/working/SIMU_SUMMER2020_1000_R50_nois1-0.0____/SIMU_SUMMER2020_1000_R50_nois1-0.0_____.PXP1.P.csv"
p_obs = "/home/psakicki/CODES/acoustyx_toolbox_2_py3/working/SIMU_SUMMER2020_3x333_x500_y500_ang0_nois1-0.0____/DF_new_format.csv"
p_bl = "/home/psakicki/CODES/acoustyx_toolbox_2_py3/working/SIMU_SUMMER2020_1000_R50_nois1-0.0____/SIMU_SUMMER2020_1000_R50_nois1-0.0_____.B.csv"

############## REAL
pZ = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/01_/PAMELI_BREST_vJupyter_2_4transducers/PAMELI_BREST_vJupyter_2_4transducers.Z.dat"
pC = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/01_/PAMELI_BREST_vJupyter_2_4transducers/PAMELI_BREST_vJupyter_2_4transducers.C.dat"
### Corrections in lever arms
p_obs = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/01_/PAMELI_BREST_vJupyter_5_4transducers_new_lever_arms/PAMELI_BREST_vJupyter_5_4transducers_new_lever_arms.csv"
p_bl  = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/01_/PAMELI_BREST_vJupyter_2_4transducers/PAMELI_BREST_vJupyter_2_4transducers_BASELINE.csv"
p_obs = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/01_/PAMELI_BREST_vJupyter_5_4transducers_new_lever_arms_GINS_PPP_coords/PAMELI_BREST_vJupyter_5_4transducers_new_lever_arms.csv"
p_obs = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/01_/PAMELI_BREST_vJupyter_7_4transducers_newRF_newDTime_BAD_BEFORE_TIMING_CORR/PAMELI_BREST_vJupyter_7_4transducers_newRF_newDTime_BAD_BEFORE_TIMING_CORR.csv"
p_obs = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/01_/PAMELI_BREST_vJupyter_6_4transducers_newRF_newDTime/PAMELI_BREST_vJupyter_6_4transducers_newRF_newDTime.csv"




p_att = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/01_/PAMELI_BREST_vJupyter_5_4transducers_new_lever_arms/PAMELI_BREST_vJupyter_5_4transducers_new_lever_arms_AttitudeInterpolator.pik"

Zinp = np.loadtxt(pZ) 
Cinp = np.loadtxt(pC)

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

DFbig['date_rec'].apply(conv.string_date2dt)

DFbl = pd.read_csv(p_bl,index_col=0)
DFbl = DFbl.drop_duplicates('dist')
DFbl = DFbl[DFbl.ID_BEA_A != DFbl.ID_BEA_B]
DFbl.reset_index(inplace=True)

############################ LOADING DATA ####################################


############################ SET ENVIRONEMENT VARIABLE #######################
cols_posi = ['N_AHD_emi', 'E_AHD_emi', 'D_AHD_emi',
             'N_AHD_rec', 'E_AHD_rec', 'D_AHD_rec']

cols_posi_emi = cols_posi[0:3]
cols_posi_rec = cols_posi[3:6]

ID_bea_list_orig = np.sort(DFbig['ID_BEA'].unique())

ID_bea_list = ID_bea_list_orig
ID_bea_list = [3,4]
ID_bea_list = [2]

n_itera = 6

DictIteraStore = dict()

p_twtt = 10**0
#####p_c    = 10**0
p_bl   = 10**-5

with_numerical_diff = 0
with_bl = False

fast_mode         = False ############ NOT PROPERLY IMPLEMENTED
snelldesc_compute = True
snelldesc_use     = False
estim_c           = False
apply_twtt_correction_factor = True

############################ SET ENVIRONEMENT VARIABLE #######################
       

############################ SET INITIAL APRI ################################
c_ref_init = equiv_straght_ray(Zinp, Cinp, 40)


Xbea_apri_all_init = np.array([[7500  , 7500  , 4010],
                               [12500 , 12500 , 3980],
                               [ 7500 , 12500 , 4030],
                               [ 12500 , 7500 , 3960]])

Xbea_apri_all_init = np.array([[7500  , 7500  , 4010],
                               [12500 , 12500 , 3980],
                               [ 12500 , 7500 , 4030],
                               [ 7500  , 12500, 3960]])


Xbea_apri_all_init = np.array([[24.31494786,   -7.20555817, 40.07],
                               [-26.20205757,  -7.53913662, 40.66],
                               [2.23923942,    22.42808402, 39.68]])
Xbea_apri_all_init = np.array([[   -7.20555817, 24.31494786, 40.07],
                               [   -7.53913662,-26.20205757, 40.66],
                               [  22.42808402,  2.23923942,  39.68]])

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
if 1:
    if apply_twtt_correction_factor:
        DFbig['TWTT_obs']          = DFbig['TWTT'] * 10**-6
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
    
    for daterec in DFbig['date_rec'].unique():
        DFepoc  = DFbig[DFbig['date_rec'] == daterec]   
        Dir_RPY = direction_vector_finder(DFepoc,lever_arm_RPY_dic,
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
        
        TWTTsort , IDsort = utils.sort_binom_list(DFepoc["TWTTraw"], DFepoc["ID_TDC"],True)
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

############################ SET APRI FOR THE 1st ITERA ######################

for i_itera in range(n_itera):
    
    print("INFO:","iteration",i_itera)

    if apply_twtt_correction_factor:
        DFbig['TWTT_obs']          = DFbig['TWTT'] * 10**-6
    DFbig['TWTT_mod']          = np.zeros(len(DFbig['TWTT']))
    DFbig['TWTT_mod_equiv']    = np.zeros(len(DFbig['TWTT']))
    DFbig['TWTT_mod_sneldesc'] = np.zeros(len(DFbig['TWTT']))

    DFbig['B_TWTT']   = np.zeros(len(DFbig['TWTT']))
    
    DFbig['A_xbea_x'] = np.zeros(len(DFbig['TWTT']))
    DFbig['A_xbea_y'] = np.zeros(len(DFbig['TWTT']))
    DFbig['A_xbea_z'] = np.zeros(len(DFbig['TWTT']))
    DFbig['A_c']      = np.zeros(len(DFbig['TWTT']))
    Xbea_apri_use     = []   

    

    for i_idbea, idbea in enumerate(ID_bea_list):
    
        print("INFO:","partial derivative for beacon",idbea)
        
        DFbea_orig = DFbig[DFbig['ID_BEA'] == idbea]
        DFbea = DFbea_orig.copy()
        
        print("INFO: ping valids",np.sum(DFbea['VALID']))
        
        xbea_apri = Xbea_apri_all[i_idbea-1]
        Xbea_apri_use.append(xbea_apri)
        
        c_bea = equiv_straght_ray(Zinp, Cinp, xbea_apri[2])
        c_used = c_ref
         
        Jac_xbea_stk = []
        Jac_c_stk    = []
        
        iiiiii=0
        Afast_stk = []

        for i_rowbea , rowbea in DFbea.iterrows():
            
            if not rowbea['VALID']:
                continue
            
            iiiiii+=1
            
            xemi = np.array(rowbea[cols_posi_emi].values,dtype=float)
            xrec = np.array(rowbea[cols_posi_rec].values,dtype=float)
            
            twtt_obs = rowbea['TWTT_obs']
                  
            ###################### MODELING & DERIVATION #######################
            
            args_for_partial_dev =  (xbea_apri[0],
                                     xbea_apri[1],
                                     xbea_apri[2],
                                     xemi,
                                     xrec,
                                     c_used)
            
            ########## Modeled observations
            twtt_mod = fct_obs_equiv(*args_for_partial_dev)
            
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

            rowbea['TWTT_mod_equiv']    = twtt_mod
            rowbea['TWTT_mod_sneldesc'] = twtt_mod_sneldesc_new[-1]
            
            if snelldesc_use:
                rowbea['TWTT_mod']          = rowbea['TWTT_mod_sneldesc'] 
            else:
                rowbea['TWTT_mod']          = rowbea['TWTT_mod_equiv'] 

            
            rowbea['B_TWTT'] = rowbea['TWTT_obs'] - rowbea['TWTT_mod']

            
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
                jac_ana_line = np.array(fct_obs_equiv_deriv_ana(*args_for_partial_dev))
                
                if not fast_mode:
                    rowbea[['A_xbea_x','A_xbea_y','A_xbea_z','A_c']] = jac_ana_line[0],jac_ana_line[1],jac_ana_line[2],jac_ana_line[3]
                else:
                    Afast_row = [jac_ana_line[0],
                                 jac_ana_line[1],
                                 jac_ana_line[2],
                                 jac_ana_line[3]]
                    
                    Afast_stk.append(Afast_row)
                    
            ######### Symbolic derivation
            #jac_sym = fct_obs_equiv_deriv_symbo(*args_for_partial_dev)
                        
            ####### Update the row
            if not fast_mode:
                DFbea.loc[i_rowbea] = rowbea
                
        if fast_mode:
            Afast = np.array(Afast_stk)
            DFbig[(DFbig["ID_BEA"] == idbea) & (DFbig["VALID"] == True)][['A_xbea_x','A_xbea_y','A_xbea_z','A_c']] = Afast
            
        ###################### MODELING & DERIVATION #######################               
        
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
            
            xbea_A = Xbea_apri_all[int(rowbl.ID_BEA_A)-1]
            xbea_B = Xbea_apri_all[int(rowbl.ID_BEA_B)-1]
            
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
    
    elif not estim_c and with_bl:
        A = np.vstack((A1,np.vstack(A_bl_stk)))
        
    elif not estim_c and not with_bl:
        A = A1
    

    
    B_twtt = np.hstack(B_twtt_stk)
    
    if with_bl:
        B = np.hstack((B_twtt,B_bl_stk))
    else:
        B = B_twtt
        
    
    
    if with_bl:
        K , Q , P = stats.weight_mat([p_twtt,p_bl],[len(B_twtt),len(B_bl_stk)])
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
    
    Xnew_reshape = Xnew.reshape(len(ID_bea_list),3)
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
    DictIteraStore_iter['Xprev']             = Xprev
    DictIteraStore_iter['Xnew']              = Xnew
    DictIteraStore_iter['Xnew_reshape']      = Xnew_reshape
    DictIteraStore_iter['BLnew']             = BLnew
    DictIteraStore_iter['BLobs']             = DFbl["dist"]

    ############################ STORE NEW VALUES ############################ 
    
    ####################### COMPUTE NATURAL RESIDUALS ########################
    B_TWTT_Stk = []
    for i_idbea, idbea in enumerate(ID_bea_list):
    
        print("INFO:","Natural residuals for beacon",idbea)
        
        DFbea = DFbig[(DFbig['ID_BEA'] == idbea) & (DFbig['VALID'])]
        
        xbea_apri = Xbea_apri_all_new[i_idbea-1]
        
        Xemi     = DFbea[cols_posi_emi].values
        Xrec     = DFbea[cols_posi_rec].values
        TWTT_obs = DFbea['TWTT'].values * 10**-6
        
         
        TWTT_mod_stk = []
        RESstk = []   
        
        c_bea = equiv_straght_ray(Zinp, Cinp, xbea_apri[2])
        
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
            twtt_mod = fct_obs_equiv(*args_for_partial_dev)
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
    c_ref = c_ref_new
    ##################### REPLACE APRIORI WITH NEW VALUES ####################                


    print("Xprev    ",i_itera,DictIteraStore[i_itera]["Xprev"])
    print("Xnew     ",i_itera,DictIteraStore[i_itera]["Xnew"])
    print("dX       ",i_itera,DictIteraStore[i_itera]["dX"])
    print("B        ",i_itera,DictIteraStore[i_itera]["B"])
    print("V        ",i_itera,DictIteraStore[i_itera]["V"])
    print("Vnat     ",i_itera,DictIteraStore[i_itera]["Vnatural"])
    print("S B**2   ",i_itera,np.sum(DictIteraStore[i_itera]["B"]**2))
    print("S V**2   ",i_itera,np.sum(DictIteraStore[i_itera]["V"]**2))
    print("S Vnat**2",i_itera,np.sum(DictIteraStore[i_itera]["Vnatural"]**2))
    print("BL new   ",i_itera,acls.BL_from_PXPlist(Xbea_apri_all_new))
    print("BL obs   ",i_itera,DictIteraStore_iter['BLobs'])
    print("Barycntr ",i_itera,acls.barycenter_calc(Xbea_apri_all_new))
    print("c_ref_old",i_itera,DictIteraStore_iter['c_ref'] )
    print("c_ref_new",i_itera,DictIteraStore_iter['c_ref_new'] )
    
    

for i in range(n_itera):
    print("B   ",i,DictIteraStore[i]["B"])
    print("V   ",i,DictIteraStore[i]["V"])
    print("Vnat",i,DictIteraStore[i]["Vnatural"])
    print("S B**2   ",i,np.sum(DictIteraStore[i]["B"]**2))
    print("S V**2   ",i,np.sum(DictIteraStore[i]["V"]**2))
    print("S Vnat**2",i,np.sum(DictIteraStore[i]["Vnatural"]**2))



if False:
    plt.clf()
    for idbea in ID_bea_list:
        DFbea = DFbig[DFbig.ID_BEA == idbea]
        BBBB  = DFbea["TWTT_mod_sneldesc"] - DFbea["TWTT_mod"]
        CCCC  = DFbea["TWTT"] * 10**-6     - DFbea["TWTT_mod"]
        DDDD  = DFbea["TWTT"] * 10**-6     - DFbea["TWTT_mod_sneldesc"]
        
        
    
        #plt.plot(BBBB,".")
        plt.plot(CCCC,"+")
        #plt.plot(DDDD,"x")


# plt.figure()

# plt.plot(V,label="V")
# plt.plot(Vnatural,label="Vnat")
# plt.legend()

# Diff1 = V - Vnatural
# Diff2 = V - np.flip(Vnatural)

# plt.figure()

# plt.plot(Diff1)
# plt.plot(Diff2)



#%%
if False:
    plt.clf()
    
    # plt.plot(DFbig['date_rec'].values , DFbig['N_GPS_rec'].values,".")
    
    DFbig.plot('date_rec','N_GPS_rec',style='.')
    DFbig.plot('date_rec','X_GPS_rec',style='.')
    DFbig.plot('date_emi','X_GPS_emi',style='r.')
    
    
    
    plt.figure()
    DFplot = DFbig[(DFbig["ID_BEA"] == 3) & (DFbig["VALID"])]
    
    Xplt = DFplot["date_rec"].values
    Yplt = DFplot["TWTT_obs"]
    
    plt.plot(Xplt,Yplt,"x")

if True:
    plt.clf()
    for idbea in ID_bea_list:
        DFbea = DFbig[DFbig.ID_BEA == idbea]
        Y  = DFbea["TWTT"] * 10**-6 * 100
        X  = DFbea["date_rec"].apply(conv.string_date2dt)
        
        N = DFbea["N_GPS_rec"]
        E = DFbea["E_GPS_rec"]
        U = DFbea["U_GPS_rec"]
        
        plt.plot(X,Y,"x")
        
        if idbea == 1:
            plt.plot(X,N,"+-")
            plt.plot(X,E,"+-")
            plt.plot(X,U,"+-")
    
    plt.clf()
    
    # plt.plot(DFbig['date_rec'].values , DFbig['N_GPS_rec'].values,".")
    
    DFbig.plot('date_rec','N_GPS_rec',style='.')
    DFbig.plot('date_rec','X_GPS_rec',style='.')
    DFbig.plot('date_emi','X_GPS_emi',style='r.')
    
    plt.figure()
    DFplot = DFbig[(DFbig["ID_BEA"] == 3) & (DFbig["VALID"])]
    
    Xplt = conv.string_date2dt(DFplot["date_rec"])
    Yplt = DFplot["TWTT_obs"]
    
    plt.suptitle("TWTT_obs")    
    plt.plot(Xplt,Yplt,"x")
    

plt.figure()
plt.suptitle("TWTT for bea 1")
DFBEA = DFbig[DFbig.ID_BEA == 1]
plt.plot(conv.string_date2dt(DFBEA["date_rec"]),DFBEA["TWTT"],".")
