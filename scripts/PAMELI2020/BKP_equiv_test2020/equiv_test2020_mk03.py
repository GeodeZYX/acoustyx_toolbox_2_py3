#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 15:18:25 2020

@author: psakicki
"""

#### Import star style
from geodezyx import *                   # Import the GeodeZYX modules
from geodezyx.externlib import *         # Import the external modules
from geodezyx.megalib.megalib import *   # Import the legacy modules names

from geodezyx import conv
from geodezyx import stats

import SSP as ssp
import raytrace as rt

import numpy as np
import pandas as pd


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


##############################################################################




############################ LOADING DATA ####################################
check_res = 1

pZ = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/01_/PAMELI_BREST_vJupyter_2_4transducers/PAMELI_BREST_vJupyter_2_4transducers.Z.dat"
pC = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/01_/PAMELI_BREST_vJupyter_2_4transducers/PAMELI_BREST_vJupyter_2_4transducers.C.dat"
Zinp = np.loadtxt(pZ) 
Cinp = np.loadtxt(pC)

# SSPpath = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/01_DATA/07-2019_IXBLUE_BREST/0_RAW_DATA/data-gaps/data/CEL_2019-07-24_08h00_utc.txt"
# SSPpath = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/01_DATA/07-2019_IXBLUE_BREST/0_RAW_DATA/data-gaps/data/CEL_2019-07-24_11h30_utc.txt"
# SSP = np.loadtxt(SSPpath)
# Zinp = SSP[:,0]
# Cinp = SSP[:,1]



RS = np.random.RandomState(42)
Cfake = Cinp + RS.rand(len(Cinp)) * 10**-5


### Corrections in lever arms
p_obs = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/01_/PAMELI_BREST_vJupyter_5_4transducers_new_lever_arms/PAMELI_BREST_vJupyter_5_4transducers_new_lever_arms.csv"
### bugs in lever arms
p_obs = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/01_/PAMELI_BREST_vJupyter_2_4transducers/PAMELI_BREST_vJupyter_2_4transducers.csv"

DFbig = pd.read_csv(p_obs,index_col=0)
DFbig_orig = DFbig.copy()
DFbig = DFbig.dropna()
DFbig['VALID'] = np.ones(len(DFbig['TWTT'])).astype(bool)

p_bl = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/01_/PAMELI_BREST_vJupyter_2_4transducers/PAMELI_BREST_vJupyter_2_4transducers_BASELINE.csv"
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

ID_bea_list = [1,2]
ID_bea_list = ID_bea_list_orig


n_itera = 6

DictIteraStore = dict()

p_twtt = 10**0
p_c    = 10**-3
p_bl   = 10**-6

with_numerical_diff = 0

############################ SET ENVIRONEMENT VARIABLE #######################


############################ SET INITIAL APRI ################################
c_ref_init = equiv_straght_ray(Zinp, Cinp, 40)
Xbea_apri_all_init = np.array([[24.31494786,   -7.20555817, 40.07],
                               [-26.20205757,  -7.53913662, 40.66],
                               [2.23923942,    22.42808402, 39.68]])
############################ SET INITIAL APRI ################################


############################ SET APRI FOR THE 1st ITERA ######################

c_ref = c_ref_init
Xbea_apri_all = Xbea_apri_all_init.copy()

############################ SET APRI FOR THE 1st ITERA ######################



for i_itera in range(n_itera):
    
    print("INFO:","iteration",i_itera)

    DFbig['TWTT_obs']          = DFbig['TWTT'] * 10**-6
    DFbig['TWTT_mod']          = np.zeros(len(DFbig['TWTT']))
    DFbig['TWTT_mod_sneldesc'] = np.zeros(len(DFbig['TWTT']))

    DFbig['B_TWTT']   = np.zeros(len(DFbig['TWTT']))
    
    DFbig['A_xbea_x'] = np.zeros(len(DFbig['TWTT']))
    DFbig['A_xbea_y'] = np.zeros(len(DFbig['TWTT']))
    DFbig['A_xbea_z'] = np.zeros(len(DFbig['TWTT']))
    DFbig['A_c']      = np.zeros(len(DFbig['TWTT']))

    for i_idbea, idbea in enumerate(ID_bea_list):
    
        print("INFO:","partial derivative for beacon",idbea)
        
        DFbea_orig = DFbig[DFbig['ID_BEA'] == idbea]
        DFbea = DFbea_orig.copy()
        
        print("INFO: ping valids",np.sum(DFbea['VALID']))
        
        xbea_apri = Xbea_apri_all[idbea-1]
        
        c_bea = equiv_straght_ray(Zinp, Cinp, xbea_apri[2])
        c_used = c_ref
         
        Jac_xbea_stk = []
        Jac_c_stk    = []

        for i_rowbea , rowbea in DFbea.iterrows():
            
            if not rowbea['VALID']:
                continue
            
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

            rowbea['B_TWTT'] = twtt_obs - twtt_mod
            
            # twtt_mod_sneldesc_leg = rt.raytrace_seek((xemi,xrec),
            #                                          xbea_apri,
            #                                          Zinp,
            #                                          Cfake,
            #                                          legacy=True,
            #                                          verbose=False,
            #                                          fulloutput=False)
            
            twtt_mod_sneldesc_new = rt.raytrace_seek((xemi,xrec),
                                                     xbea_apri,
                                                     Zinp,
                                                     Cfake,
                                                     legacy=False,
                                                     verbose=False,
                                                     fulloutput=False)
            
            rowbea['TWTT_mod']          = twtt_mod
            rowbea['TWTT_mod_sneldesc'] = twtt_mod_sneldesc_new[-1]
            
            
            
            ########## derivative of position
            if with_numerical_diff:
                ######### Numerical derivation
                jac_xbea_line = []
                for i in range(3):
                    jac_xbea = stats.partial_derive(fct_obs_equiv, i,
                                                    args_f = args_for_partial_dev)
                    jac_xbea_line.append(jac_xbea)
                    
                rowbea['A_xbea_x'] = jac_xbea_line[0]
                rowbea['A_xbea_y'] = jac_xbea_line[1]
                rowbea['A_xbea_z'] = jac_xbea_line[2]
                            
                ########## derivative of c
                jac_c = stats.partial_derive(fct_obs_equiv, 'c_bea_in',
                                             args_f = args_for_partial_dev)    
                rowbea['A_c'] = jac_c
            else:
                ######### Analytic derivation
                jac_ana_line = np.array(fct_obs_equiv_deriv_ana(*args_for_partial_dev))
                rowbea['A_xbea_x'] = jac_ana_line[0]
                rowbea['A_xbea_y'] = jac_ana_line[1]
                rowbea['A_xbea_z'] = jac_ana_line[2]
                rowbea['A_c']      = jac_ana_line[3]

                
            ######### Symbolic derivation
            #jac_sym = fct_obs_equiv_deriv_symbo(*args_for_partial_dev)
                        
            ####### Update the row
            DFbea.loc[i_rowbea] = rowbea
            
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
        
        
        

    A_bl_stk = []
    B_bl_stk = []
    for inat_rowbl , (i_rowbl , rowbl) in enumerate(DFbl.iterrows()):
        A_bl_line = np.zeros(3*len(ID_bea_list) + 1)
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
    A3 = np.column_stack((A1,A2))
    
    A = np.vstack((A3,*A_bl_stk))
    
    B_twtt = np.hstack(B_twtt_stk)
    
    B = np.hstack((B_twtt,B_bl_stk))
    
    K , Q , P = stats.weight_mat([p_twtt,p_bl],[len(B_twtt),len(B_bl_stk)])
    
    ################# BUILD FINAL DESIGN MATRIX AND VECTOR ###################


    ############################# INVERSION ##################################    
    At   = A.T
    
    N    = At.dot(P).dot(A)
    Ninv = scipy.linalg.inv(N)
    
    dX   = Ninv.dot(At).dot(P).dot(B)
    ############################# INVERSION ##################################    


    ############################# GET NEW VALUES #############################        
    X = np.hstack(Xbea_apri_all)
    X = np.append(X,c_used)
    
    Xnew = X + dX

    Xbea_apri_all_new = Xnew[:-1].reshape(3,3)
    c_ref_new = Xnew[-1]
    
    V = B - A.dot(dX)
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
    DictIteraStore_iter['Xnew']              = Xnew


    ############################ STORE NEW VALUES ############################ 
    
    ####################### COMPUTE NATURAL RESIDUALS ########################
    B_TWTT_Stk = []
    for i_idbea, idbea in enumerate(ID_bea_list):
    
        print("INFO:","Natural residuals for beacon",idbea)
        
        DFbea = DFbig[(DFbig['ID_BEA'] == idbea) & (DFbig['VALID'])]
        
        xbea_apri = Xbea_apri_all_new[idbea-1]
        
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


    print("Xnew     ",i_itera,DictIteraStore[i_itera]["Xnew"])
    print("dX       ",i_itera,DictIteraStore[i_itera]["dX"])
    print("B        ",i_itera,DictIteraStore[i_itera]["B"])
    print("V        ",i_itera,DictIteraStore[i_itera]["V"])
    print("Vnat     ",i_itera,DictIteraStore[i_itera]["Vnatural"])
    print("S B**2   ",i_itera,np.sum(DictIteraStore[i_itera]["B"]**2))
    print("S V**2   ",i_itera,np.sum(DictIteraStore[i_itera]["V"]**2))
    print("S Vnat**2",i_itera,np.sum(DictIteraStore[i_itera]["Vnatural"]**2))
    

np.linalg.norm(Xbea_apri_all_new[0,:] - Xbea_apri_all_new[1,:])
np.linalg.norm(Xbea_apri_all_new[0,:] - Xbea_apri_all_new[2,:])
np.linalg.norm(Xbea_apri_all_new[1,:] - Xbea_apri_all_new[2,:])



for i in range(n_itera):
    print("B   ",i,DictIteraStore[i]["B"])
    print("V   ",i,DictIteraStore[i]["V"])
    print("Vnat",i,DictIteraStore[i]["Vnatural"])
    print("S B**2   ",i,np.sum(DictIteraStore[i]["B"]**2))
    print("S V**2   ",i,np.sum(DictIteraStore[i]["V"]**2))
    print("S Vnat**2",i,np.sum(DictIteraStore[i]["Vnatural"]**2))




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

