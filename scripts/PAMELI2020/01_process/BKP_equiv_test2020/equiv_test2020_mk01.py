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


p_obs = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/01_/PAMELI_BREST_vJupyter_2_4transducers/PAMELI_BREST_vJupyter_2_4transducers.csv"

DFbig = pd.read_csv(p_obs,index_col=0)
DFbig_orig = DFbig.copy()
DFbig = DFbig.dropna()
DFbig['VALID']   = np.ones(len(DFbig['TWTT'])).astype(bool)


############################ LOADING DATA ####################################


############################ SET ENVIRONEMENT VARIABLE #######################
cols_posi = ['N_AHD_emi', 'E_AHD_emi', 'D_AHD_emi',
             'N_AHD_rec', 'E_AHD_rec', 'D_AHD_rec']

cols_posi_emi = cols_posi[0:3]
cols_posi_rec = cols_posi[3:6]

ID_bea_list = np.sort(DFbig['ID_BEA'].unique())

n_itera = 15

DictIteraStore = dict()

p_twtt = 10**-6
p_c    = 1


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

    DFbig['TWTT_mod'] = np.zeros(len(DFbig['TWTT']))
    DFbig['B_TWTT']   = np.zeros(len(DFbig['TWTT']))
    
    
    for i_idbea, idbea in enumerate(ID_bea_list):
    
        print("INFO:","partial derivative for beacon",idbea)

        
        DFbea_orig = DFbig[DFbig['ID_BEA'] == idbea]
        DFbea = DFbea_orig.copy()
        
        xbea_apri = Xbea_apri_all[idbea-1]
        
        Xemi     = DFbea[cols_posi_emi].values
        Xrec     = DFbea[cols_posi_rec].values
        TWTT_obs = DFbea['TWTT'].values * 10**-6
        
         
        TWTT_mod_stk = []
        RESstk = []   
        
        c_bea = equiv_straght_ray(Zinp, Cinp, xbea_apri[2])
        
        c_used = c_ref
         
        Jac_xbea_stk = []
        Jac_c_stk    = []

        ZIP = list(zip(Xemi,Xrec,TWTT_obs))    
        for xemi,xrec,twtt in ZIP:
    
            ###################### MODELING & DERIVATION #######################
            
            args_for_partial_dev =  (xbea_apri[0],
                                     xbea_apri[1],
                                     xbea_apri[2],
                                     xemi,
                                     xrec,
                                     c_used)
            
            
            ########## Modeled observations
            twtt_mod = fct_obs_equiv(*args_for_partial_dev)
            TWTT_mod_stk.append(twtt_mod)
            
            ########## derivative of position
            jac_xbea_line = []
            for i in range(3):
                jac_xbea = stats.partial_derive(fct_obs_equiv, i,
                                                args_f = args_for_partial_dev)
                jac_xbea_line.append(jac_xbea)
    
            if with_numerical_diff:
                Jac_xbea_stk.append(jac_xbea_line)
                        
            ########## derivative of c
            jac_c = stats.partial_derive(fct_obs_equiv, 'c_bea_in',
                                         args_f = args_for_partial_dev)    

            if with_numerical_diff:
                Jac_c_stk.append(jac_c)

            ######### Analytic derivation
            jac_ana = np.array(fct_obs_equiv_deriv_ana(*args_for_partial_dev))

            ######### Symbolic derivation
            #jac_sym = fct_obs_equiv_deriv_symbo(*args_for_partial_dev)

            if not with_numerical_diff:
                Jac_xbea_stk.append(jac_ana[:3])
                Jac_c_stk.append(jac_ana[-1])
            
            ###################### MODELING & DERIVATION #######################
    
        ### Stacking the values in Matrices/Vectors
        TWTT_mod     = np.array(TWTT_mod_stk)
        B_TWTT       = (TWTT_obs - TWTT_mod) 
        Jac_xbea_mat = np.vstack(Jac_xbea_stk)
        
        ### Filling the Big DataFrame with modeled values and derivatives
        Bool_good_rows = ((DFbig['ID_BEA'] == idbea))
        
        DFbig.loc[Bool_good_rows,'TWTT_mod'] = TWTT_mod
        DFbig.loc[Bool_good_rows,'B_TWTT']   = B_TWTT
        DFbig.loc[Bool_good_rows,'A_xbea_x'] = Jac_xbea_mat[:,0]
        DFbig.loc[Bool_good_rows,'A_xbea_y'] = Jac_xbea_mat[:,1]
        DFbig.loc[Bool_good_rows,'A_xbea_z'] = Jac_xbea_mat[:,2]
        DFbig.loc[Bool_good_rows,'A_c']      = Jac_c_stk
        
    A_c_stk    = []
    A_xbea_stk = []
    B_stk      = []
    
    for i_idbea, idbea in enumerate(ID_bea_list):
        
        Bool_idbea = (DFbig['ID_BEA'] == idbea)
        
        B_stk.append(DFbig.loc[Bool_idbea,'B_TWTT'].values)
        
        A_c_stk.append(DFbig.loc[Bool_idbea,'A_c'].values)
        A_xbea_stk.append(DFbig.loc[Bool_idbea,['A_xbea_x',
                                                'A_xbea_y',
                                                'A_xbea_z']].values)    
    
    ################# BUILD FINAL DESIGN MATRIX AND VECTOR ###################
    A1 = scipy.linalg.block_diag(*A_xbea_stk)
    A2 = np.hstack(A_c_stk)
    A = np.column_stack((A1,A2))
    
    B = np.hstack(B_stk)
    
    K , Q , P = stats.weight_mat([p_twtt],[len(B)])
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
    
    V = B - A.dot(Xnew)
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
    ############################ STORE NEW VALUES ############################ 
    
    ####################### COMPUTE NATURAL RESIDUALS ########################
    B_TWTT_Stk = []
    for i_idbea, idbea in enumerate(ID_bea_list):
    
        print("INFO:","Natural residuals for beacon",idbea)
        
        DFbea = DFbig[DFbig['ID_BEA'] == idbea]
        
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
        for xemi,xrec,twtt in ZIP:
        
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

for i in range(n_itera):
    print("B   ",i,DictIteraStore[i]["B"])
    print("V   ",i,DictIteraStore[i]["V"])
    print("Vnat",i,DictIteraStore[i]["Vnatural"])






# plt.figure()

# plt.plot(V,label="V")
# plt.plot(Vnatural,label="Vnat")
# plt.legend()

# Diff1 = V - Vnatural
# Diff2 = V - np.flip(Vnatural)

# plt.figure()

# plt.plot(Diff1)
# plt.plot(Diff2)




