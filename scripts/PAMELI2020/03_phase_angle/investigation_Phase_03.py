#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 12:05:36 2020

@author: psakicki

mk 03 => find an optimal labmda and c
"""

#### Import star style
from geodezyx import *                   # Import the GeodeZYX modules
from geodezyx.externlib import *         # Import the external modules
from geodezyx.megalib.megalib import *   # Import the legacy modules names

from pygoat.seafloorpos import ixblue_process_fcts as ibpf

from scipy.spatial.transform import Rotation



######################################## FCT DEF ##############################

def theta_ambiguities_finder(dPhi_in, 
                             lambda_in,
                             D_in,
                             ambig_range=(-4,5),
                             degrees_out=True):
    
    Cos_ab_stk, Theta_ab_stk = [],[]
    
    K_range = np.arange(*ambig_range)
    
    for k in K_range:
        dPhi_ambi = dPhi_in + 2*np.pi*k
        cos_ab   = (dPhi_ambi * lambda_in) / (2*np.pi*D_in)
        theta_ab = np.arccos(cos_ab)
        
        Cos_ab_stk.append(cos_ab)
        Theta_ab_stk.append(theta_ab)
        
    if degrees_out:
        Theta_ab_stk = np.rad2deg(Theta_ab_stk)
    
    return np.array(Cos_ab_stk), np.array(Theta_ab_stk) , K_range

######################################## FCT DEF ##############################

p_obs = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/01_/PAMELI_BREST_vJupyter_6_4transducers_newRF_newDTime/PAMELI_BREST_vJupyter_6_4transducers_newRF_newDTime.csv"

path_obs = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/02_/"
p_obs = path_obs + "PAM_BST_v221_m2507-1109_m5_d02_bea2_ITRF14_RTKLIB/PAM_BST_v221_m2507-1109_m5_d02_bea2_ITRF14_RTKLIB.csv"
p_obs = path_obs + "PAM_BST_v231_m2507-1245_m4_d06_bea3_ITRF14_RTKLIB/PAM_BST_v231_m2507-1245_m4_d06_bea3_ITRF14_RTKLIB.csv"
p_obs = path_obs + "PAM_BST_v211_m2507-1158_m3_d04_bea1_ITRF14_RTKLIB/PAM_BST_v211_m2507-1158_m3_d04_bea1_ITRF14_RTKLIB.csv"

p_obs = path_obs + "PAM_BST_v211_m2507-1158_m3_d04_bea1_RGF93_RTKLIB/PAM_BST_v211_m2507-1158_m3_d04_bea1_RGF93_RTKLIB.csv"

############################## PREPARE DFbig ##################################
DFbig = pd.read_csv(p_obs,index_col=0)
DFbig_orig = DFbig.copy()
DFbig = DFbig.dropna()
DFbig['VALID'] = np.ones(len(DFbig['TWTT'])).astype(bool)
DFbig['TWTTorig'] = DFbig['TWTT'].values
DFbig["Phase_raw"] = DFbig["Phase_raw"]  + np.pi
DFbig.rename(columns={'Phase_raw':'Ph_raw'},inplace=True)

DFgroup = DFbig.groupby("date_emi")
############################## PREPARE DFbig ##################################

################################ PRINT OPTIONS ################################
pd.options.display.float_format = '{:.3f}'.format
pd.set_option("display.max_columns",10)
################################ PRINT OPTIONS ################################

################################ LEVER ARM DEF ################################
lever_arm_RPY_dic = dict()
# !!!!! CHECK LEVER ARMS and DIRECTIONs !!!!!
# -0.039 m / +0.003 m / +1.481 m
lever_arm_RPY_dic["GPS"]  = np.array([-0.039,+0.003,-1.481])
lever_arm_RPY_dic["AHD"]  = np.array([0.,0.,0.3592])

lever_arm_RPY_dic["TDC1"] = np.array([ 0.10735,0.,0.5095])
lever_arm_RPY_dic["TDC2"] = np.array([-0.10735,0.,0.5095])
lever_arm_RPY_dic["TDC3"] = np.array([0.,-0.10735,0.5733])
lever_arm_RPY_dic["TDC4"] = np.array([0., 0.10735,0.5733])
################################ LEVER ARM DEF ################################

################################ APRIORI DICT  ################################
ApriDict = dict()
### Index is the Beacon
# Apri iXBlue
ApriDict[1] = [-7.1674,24.2982,40.0700]
ApriDict[2] = [-7.5010,-26.2188,40.6600]
ApriDict[3] = [22.4663,2.2225,39.6800]

ApriDict[1] = [-7.0959, 24.3475,39.8974]
ApriDict[2] = [-7.3685,-26.1848,40.2508]
################################ APRIORI DICT  ################################

################################ CONSTANT DEF #################################
cref    = 1510
cref    = 1450
cref    = 1500
cref    = 1535

freqref = 26000
ncy     = "ncy_i"

lambdaa = ( cref / freqref )
################################ CONSTANT DEF #################################




def cost_fct(args_in,Ambig_stk=None,short_output=True,verbose=False):
    
    if Ambig_stk:
        estimate_ambig = False
        Ambig_X_stk,Ambig_Y_stk = Ambig_stk
    else:
        estimate_ambig = True
        Ambig_X_stk,Ambig_Y_stk = [],[]
    
    c_in,freq_in = args_in
    
    DiffStk = []
    
    iepoc = -1
    
    
    DFab_stk = []
    
    for epoc , DFepoc_orig in DFgroup:
        
        iepoc+=1
        
        if verbose:
            print('**********************************************************************************************')
        
        #'2019-07-25 09:59:30.729405'
        if epoc != '2019-07-25 10:14:04.450046' and 0:
            continue
        
        if len(DFepoc_orig) != 4:
            continue
        
        DFepoc = DFepoc_orig.copy()
        DFepoc.sort_values(by="TWTT",inplace=True)
    
        DFtwtt = ibpf.generate_DFtwtt(DFepoc)
        
        vNED_ixb = DFepoc[["vN_nrm","vE_nrm","vD_nrm"]].values[0]
        vRPY_ixb = DFepoc[["vR_nrm","vP_nrm","vY_nrm"]].values[0]
        
        NED_src = DFepoc[['N_AHD_rec', 'E_AHD_rec', 'D_AHD_rec']].mean().values
        NED_bea = np.array(ApriDict[DFepoc["ID_BEA"].unique()[0]])
        
        vNED_nat = NED_bea - NED_src
        vNED_nat = vNED_nat / np.linalg.norm(vNED_nat)
        
        ###### Manage rotation
        Euler_all       = DFepoc[['head_rec','pitc_rec','roll_rec']].mean().values  
        
        Euler_head_only = Euler_all.copy()    
        Euler_head_only[1] = 0
        Euler_head_only[2] = 0
    
        Euler_pi_ro_only = DFepoc[['head_rec','pitc_rec','roll_rec']].mean().values         
        Euler_pi_ro_only[0] = 0
        
        Euler_ro_pi_only = DFepoc[['head_rec','roll_rec','pitc_rec']].mean().values          
        Euler_ro_pi_only[0] = 0
        
        pitch = np.rad2deg(Euler_all[1])
        roll  = np.rad2deg(Euler_all[2])
        if verbose:
            print("%%% pitch,roll",pitch,roll)
    
        Rot_all        = Rotation.from_euler("zyx",Euler_all,       degrees=False)
        Rot_head_only  = Rotation.from_euler("zyx",Euler_head_only, degrees=False)
        Rot_pi_ro_only = Rotation.from_euler("zyx",Euler_pi_ro_only,degrees=False)
        Rot_ro_pi_only = Rotation.from_euler("zxy",Euler_ro_pi_only,degrees=False)
        
        Rot_head_only.as_euler('zyx', degrees=True)
        
        vNED_frm_RPY = Rot_head_only.apply(vRPY_ixb)
        vRPY_frm_NED = Rot_head_only.inv().apply(vNED_ixb)
        vRPY_nat     = Rot_head_only.inv().apply(vNED_nat)
    
        vGAP_frm_RPY_ropi = Rot_ro_pi_only.inv().apply(vRPY_ixb)   
        vGAP_frm_RPY_piro = Rot_pi_ro_only.inv().apply(vRPY_ixb)   
    
        vGAP_frm_NED = Rot_all.inv().apply(vNED_ixb)   
    
        if verbose:  
            print("vRPY iXBlue           ",vRPY_ixb)
            print("vRPY from vNED PTSA   ",vRPY_frm_NED) 
            print("vRPY Natural SRC > BEA",vRPY_nat) 
            print("###########################################################")    
            print("vNED iXBlue           ",vNED_ixb)
            print("vNED iXBlue from vRPY ",vNED_frm_RPY) 
            print("vNED Natural SRC > BEA",vNED_nat) 
            print("###########################################################")    
            print("vGAP iXBlue from vNED ",vGAP_frm_NED) 
            print("vGAP iXBlue from vRPY ropi",vGAP_frm_RPY_ropi) 
            print("vGAP iXBlue from vRPY piro",vGAP_frm_RPY_piro) 
            
            print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")    
        
        DFab = pd.DataFrame()
        
        for iTDCpair,TDCpair in enumerate(((1,2),(3,4))):
            
            if TDCpair == (1,2):
                TDC_pos_a = lever_arm_RPY_dic["TDC2"] 
                TDC_pos_b = lever_arm_RPY_dic["TDC1"]
                xy ="X"
            else:
                TDC_pos_a = lever_arm_RPY_dic["TDC3"] 
                TDC_pos_b = lever_arm_RPY_dic["TDC4"]    
                xy = "Y"             
                
            vTDCab  = TDC_pos_b - TDC_pos_a
    
            d_TDCab = np.linalg.norm(vTDCab)
            
            dPhi = DFtwtt.loc[TDCpair[1],'Ph_raw'] - DFtwtt.loc[TDCpair[0],'Ph_raw']
            
            
            lambdaa_in = ( c_in / freq_in )
            
            Cos_ab , Theta_ab , K_ab = theta_ambiguities_finder(dPhi,
                                                                lambdaa_in,
                                                                d_TDCab)
            DFab["k" + xy] = K_ab
            DFab["c" + xy] = Cos_ab
            DFab["a" + xy] = Theta_ab
    
        DFab["epoc"] = epoc
        
    
        cGAP_frm_NED , aGAP_frm_NED = ibpf.ned2dir_cos(vGAP_frm_NED,True)
        cRPY_ixb     , aRPY_ixb     = ibpf.ned2dir_cos(vRPY_ixb,True)
    
        if verbose:
            print("aGAP from vGAP            ",aGAP_frm_NED) 
            print("aRPY iXBlue               ",aRPY_ixb) 
            print("aGAP - aRPY               ",aGAP_frm_NED - aRPY_ixb)
            print(DFab)
        
        DFaX = DFab.set_index("kX").aX
        DFaY = DFab.set_index("kY").aY
        
        if estimate_ambig:
            ambig_X = (DFaX - aGAP_frm_NED[0]).abs().idxmin()
            ambig_Y = (DFaY - aGAP_frm_NED[1]).abs().idxmin()
            DiffStk.append(DFaX[ambig_X] - aGAP_frm_NED[0])
            DiffStk.append(DFaY[ambig_Y] - aGAP_frm_NED[1])
            Ambig_X_stk.append(ambig_X)
            Ambig_Y_stk.append(ambig_Y)
        else:
            ambig_X = Ambig_X_stk[iepoc]
            ambig_Y = Ambig_Y_stk[iepoc]
            DiffStk.append(DFaX[ambig_X] - aGAP_frm_NED[0])
            DiffStk.append(DFaY[ambig_Y] - aGAP_frm_NED[1])            
            
        DFab["aGAPX"] = aGAP_frm_NED[0]
        DFab["aGAPY"] = aGAP_frm_NED[1]

        DFab_stk.append(DFab)


    outsum = np.sum(np.square(DiffStk))
    print("OUT SUM:",      outsum,args_in)
    print("checksum ambig",np.sum(Ambig_X_stk),np.sum(Ambig_Y_stk))

    if short_output:
        return outsum
    else:
        DFab_big = pd.concat(DFab_stk)
        OUT = (np.array(Ambig_X_stk),np.array(Ambig_Y_stk)) , outsum , DFab_big
        return OUT


Nfeval = 1

def callbackF(Xi):
    global Nfeval
    print('Iteration {0:4d}'.format(Nfeval))
    Nfeval += 1
    
args_in1         = (cref,freqref)
Ambig_stk1, out , DFab_big  = cost_fct(args_in1,
                            verbose=False,
                            short_output=False)

args_in2         = (1,freqref)
Ambig_stk2, out , DFab_big  = cost_fct(args_in2,
                            Ambig_stk1,
                            verbose=False,
                            short_output=False)

args_in3 = (1341.66420787, 25999.9746865)
args_in3 = (1249.1092347, 25996.82361177)

Ambig_stk3, out , DFab_big = cost_fct(args_in3,
                            verbose=False,
                            short_output=False)

Ambig_stk_iter = Ambig_stk3
args_in_iter   = args_in3

if 0:
    scipy.optimize.minimize(cost_fct,args_in_iter,
                            (Ambig_stk_iter,True),
                            options={'disp':True},
                            callback=callbackF,
                            method='Powell')
    




                
            
            
            
            
    
            
        
    
    
            
        
            
        


    
    

    
    
        



