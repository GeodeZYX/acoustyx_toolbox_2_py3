#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 16:14:46 2020

@author: psakicki
"""

#### Import star style
from geodezyx import *                   # Import the GeodeZYX modules
from geodezyx.externlib import *         # Import the external modules
from geodezyx.megalib.megalib import *   # Import the legacy modules names

from pygoat.seafloorpos import ixblue_process_fcts as ibpf

from scipy.spatial.transform import Rotation


p_obs = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/01_/PAMELI_BREST_vJupyter_6_4transducers_newRF_newDTime/PAMELI_BREST_vJupyter_6_4transducers_newRF_newDTime.csv"

path_obs = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/02_/"
p_obs = path_obs + "PAM_BST_v221_m2507-1109_m5_d02_bea2_ITRF14_RTKLIB/PAM_BST_v221_m2507-1109_m5_d02_bea2_ITRF14_RTKLIB.csv"
p_obs = path_obs + "PAM_BST_v231_m2507-1245_m4_d06_bea3_ITRF14_RTKLIB/PAM_BST_v231_m2507-1245_m4_d06_bea3_ITRF14_RTKLIB.csv"
p_obs = path_obs + "PAM_BST_v211_m2507-1158_m3_d04_bea1_ITRF14_RTKLIB/PAM_BST_v211_m2507-1158_m3_d04_bea1_ITRF14_RTKLIB.csv"

p_obs = path_obs + "/PAM_BST_v211_m2507-1158_m3_d04_bea1_RGF93_RTKLIB/PAM_BST_v211_m2507-1158_m3_d04_bea1_RGF93_RTKLIB.csv"


ApriDict = dict()

### Index is the Beacon
# Apri iXBlue
ApriDict[1] = [-7.1674,24.2982,40.0700]
ApriDict[2] = [-7.5010,-26.2188,40.6600]
ApriDict[3] = [22.4663,2.2225,39.6800]

ApriDict[1] = [-7.0959, 24.3475,39.8974]
ApriDict[2] = [-7.3685,-26.1848,40.2508]

pout = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/02_/PAM_BREST_AmbigCorr02/"
ambig_corr_suffix = "_AmbigCorr02.csv"

# ### Boxin Beacon 2
# p_gaps_02 = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/01_DATA/07-2019_IXBLUE_BREST/0_RAW_DATA/data-gaps/data - jeudi 25 juillet 2019/02 - Repeater Gaps - 11h09 - box in TP2.txt"
# ### Boxin Beacon 3
# p_gaps_03 = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/01_DATA/07-2019_IXBLUE_BREST/0_RAW_DATA/data-gaps/data - jeudi 25 juillet 2019/06 - Repeater Gaps - 12h45 - box in TP3.txt"
# ### Boxin Beacon 1
# p_gaps_01 = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/01_DATA/07-2019_IXBLUE_BREST/0_RAW_DATA/data-gaps/data - jeudi 25 juillet 2019/04 - Repeater Gaps - 11h58 - box in TP1.txt"
# # p_obs = p_gaps_01

DFbig = pd.read_csv(p_obs,index_col=0)
DFbig_orig = DFbig.copy()
DFbig = DFbig.dropna()
DFbig['VALID'] = np.ones(len(DFbig['TWTT'])).astype(bool)
DFbig['TWTTorig'] = DFbig['TWTT'].values
DFbig["Phase_raw"] = DFbig["Phase_raw"]  + np.pi
DFbig.rename(columns={'Phase_raw':'Ph_raw'},inplace=True)


pd.options.display.float_format = '{:.3f}'.format
pd.set_option("display.max_columns",10)

lever_arm_RPY_dic = dict()
# !!!!! CHECK LEVER ARMS and DIRECTIONs !!!!!
# -0.039 m / +0.003 m / +1.481 m
lever_arm_RPY_dic["GPS"]  = np.array([-0.039,+0.003,-1.481])
lever_arm_RPY_dic["AHD"]  = np.array([0.,0.,0.3592])

lever_arm_RPY_dic["TDC1"] = np.array([ 0.10735,0.,0.5095])
lever_arm_RPY_dic["TDC2"] = np.array([-0.10735,0.,0.5095])
lever_arm_RPY_dic["TDC3"] = np.array([0.,-0.10735,0.5733])
lever_arm_RPY_dic["TDC4"] = np.array([0., 0.10735,0.5733])

DFgroup = DFbig.groupby("date_emi")

cref    = 1510
cref    = 1535
cref    = 1450
cref    = 1500


freqref = 26000
ncy     = "ncy_i"

lambdaa = ( cref / freqref )

NCycle_stack = []

StkDFtwtt = []


def ambiguities_estime(DFepocin,Sin):
        
    if 1:
        Dir_RPY = ibpf.direction_vector_finder_lsq(Sin,
                                                   DFepocin,
                                                   lever_arm_RPY_dic,
                                                   cin=cref,
                                                   twtt_key='TWTT',
                                                   with_estim_c=False)
    else:
        Dir_RPY = ibpf.direction_vector_finder_minimize(Sin,
                                                        DFepocin,
                                                        lever_arm_RPY_dic,
                                                        cin=cref,
                                                        twtt_key='TWTT')
    
    DFtwtt = DFepocin[["TWTT","Ph_raw"]].copy()
        
    DFtwtt["dP"] = np.nan
    DFtwtt.index = DFepocin["ID_TDC"]
           
    DFtwtt_min = DFtwtt[DFtwtt.TWTT == DFtwtt.TWTT.min()]
    idx_ref    = DFtwtt_min.index[0]
    phase_ref = DFtwtt.loc[idx_ref,"Ph_raw"]    
    twtt_ref  = DFtwtt.loc[idx_ref,"TWTT"]    
    
    for idx in DFtwtt.index:
        P = reffram.project_point_on_plan(Dir_RPY,
                                          lever_arm_RPY_dic["TDC" + str(idx_ref)],
                                          lever_arm_RPY_dic["TDC" + str(idx)])
        
        # if idx_ref in (3,4):
        #     print("AAAAAAAA",idx_ref,idx,P[2])
        
        dP = np.linalg.norm(P)
        DFtwtt.loc[idx,"dP"] = dP
        DFtwtt.loc[idx,"ncy_f"] = dP/lambdaa
        
        
        ############### here is the selection of the ambig_search_mode
        ambig_search_mode = 1
        ### mode 1: round
        ### mode 2: more complex, find the closest one to the rough time distances
        
        if ambig_search_mode == 1:
            DFtwtt.loc[idx,"ncy_i"] = np.round(dP/lambdaa)
            
        elif ambig_search_mode == 2:
            ncy_floor = np.floor(dP/lambdaa)
            ncy_ceil  = np.ceil(dP/lambdaa)
            
            phi_floor = DFtwtt["Ph_raw"].loc[idx] + ncy_floor*2*np.pi - phase_ref
            phi_ceil  = DFtwtt["Ph_raw"].loc[idx] + ncy_ceil*2*np.pi  - phase_ref
            
            dP_floor = lambdaa * (phi_floor/(2*np.pi)) 
            dP_ceil  = lambdaa * (phi_ceil/(2*np.pi)) 
            
            dP_raw   = DFtwtt.loc[idx,"dP"] 
            
            if np.abs(dP_raw - dP_floor) < np.abs(dP_raw - dP_ceil):
                DFtwtt.loc[idx,"ncy_i"] = ncy_floor
            else:
                DFtwtt.loc[idx,"ncy_i"] = ncy_ceil                
                
        ############### end of the selection of the ambig_search_mode
        
    DFtwtt["Ph2"]   = DFtwtt["Ph_raw"] + DFtwtt[ncy]*2*np.pi - phase_ref
    DFtwtt["ncy2"]  =  DFtwtt["Ph2"] / (2*np.pi)
    
    DFtwtt["TWTT2"] = twtt_ref + DFtwtt["ncy2"] * (1/freqref) * 10**6
    DFtwtt["dP2"]   = lambdaa * DFtwtt["ncy2"] 
    DFtwtt["G"]  = np.abs(DFtwtt["dP2"] - DFtwtt["dP"]) < lambdaa/2
    
    DFtwtt["G"]  = DFtwtt["G"].apply(int)
    
    ###########################################################################
    ##################### estimate an offset for time recallage the TWTT time
    
    y1out = DFepocin["TWTTorig"].values   
    y2out = DFtwtt["TWTT2"].values
    
    ####### MLE
    
    guess = np.array([1,1])   

    if 1:
        def cost_fct_MLERegression(params):
            offset , sd = params[0], params[1] # inputs are guesses at our parameters
            #yhat = intercept + beta*x # predictions# next, we flip the Bayesian question
            y    = y1out
            yhat = y2out + offset
            # compute PDF of observed values normally distributed around mean (yhat)
            # with a standard deviation of sd
            negLL = -np.sum( scipy.stats.norm.logpdf(y, loc=yhat, scale=sd) )# return negative LL
            return(negLL)
    
        
        results1 = scipy.optimize.minimize(cost_fct_MLERegression, guess,
                       method = 'Nelder-Mead')
        print( "recallage offset MLE" , results1.x[0] )

    ####### LSQ
    if 1:
        def cost_fct_LSQ(offset,y1in,y2in):
            #val = np.sum(np.abs(y1in - y2in + offset))
            val = np.sum((y1in - (y2in + offset))**2)
            return val
        
        results2 = scipy.optimize.minimize(cost_fct_LSQ, guess[0],
                           (y1out,y2out) , 
                           method = 'Powell')
        print( "recallage offset LSQ" , results2.x[0] )
        
    offset = results1.x[0]    
    
    DFtwtt["TWTT2"] = DFtwtt["TWTT2"] + offset
    ##### NB : MLE and LSQ are equivalent
    
    ##################### estimate an offset for time recallage the TWTT time
    ###########################################################################    
    return DFtwtt

#### '2019-07-25 12:57:55.474306'
DFepocStk = []

for epoc , DFepoc_orig in DFgroup:
    
    print('**********************************************************************************************')

    
    #'2019-07-25 09:59:30.729405'
    if epoc != '2019-07-25 10:14:04.450046' and 0:
        continue
    
    if len(DFepoc_orig) != 4:
        continue
    
    DFepoc = DFepoc_orig.copy()
    
    DFepoc.sort_values(by="TWTT",inplace=True)
        
    vNED = DFepoc[["vN_nrm","vE_nrm","vD_nrm"]].values[0]
    vRPY = DFepoc[["vR_nrm","vP_nrm","vY_nrm"]].values[0]
        
    NED_src = DFepoc[['N_AHD_rec', 'E_AHD_rec', 'D_AHD_rec']].mean().values
    NED_bea = np.array(ApriDict[DFepoc["ID_BEA"].unique()[0]])
    
    VC_NED = NED_bea - NED_src
    VC_NED = VC_NED / np.linalg.norm(VC_NED)
    


    Euler_all       = DFepoc[['head_rec','pitc_rec','roll_rec']].mean().values  
    
    Euler_head_only = DFepoc[['head_rec','pitc_rec','roll_rec']].mean().values    
    Euler_head_only[1] = 0
    Euler_head_only[2] = 0

    Euler_pi_ro_only = DFepoc[['head_rec','pitc_rec','roll_rec']].mean().values    
    Euler_pi_ro_only[0] = 0
    
    Euler_ro_pi_only = DFepoc[['head_rec','roll_rec','pitc_rec']].mean().values    
    Euler_ro_pi_only[0] = 0
    
    pitch = np.rad2deg(Euler_all[1])
    roll  = np.rad2deg(Euler_all[2])
    print("%%% pitch,roll",pitch,roll)

    Rot_all        = Rotation.from_euler("zyx",Euler_all,degrees=False)
    Rot_head_only  = Rotation.from_euler("zyx",Euler_head_only,degrees=False)
    Rot_pi_ro_only = Rotation.from_euler("zyx",Euler_pi_ro_only,degrees=False)
    Rot_ro_pi_only = Rotation.from_euler("zyx",Euler_ro_pi_only,degrees=False)

    if False:
        Rot = Rot_all
    else:
        Rot = Rot_head_only
    
    VC_RPY  = Rot.inv().apply(VC_NED)
    VC_RPY2 = Rot.inv().apply(vNED)
    VC_NED2 = Rot.apply(vRPY)

    #VC_RPY = Rot.apply(VC_NED)
    
    vRPYin = VC_RPY
    vRPYin = vRPY



    print("vNED PTSA          ",vNED)
    print("vNED from vRPY PTSA",VC_NED2) 
    print("vNED Vect SRC > BEA",VC_NED) 
    print("###########################################################")
    print("vRPY PTSA          ",vRPY)
    print("vRPY from vNED PTSA",VC_RPY2) 
    print("vRPY Vect SRC > BEA",VC_RPY) 
    
    Sinin = [-0.83846442 , 0.49562207, -0.22657488]
    Sinin = vRPY
    SitAziin = ibpf.ned2site_azim(Sinin)    



    if 0:
        print("==================== DV LSQ ==========================")
        Dir_RPY = ibpf.direction_vector_finder_lsq(Sinin,
                                               DFepoc,
                                               lever_arm_RPY_dic,
                                               cin=cref,
                                               twtt_key='TWTT')
        Dir_NED = Rot.inv().apply(Dir_RPY)

        
    if 0:
        print("==================== DV MINI ==========================")
        for meth in ('Powell',"CG",'Nelder-Mead'):
            print("---------------------------------------------------")
            RES = scipy.optimize.minimize(ibpf.direction_vector_finder_cost_fct,
                                          Sinin,
                                          args=(DFepoc,lever_arm_RPY_dic,cref),
                                          method=meth)
            
            RES_RPY = Rot_pi_ro_only.apply(RES.x)
            RES_RPY2 = Rot_ro_pi_only.apply(RES.x)
            RES_NED = Rot_all.apply(RES.x)

            print(RES)
            
            print(RES_RPY)
            print(RES_RPY2)
            print(RES_NED)
            



    if 0:
        print("==================== SA LSQ ==========================")        
        SitAziout = ibpf.siteazim_finder_lsq(SitAziin[:2],
                                             DFepoc,
                                             lever_arm_RPY_dic,
                                             cin=cref,
                                             twtt_key='TWTT',
                                             with_estim_c=False)
        
                

    if 0:
        print("==================== SA MINI ==========================")        
        for meth in ('Powell',"CG",'Nelder-Mead'):
            print("---------------------------------------------------")
            RES = scipy.optimize.minimize(ibpf.siteazim_finder_cost_fct,
                                          SitAziin[:2],
                                          args=(DFepoc,lever_arm_RPY_dic,cref),
                                          method=meth)
            print(RES)
            
            ibpf.site_azim2ned(list(RES.x) + [1])
            

            
    if 0:
        R = ibpf.direction_vector_finder_cost_fct(RES.x,
                                              DFepoc,
                                              lever_arm_RPY_dic,
                                              cin=cref,
                                              id_tdc_ref = 1,
                                              twtt_key='TWTT',
                                              debug=True,
                                              short_output=True)
        
        
    if 1:
        print("==================== DOA DPHI ==========================")
        DFtwtt = idpf.generate_DFtwtt(DFepoc)
        
        print("************ EPOC",epoc)
        
        Theta_ab_calc_save,Cos_ab_calc_save = [],[]
        
        Cos_ab_dict_save = dict()
        
        for iTDCpair,TDCpair in enumerate(((1,2),(3,4))):
            if TDCpair == (1,2):
                TDC_pos_a = lever_arm_RPY_dic["TDC2"] 
                TDC_pos_b = lever_arm_RPY_dic["TDC1"] 
            else:
                TDC_pos_a = lever_arm_RPY_dic["TDC3"] 
                TDC_pos_b = lever_arm_RPY_dic["TDC4"]                 
            
            vTDCab = TDC_pos_b - TDC_pos_a
            vTDCab = vTDCab

            if DFtwtt.loc[TDCpair[0],'TWTT'] < DFtwtt.loc[TDCpair[1],'TWTT']:
                print("§§§ order 1")
                idtdc_1st,idtdc_2nd = TDCpair[0],TDCpair[1]
            else:
                print("§§§ order 2")
                idtdc_1st,idtdc_2nd = TDCpair[1],TDCpair[0]
                
            idtdc_1st,idtdc_2nd = TDCpair[0],TDCpair[1]        
            
            dPhi = DFtwtt.loc[idtdc_2nd,'Ph_raw'] - DFtwtt.loc[idtdc_1st,'Ph_raw']
            
            if False:
                print("WARN: flip Phi !!!!!")
                dPhi = -dPhi
            
            D = np.linalg.norm(vTDCab)
            
            vRPY_4_DOA = VC_RPY
            theta_ab_true = np.arccos(np.dot(vRPY_4_DOA,vTDCab) / (np.linalg.norm(vRPY_4_DOA) * np.linalg.norm(vTDCab)))
            #print("theta_ab_true",idtdc_1st,idtdc_2nd,theta_ab_true,np.rad2deg(theta_ab_true))
            
            Theta_ab_stk, Cos_ab_stk = [],[]
            K_range = np.arange(-4,5)
            for k in K_range:
                dPhi_ambi = dPhi + 2*np.pi*k
                cos_ab   = (dPhi_ambi * lambdaa) / (2*np.pi*D)
                theta_ab = np.arccos(cos_ab)
                
                Cos_ab_stk.append(cos_ab)
                Theta_ab_stk.append(theta_ab)
                #print("theta_ab     ",k,idtdc_1st,idtdc_2nd,dPhi_ambi,theta_ab,np.rad2deg(theta_ab))
            
            idxmin = np.nanargmin(np.abs(np.array(Theta_ab_stk) - theta_ab_true))
            
            theta_ab_calc = np.rad2deg(Theta_ab_stk[idxmin])
            
            Cos_ab_dict_save[iTDCpair] = Cos_ab_stk
            
            fct_conv = lambda x: x
            fct_conv = np.rad2deg
            
            print("n ambi,θXB,θPS,dθ", TDCpair,
                             K_range[idxmin],
                             fct_conv(theta_ab_true),
                             fct_conv(Theta_ab_stk[idxmin]),
                             fct_conv(theta_ab_true)-fct_conv(Theta_ab_stk[idxmin]))
            
            Theta_ab_calc_save.append(theta_ab_calc)
            Cos_ab_calc_save.append(Cos_ab_stk[idxmin])
            
            
        #for alpha,beta in itertools.product(Cos_ab_dict_save[0],Cos_ab_dict_save[1]):
        alpha,beta = Cos_ab_calc_save
            
        gamma = np.sqrt(1 - alpha**2 - beta**2)
        theta_ab_z = np.arccos(gamma)
        
        vRPY_estim = np.array([alpha,beta,gamma])
        
        vRPY_estim1 = Rot_pi_ro_only.apply(vRPY_estim)
        vNED_estim2 = Rot_all.apply(vRPY_estim)
        vRPY_estim3 = Rot_ro_pi_only.apply(vRPY_estim)
        
        
        print("vRPY estim         ",vRPY_estim,vRPY_estim-vRPY)         
        print("vRPY estim pi/ro   ",vRPY_estim1,vRPY_estim1-vRPY) 
        print("vRPY estim ro/pi   ",vRPY_estim3,vRPY_estim3-vRPY) 
        
        print("vNED estim all     ",vNED_estim2) 

        
  
                
    if 0:
        for iii in range(1):
        
            DFtwtt_step1 = ambiguities_estime(DFepoc,Sin = vRPYin)
            print(str(iii) * 80)
            print(DFtwtt_step1)  
            print(str(iii) * 80)
        
            DFepoc["TWTT"] = DFtwtt_step1["TWTT2"].values
            
            if DFtwtt_step1["G"].sum() == 4:
                DFepoc["VALID"]  = True
                break
            else:
                DFepoc["VALID"]  = False

        DFepocStk.append(DFepoc)        
        StkDFtwtt.append(DFtwtt_step1)
    

DFepocNew1 = pd.concat(DFepocStk)
DFepocNew  = DFepocNew1[DFepocNew1["VALID"]]

DF3       = DFepocNew.groupby("date_emi")



solve_rate = DFepocNew1["VALID"].sum() / len(DFepocNew1["VALID"])
print("solve_rate",solve_rate)

# for epoc,DFepoc_c_empi in DF3:
#     DFF1 = ibpf.c_empiric(DFepoc_c_empi,lever_arm_RPY_dic,twtt_key="TWTTorig")
#     DFF2 = ibpf.c_empiric(DFepoc_c_empi,lever_arm_RPY_dic,twtt_key="TWTT")
#     print(pd.concat((DFF1,DFF2),axis=1))

utils.create_dir(pout)

DFepocNew.to_csv(pout + os.path.basename(p_obs)[:-4] + ambig_corr_suffix)

DFtwttBig = pd.concat(StkDFtwtt)
DFtwttBig["G"].sum()




        
NCycle = (DFtwttBig["ncy_f"] - np.floor(DFtwttBig["ncy_f"]))
    
# NCycle = np.hstack(NCycle_stack)
NCycle = NCycle[NCycle != 0]
plt.hist(NCycle,bins=100)
    

    
    # dphi_12 = phi_2 - phi_1
    # dphi_34 = phi_4 - phi_3
    
    # d12 = 0.10735*2
    # d34 = 0.10735*2
    
    # n = 1
    
    # dphi_12_n = (dphi_12 + n*2*np.pi)
    # dphi_34_n = (dphi_34 + n*2*np.pi)
    
    # theta_x = np.arccos((lambdaa * dphi_12_n) / (2 * np.pi * d12))
    # theta_y = np.arccos((lambdaa * dphi_34_n) / (2 * np.pi * d34))
    
    # R = 40
    
    # xu = (R * lambdaa * dphi_12_n)/(2 * np.pi * d12)
    # yu = (R * lambdaa * dphi_12_n)/(2 * np.pi * d34)
    
    # llll = np.linalg.norm([xu,yu])
    
#   print(theta_x,theta_y,xu,yu,llll)











