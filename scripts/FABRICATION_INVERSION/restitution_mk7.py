#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May  4 15:30:07 2017

@author: psakicki
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 13:04:35 2015

@author: psakicki

v4 : en mode mono Z , on injecte le dZi comme observable
MK6v2 : implementation du sigma_batch_mode (160616)
MK6v3 : fusion avec GEODESEA
MK7   : introduction des fichiers de config
"""
# =========== Scientific librairies
import matplotlib
#matplotlib.use('Agg') # Must be before imported matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy import interpolate
import scipy.optimize as optimize
import scipy.stats
import scipy.sparse as sprs
import sympy
import pandas as pd

# =========== General librairies
import sys
from tempfile import mkdtemp
import time
import datetime as dt
import dateutil
import os
import inspect
import itertools
import multiprocessing as mp
import glob
import platform
import shutil
import subprocess
import dateutil.parser
import collections
import re
import copy
import tabulate
from configparser import SafeConfigParser
from configparser import ConfigParser

# =========== Specifics (Geodesy/Acoustics) librairies

#### Import star style
from geodezyx import *                   # Import the GeodeZYX modules
from geodezyx.externlib import *         # Import the external modules
from geodezyx.megalib.megalib import *   # Import the legacy modules names


import acouclass as acls
import raytrace as rt
import SSP as ssp




##### ================= set the path of the config file here =================
configfile_path = '/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/config_files/GEODESEA.cfg'
configfile_path = '/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/config_files/CANOPUS.cfg'
configfile_path = "/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/config_files/Detection15mai_CANOPUS_1505.cfg"
configfile_path = '/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/config_files/FINALc_Simulation_exemple.cfg'
configfile_path = '/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/config_files/Simulation_1709A_WithorWithoutLateralGradient.cfg'
configfile_path = '/home/adminuser/Documents/CODES/acoustyx_toolbox_2/config_files/Detection15mai_CANOPUS_1505_RealCampaign_exemple.cfg'

configfile_path = '/home/psakicki/CODES/acoustyx_toolbox_2_py3/config_files/Detection15mai_CANOPUS_1505_RealCampaign_exemple.cfg'
configfile_path = '/home/psakicki/CODES/acoustyx_toolbox_2_py3/config_files/PAMELI_BREST_wrong_lever_arm_5cm.cfg'
configfile_path = '/home/psakicki/CODES/acoustyx_toolbox_2_py3/config_files/PAMELI_BREST_wrong_lever_arm.cfg'
configfile_path = '/home/psakicki/CODES/acoustyx_toolbox_2_py3/config_files/PAMELI_BREST.cfg'

###############################################################################


############# SETING GENERIC VALUES to avoid editor errors
multiexp_mode , batch_mode    = 0 , 0
purge_prelim_mode , purge_mode = 0 , 0
gene_path , expprefix , exp , expdirsuffix = '' , '' , '' , ''
alternat_ssp_path = ''
smart_SSP_Sdic_path = ''
force_mode_reading = 0
mnop_files_ext = ''
with_zmaster = 0
sigma_batch_mode = 0
plot_mode , force_mode = 0 , 0
with_smart_ssp = 0
jk_rand_seed = 42
sigma_defo_zmaster = 42
sigma_defo_asm = 42
sigma_defo_bl  = 42
nbproc = 0
iitermax = 0
cleaning_coef = 0
dX_sum_stop  = 0
time_window_mode = ''
with_time_window = 0
start,end = [] , []
start_end_multi = []
exp_custom_name = ''

with_barycenter            = 0 # Use this mode for a direct determination of the beacons array barycenter coordinate and associated formal sigmas
with_bl                    = 0 # If avaiable, use distance information between the seafloor beacon
                               # data are stocked in a .B.dat file
                               
with_cleaning              = 1 # Activate cleaning with resdual standard dev : keeping the ping if  v_ping < cleaning_coef * std(V_ping) (check the value of the cleaning_coef variable ! usually 3s)
with_v_4_p_reinject        = 1 # use residuals of the previous iteration as weights for the new one

with_mono_z                = 0 # Estimate the same Z (depth) for all the beacons
with_zmaster               = 0 # Use depth differences between beacons, and restim them (experimental)
with_dzcst                 = 0 # Use depth differences between beacons, as constant (experimental)
                               # THOSE 3 OPTIONS ARE ONLY FOR SIMULATIONS
                               # the ref beacon is the first one
                               # Depth informations are in the .O/P.dat file (containing the pings)
                               # and is the apriori depth of each beacon (line 'pxp_coords')

# Some modes are impossible :
#     with_mono_z and with_zmaster 
#     with_zmaster and not with_bl
#     with_dzcst and with_zmaster

with_noising_apriori       = 0 # This part shall be activated only for simulation, it activates or not the [sigmas_noise_apriori] part, in order to perturbate the apriori

with_time_window           = 0 # apply a timewindow, using start/end or start_end_multi (see below)
with_specific_timwin_bool  = 0 # discontinued

with_jackknife             = 0 # Apply a Jack Knife, i.e. picking randomly some pings in a full set. The behavior of the Jack Knife can be changed by the variable jk_rand_seed (an int which controls the random selection), and the keep_ratio (ranging between  0 & 1)
with_invert_jk             = 0 # Invert the selection done by the Jack Knife

with_fb_mode               = 1 # Activate the forward/backward mode : using a different raytracing for the surface=>seafloor & the seafloor=>surface shoot (this mode is more accurate so shall be ON) 
with_old_style             = 0 # If True, select the F or the B as the time of ping, If false time of ping = TWTT/2 (only for debug), WORKS ONLY IF with_fb_mode == False
with_backward_select       = 0 # use the backward raytracing (only for debug), WORKS ONLY IF with_fb_mode == False and_with_old_style == True !!!!
# In fact with_fb_mode take as input the Two way travel time (TWTT/2), and not the TWTT/2
# with_fb_mode is False, then it is the former processing way using TWTT/2, and with_old_style allows to get TWTT/2 in data with specific keywords
# Nevertheless, with_fb_mode = True shall be used now

with_ssp_bilin             = 0 # Use a bilinear Sound Speed (experimental & unstable)
with_decimate_ssp          = 1 # Decimate the SSP using the IAPSO 1936 convention (enhance calc speed, so recommended)
with_alternat_ssp          = 0 # Using an alternative SSP, defined by the alternat_ssp_path
with_smart_ssp             = 0 # use SSPs in a dict (specified by) (very specific case)
with_munk_ssp              = 1 # use the munk generic sound speed profile instead of the real one
with_print_simu_infos      = 0

# ===========================================================================
# ===================      CONFIG FILE READING      =========================
# ===========================================================================

parser = ConfigParser()
parser.read(configfile_path)    

ConfigDic = parser._sections

########### BRUTFORCE READING  ###########
# Brutforce reading of the config file Dic, in order to load all the parameters
print('============================================')
print("INFO : READING THE CONFIG FILE")
print("")
for ksec , Dicsec in list(ConfigDic.items()):
    for kparam , valparam in  list(Dicsec.items()):
        if kparam == '__name__':
            continue
        else:
            print(kparam + ' = ' + valparam)
            exec(kparam + ' = ' + valparam)
print("")
print("INFO : END OF  CONFIG FILE READING")
print('============================================')
print("")


########### PROCESSING OPTIONS ###########
# reading of the config file Dic, processing_options section,
# in order to generate the booldic_lis
# for differents processing_options

no_variable_dict = dict() 
variable_dict    = dict() 

for kopts , vopts in list(ConfigDic['processing_options'].items()):
    if kopts == '__name__':
        continue
    else:
        vopts = eval(vopts)
        if genefun.is_iterable(vopts):
            variable_dict[kopts] = vopts
        else:
            no_variable_dict[kopts] = vopts

bool_kw_lis = list(variable_dict.keys())
                
booldic_lis = geok.kwargs_for_jacobian(no_variable_dict,
                                       variable_dict)

########### SIGMA READING ###########
# reading of the config file Dic, sigmas section,
# in order to generate the booldic_lis
# for differents sigmas

no_variable_sigmas_dict = dict() 
variable_sigmas_dict    = dict() 

for sigmas in ('sigmas_noise_apriori','sigmas_weights'):
    for ksigmas , vsigmas in list(ConfigDic[sigmas].items()):
        if ksigmas == '__name__':
            continue
        else:
            vsigmas = eval(vsigmas)
            if genefun.is_iterable(vsigmas):
                variable_sigmas_dict[ksigmas]    = vsigmas
            else:
                no_variable_sigmas_dict[ksigmas] = vsigmas
                
sigmadic_lis = geok.kwargs_for_jacobian(no_variable_sigmas_dict,
                                       variable_sigmas_dict)

# ===========================================================================
# =================== END OF CONFIG FILE READING    =========================
# ===========================================================================


############# LOADING THE EXPERIENCES
if multiexp_mode == 'single':
    exp_lis = [exp]
elif multiexp_mode == 'list':
    exp_lis = exp_lis
elif multiexp_mode == 'prefix':
    exp_lis   = glob.glob(os.path.join(gene_path,expprefix + '*'))
    exp_lis   = list(reversed(sorted([os.path.basename(e) for e in exp_lis])))


############# CLEANING ALL THE FOLDERS IN ONE COMMAND
if purge_prelim_mode:
    print("WARN : purge mode !!! (preliminary), all logfiles will be removed !")
    print("INFO : you have 10sec to cancel")
    print("INFO : Ctrl+C or close the console to cancel")
    time.sleep(12)
    
    for expp in exp_lis:
        for ext in ("/*.log","/*.exp","/*.png","/*.V",
                    "/*.smartV","/*.pdf"):
            rmlis = glob.glob(gene_path + '/' + expp + ext)
            if len(rmlis) == 0:
                print('purge prelim. : remove list empty !!!')
            for f in rmlis:
                print('INFO : remove : ' , f)
                os.remove(f)       
                
############# TIME WINDOW MANAGEMENT

if not with_time_window:
    startend_lis = [[[1970,1,1,1,1,1],[2038,1,1,1,1,1]]]
else:
    if time_window_mode == 'single':
        startend_lis = [[start,end]]
    else:
        startend_lis = start_end_multi
    
        
iterlist = itertools.product(exp_lis,booldic_lis,sigmadic_lis,startend_lis)
              
# no interactive plot
plt.ioff()

for exp , booldic , sigmadic , startend in  iterlist:
    
    print("##### EXPERIENCE : " , exp_lis.index(exp) + 1 , '/' , len(exp_lis)) 
    
    expini   = exp
    exp_path = os.path.join(gene_path,expini + expdirsuffix)
    exp_name = os.path.basename(expini)
    
    ############# READING THE DATA
    try:
        bigdico  = acls.give_me_the_path(exp_path,exp_name) #,[1,3,5])
    except Exception as e:
        if force_mode_reading:
            print(e)
            print("WARN : error while reading but force_mode_reading is ON")
            continue
        else:
            raise e
    
#    if 'OPERA' in exp_path and 'GEODESEA' in exp_path:
#        MNOPfile = 'O'
#    elif 0:
#        MNOPfile = 'M'  
#        MNOPfile = 'P'
#    else:
    
    ############# CHOOSING THE EXTENSION
    if mnop_files_ext  == 'auto':
        MNOPfile = acls.auto_find_extension_MNOP_files(exp_path)
    elif mnop_files_ext not in ('M', 'N', 'O' , 'P'):
        print('ERR : file extension is not auto nor M, N, O, P')
        MNOPfile = None
    else:
        MNOPfile = mnop_files_ext

    ############# CLEANING THE FOLDERS ONE BY ONE
    if purge_mode and not purge_prelim_mode:
        print("WARN : purge mode !!! all logfiles will be removed !")
        print("INFO : you have 5sec to cancel")
        print("INFO : Ctrl+C or close the console to cancel")
        time.sleep(8)
        
        for ext in ("/*.log","/*.exp","/*.png","/*.V"):
            rmlis = glob.glob(exp_path + ext)
            for f in rmlis:
                print('INFO : remove : ' , f)
                os.remove(f)
    
    try:
        if batch_mode:
            genefun.eval_a_dict(booldic,globals(),verbose=False)
            if with_zmaster:
                with_barycenter = 1
            bool_4_exp = [str(int(eval(ee))) for ee in bool_kw_lis]
            bool_4_exp_str = ''.join(bool_4_exp)
        else:
            bool_4_exp_str = ''
            
        if 'V_4_P_reinject' in list(globals().keys()):
            del V_4_P_reinject
        
        #np.set_printoptions(threshold=np.nan)
        
        
        timestamp_start = gf.get_timestamp()
        
        exp_path_output = os.path.join(exp_path,
                                       exp_name + '_' + bool_4_exp_str + '_' + timestamp_start + '_' + exp_custom_name)
        gf.create_dir(exp_path_output)
        
        # Start to write the  logfile
        F = genefun.Tee_frontend(exp_path_output , exp_name , 
                                 bool_4_exp_str + '_' + timestamp_start + '_' + exp_custom_name , 
                                 print_timestamp=False )
        # here is the path basis for all the files which will be writen
        protopath = os.path.join( exp_path_output , 
                                 exp_name + '_' + bool_4_exp_str + '_' + timestamp_start + '_' + exp_custom_name + '_')

  
        # ==============================================
        # PARAMETERS INITIALISATION
        # ==============================================
        
        nPXP_proto = len(bigdico[MNOPfile])
        idPXP_lis = list(bigdico[MNOPfile].keys())
        
        
        # PARAMETRE BRUIT PXP / Z
        #sigma_pxp_apri = 10
        #kmeta_pxp_apri = 10 # le K normal c'est l'ID du PXP, le meta c'est pour tous
        #                   # le Kmeta sera un multiple de 10 par convention
        #
        #sigma_z_apri  = 1 # Ne sert qu'en z mono, pour l'apriori du zmono ...
        #k_z_apri      = 42 + 10
        #
        #k_dz_apri     = 1789 + 10
        #sigma_dz_apri = 10**-6
        #     
        #k_bl_apri     = 852369 + 10
        #sigma_bl_apri = 10**-2
        #
        ## PARAMETRES ANNEXES
        #kriter    = 10**7
        #kriterold = 10**10
        #iiter = 3
        #h = 0 
        ## for jackknife
        #keep_ratio   = 0.5
        #keep_ratio   = 0.1
        #jk_rand_seed = 1452 Le jk_seed a été ramméné avant l'interprétation du dico
        
        if sigma_batch_mode:

            if not with_noising_apriori:
                k_z_apri      = 0
                k_dz_apri     = 0
                k_bl_apri     = 0            
            # Construction d'un Mersene Twister propre a chaque sigma_apri
            else:
                k_z_apri      *= np.abs(np.round(np.log10(sigma_z_apri))  + 11)
                k_dz_apri     *= np.abs(np.round(np.log10(sigma_dz_apri)) + 22)
                k_bl_apri     *= np.abs(np.round(np.log10(sigma_bl_apri)) + 33)
                k_z_apri      = int(k_z_apri)
                k_dz_apri     = int(k_dz_apri)
                k_bl_apri     = int(k_bl_apri)


        PXPold = np.array([9999999,9999999,9999999] * nPXP_proto) 
        PXPnew_stk = []
        ObsASM_lis = []
        Xbato_lis  = []
        TTT_lis    = []

        PXP_lis     = []
        PXPapri_lis = []
        
        PXPapri0_lis = []

        fuv = 1

        expdic = collections.OrderedDict()
        iiter = 3
        iterdic = acls.get_iterdic_in_expdic(expdic,iiter)
        
        A , N , Ninv , P , B = [],[],[],[],[]
    
        # PARAMETRES IMPORTANTS

        #if not batch_mode:
        #    if MNOPfile != 'P': 
        #        with_fb_mode         = 0 # True
        #        with_old_style       = 1 # select the F or the B as the time of ping, else TWTT/2               
        #        with_backward_select = 0 # WORKS ONLY IF with_fb_mode == False and_with_old_style == True !!!!
        #    else:
        #        with_fb_mode         = 1 # True
        #        with_old_style       = 0 # select the F or the B as the time of ping, else TWTT/2               
        #        with_backward_select = 0 # WORKS ONLY IF with_fb_mode == False and_with_old_style == True !!!!
        #    

        if with_zmaster:
            with_barycenter = 1

        if with_bl:
            BL = bigdico['B']['d']
        
        if with_ssp_bilin:
            Z = bigdico['2lin_Z']['d']
            C = bigdico['2lin_C']['d'] 
        else:
            Z = bigdico['Z']['d']
            C = bigdico['C']['d'] 
                                
        if with_alternat_ssp:
            alternat_SSP = np.loadtxt(alternat_ssp_path)
            Z = alternat_SSP[:,0]
            C = alternat_SSP[:,1]
    
        if with_munk_ssp:
            print('INFO : with_munk_ssp ON, you will use a Munk profile')
            Z,C = ssp.munk(6000,10)
            
            #C = np.ones(len(Z)) * 1550. + np.random.randn(len(Z)) * 10.
            #raise Exception
            
            # some  non significative random behavior to avoid error 
        RS = np.random.RandomState(seed=42)
        C = C + RS.randn(len(C)) * 10**-4  *0.
        print("C mean",np.mean(C))
        
        Z2, C2 = ssp.SSP_extrapolate(Z, C, 80, 1)

        
        
            
        if with_decimate_ssp:
            print("decimating SSP")
            Zfull , Cfull = Z , C
            Z , C = ssp.SSP_light(Z,C)
            
        if not with_noising_apriori:
            print('INFO : with_noising_apriori is OFF, noise on apriori will be set at 0')
            # in practice, set @ 10**-15 because of a log used after
            sigma_pxp_apri = 10**-12
            sigma_dz_apri  = 10**-12
            sigma_bl_apri  = 10**-12
            sigma_z_apri   = 10**-12
        
        if with_mono_z and with_zmaster:
            print("ERR : with_mono_z and with_zmaster together, impossible combination, skipping ...")
            raise Exception

        if with_zmaster and not with_bl:
            print("ERR : with_zmaster and not with_bl, impossible combination, skipping ...")
            raise Exception

        if with_dzcst and with_zmaster:
            print("ERR : with_dzcst and with_zmaster, impossible combination, skipping ...")
            raise Exception

        if not with_jackknife and with_invert_jk:
            print("ERR : not with_jackknife and with_invert_jk, useless combination, skipping ...")
            raise Exception

        if not with_time_window and with_specific_timwin_bool:
            print("ERR : not with_time_window and with_specific_timwin_bool, useless combination, skipping ...")
            raise Exception

        R_z = np.random.RandomState(k_z_apri)
        err_z = R_z.randn(1) * sigma_z_apri
        i_pxp_master = 0

        for ipxp,Mpxp in sorted(bigdico[MNOPfile].items()):
            
            Mdata = Mpxp['d']
            Mcomm = Mpxp['c']
            
            if MNOPfile == 'P':
                Mtab    = Mpxp['t']  
                
            if MNOPfile == 'N':
                Xbato = list(Mdata[:,1:4])
            elif MNOPfile == 'M':
                Xbato = list(Mdata[:,:3])
            elif MNOPfile == 'O':
                Xbato_f = list(Mdata[:,1:4]  )
                Xbato_b = list(Mdata[:,12:15])
                Xbato   = list(zip(Xbato_f , Xbato_b))
            elif MNOPfile == 'P':
                if with_fb_mode:
                    Xbato_f = list(np.column_stack((Mtab['X_emi_noise'],Mtab['Y_emi_noise'],Mtab['Z_emi_noise'])))
                    Xbato_b = list(np.column_stack((Mtab['X_rec_noise'],Mtab['Y_rec_noise'],Mtab['Z_rec_noise'])))
                    Xbato   = list(zip(Xbato_f , Xbato_b))
                else:
                    if with_backward_select:
                        Xbato = list(np.column_stack((Mtab['X_rec_noise'],Mtab['Y_rec_noise'],Mtab['Z_rec_noise'])))
                    else:
                        Xbato = list(np.column_stack((Mtab['X_emi_noise'],Mtab['Y_emi_noise'],Mtab['Z_emi_noise'])))

            Xbato_lis.append(Xbato)

            if MNOPfile == 'N':
                ObsASM_load = Mdata[:,-1]
            elif MNOPfile == 'M':
                ObsASM_load = Mdata[:,3]
            elif MNOPfile == 'O':
                ObsASM_load = Mdata[:,-1]
            elif MNOPfile == 'P':
                if with_fb_mode:
                    ObsASM_load = np.array(Mtab['TWTT_noise'])
                elif with_old_style:
                    if with_backward_select:
                        ObsASM_load = np.array(Mtab['T_rec_noise'])
                    else:
                        ObsASM_load = np.array(Mtab['T_emi_noise'])
                else:
                    ObsASM_load = np.array(Mtab['TWTT_noise'] * .5)

            ObsASM_lis.append(ObsASM_load)  
            if MNOPfile == 'N':
                TTT_lis.append(Mdata[:,0])
            elif MNOPfile == 'O':
                TTT_lis.append((Mdata[:,0] , Mdata[:,11]))
            elif MNOPfile == 'P':
                if with_fb_mode:
                    TTT_lis.append((Mtab['E_emi_noise'] , Mtab['E_rec_noise']))
                else:
                    if with_backward_select:
                        TTT_lis.append(Mtab['E_rec_noise'])
                    else:
                        TTT_lis.append(Mtab['E_emi_noise'])

            PXP = acls.pxp_string_2_array(Mcomm['pxp_coords'])
            PXP_lis.append(PXP)
            nPXP = len(PXP_lis)

            k_pxp_apri = ipxp + kmeta_pxp_apri
            R_pxp_apri = np.random.RandomState(k_pxp_apri)
            if with_mono_z:
                PXPapri_mono = PXP + np.concatenate((R_pxp_apri.randn(3)[0:2] * sigma_pxp_apri , err_z))
            else:
                if MNOPfile in ('N','O'): 
                    PXPapri_mono = PXP
                elif MNOPfile == 'M':
                    PXPapri_mono = PXP + R_pxp_apri.randn(3) * sigma_pxp_apri
                elif MNOPfile == 'P':
                    print('This par must be splited ...')
                    PXPapri_mono = PXP + R_pxp_apri.randn(3) * sigma_pxp_apri

            PXPapri_lis.append(PXPapri_mono)
            
            PXPapri0_lis.append(PXPapri_mono)

            
        if with_mono_z:
            for i,pxp in enumerate(PXPapri_lis):
                if i_pxp_master != i:
                    pxp[2] = PXPapri_lis[i_pxp_master][2] 
            
        PXPZ_lis       = np.array(PXP_lis)[:,-1]
        PXPdZ_lis_true = PXPZ_lis - PXPZ_lis[i_pxp_master] 
        
        noise4PXPdZ = np.random.RandomState(k_dz_apri).randn(nPXP) * sigma_dz_apri
        PXPdZ_lis   = noise4PXPdZ +  PXPdZ_lis_true  # <= par defaut c'est le PXP1 qui est maitre
        PXPdZ_lis[i_pxp_master] = 0
        PXPtrue_arr = np.array(PXP_lis)
        shape_arr   = PXPtrue_arr.shape
        PXPapri0_arr = np.array(PXPapri_lis)
        PXPapri0_vct = genefun.vectorialize(PXPapri0_arr)
        PXPapri = PXPapri0_vct
        
        if with_bl:
             noise4BLs = np.random.RandomState(k_bl_apri).randn(nPXP , nPXP) * sigma_bl_apri
             noise4BLs = np.zeros((nPXP,nPXP)) + np.triu(noise4BLs,1) + np.triu(noise4BLs,1).T
             BLnoised  = BL + noise4BLs
                                
        if with_barycenter:
            PXPbary0      = acls.barycenter_calc(PXPapri0_arr)
            PXPbary       = PXPbary0
            dPXPapri0_arr = PXPapri0_arr - PXPbary0
            dPXPapri0_vct = genefun.vectorialize(dPXPapri0_arr)
            dPXPapri      = dPXPapri0_vct
            dPXPapri0     = dPXPapri0_vct
            dPXPapri_lis  = list(dPXPapri0_arr)   
        
        if with_zmaster:
            # dPXP p/r au BARYCENTRE
            PXPref0       = np.array(PXPbary0)
            #PXPref0[2]    = PXPapri_lis[i_pxp_master][2] Cette ligne plantait tout (160610)
            PXPref        = PXPref0
            dPXPapri0_arr = PXPapri0_arr - PXPref0
            dPXPapri0_vct = genefun.vectorialize(dPXPapri0_arr)
            dPXPapri0     = dPXPapri0_vct
            dPXPapri      = dPXPapri0_vct
            dPXPapri_lis  = list(dPXPapri0_arr)  
     
        if with_time_window:
            ObsASM_lis_orig = list(ObsASM_lis)
            Xbato_lis_orig  = list(Xbato_lis)
            ObsASM_lis  = []
            Xbato_lis   = []  
            TTT_lis_4_resid = []

            start , end = startend

            start = geok.tup_or_lis2dt(start)
            end   = geok.tup_or_lis2dt(end)
            
            ssss  = start
            eeee  = end
            
            ssssp = geok.dt2posix(start)
            eeeep = geok.dt2posix(end)
            
            if with_smart_ssp:
                Sdic = genefun.pickle_loader(smart_SSP_Sdic_path)
                ref_time = ssss + ((eeee - ssss)/2)
                print('INFO : smart SSP time' , ref_time)
                tnear , inear = gf.find_nearest(list(Sdic.keys()),ref_time)
                Z , C = Sdic[tnear][:,0] , Sdic[tnear][:,1]
                Z,C = ssp.SSP_extrapolate(Z,C,3000,50)

                
            for ttt,xbato in zip(TTT_lis,Xbato_lis_orig):
                if MNOPfile != 'O':
                    _ , outxbato  = geok.time_win_basic(ssssp,eeeep,ttt,xbato)
                else:
                    _ , outxbato  = geok.time_win_basic(ssssp,eeeep,ttt[-1],xbato)
                    outxbato = [tuple(e2) for e2 in outxbato]
                Xbato_lis.append(outxbato)
            for iiiddd , (ttt,obsasm) in enumerate(zip(TTT_lis,ObsASM_lis_orig)):
                if MNOPfile != 'O':
                    epocoutobsasm , outobsasm = geok.time_win_basic(ssssp,eeeep,ttt,obsasm)
                else:
                    epocoutobsasm , outobsasm = geok.time_win_basic(ssssp,eeeep,ttt[-1],obsasm)
                    
                print("avant/apres time win. " , len(obsasm) , len(outobsasm))
                if len(outobsasm) == 0:
                    print("ERR : no more valid pings after windowing, check the time window !!!!")
                    
                ObsASM_lis.append(outobsasm)
                TTT_lis_4_resid.append(epocoutobsasm)
                
        else: # not with time win
            TTT_lis_4_resid = []
            for iiiddd in range(len(TTT_lis)):
                TTT_lis_4_resid.append(TTT_lis[iiiddd][0])

            
        ObsASMgoodbool = []
        for obsasm in ObsASM_lis:
            if not with_jackknife:
                ObsASMgoodbool.append( np.array(len(obsasm) * [True] )) 
            else:
                RState = np.random.RandomState(jk_rand_seed)
                rand_vals = RState.rand(len(obsasm))
                bool_arr  = rand_vals < keep_ratio 
                if with_invert_jk:
                    bool_arr = np.logical_not(bool_arr)
                ObsASMgoodbool.append( bool_arr )
           
        A_stk = []
   
        # ===============================
        # BOUCLE D'INVERSION
        # ===============================
        print("START")
        print("===================== START =====================")
        start = genefun.get_timestamp(0)
        acls.print_n_dicwrit( "start" , str(start) , iterdic , 1 )
        acls.print_n_dicwrit("name" , exp_name , iterdic , 1)
        acls.print_n_dicwrit( "path" , exp_path , iterdic , 1)
        acls.print_n_dicwrit( "plateform" , gf.get_computer_name() , iterdic , 1)
        acls.print_n_dicwrit( "version python" , sys.version , iterdic , 1)
        acls.print_n_dicwrit( "version numpy" , np.version.version , iterdic , 1)
        acls.print_n_dicwrit( "version scipy" , scipy.version.version, iterdic , 1)
        acls.print_n_dicwrit("batch_mode" , bool(batch_mode) , iterdic , 1)
        acls.print_n_dicwrit("force_mode" , bool(force_mode) , iterdic , 1)

        with_var_lis = sorted([e for e in [e for e in list(globals().keys()) if type(e) is str] if e.startswith('with_')])
        
        for withvar in with_var_lis:
            acls.print_n_dicwrit( withvar , bool(globals()[withvar]) , iterdic,1)

        acls.print_n_dicwrit("weight ASM",sigma_defo_asm     , iterdic )
        acls.print_n_dicwrit("weight BL", sigma_defo_bl      , iterdic ) 
        acls.print_n_dicwrit("weight dZ", sigma_defo_zmaster , iterdic ) 
            
        acls.print_n_dicwrit("noise PXP apriori"     , sigma_pxp_apri  , iterdic )
        acls.print_n_dicwrit("seed noise PXP apriori", kmeta_pxp_apri  , iterdic )
        acls.print_n_dicwrit("noise dZ"              ,  sigma_dz_apri  , iterdic )
        acls.print_n_dicwrit("seed noise dZ"         , k_dz_apri       , iterdic )
        acls.print_n_dicwrit("noise  Z"              , sigma_z_apri    , iterdic ) 
        acls.print_n_dicwrit("seed noise Z"          , k_z_apri        , iterdic ) 
        acls.print_n_dicwrit("noise BL"              , sigma_bl_apri   , iterdic )
        acls.print_n_dicwrit("seed noise BL"         , k_bl_apri       , iterdic )

        if with_time_window:
            acls.print_n_dicwrit( "time window start" , (ssss), iterdic,1)
            acls.print_n_dicwrit( "time window end  " , (eeee), iterdic,1)
        if batch_mode:
            acls.print_n_dicwrit( "variables params" , bool_kw_lis , iterdic,1)
            acls.print_n_dicwrit( "var. params of the exp." , booldic , iterdic,1)
            
        if with_jackknife:
            acls.print_n_dicwrit( "jackknife inverted" , with_invert_jk , iterdic,1)
            acls.print_n_dicwrit( "jackknife random seed" , jk_rand_seed , iterdic,1)
            acls.print_n_dicwrit( "keep_ratio" , keep_ratio , iterdic,1)

        print("")
        
        iiter = 1
        
        # stop criteras 'kriter(old)' are defined as big enought values
        # for the initialisation and will be replaced after the 1st iteration
        kriter    = 10**7
        kriterold = 10**10
        
#       while np.linalg.norm(PXPapri - PXPold) > 5 * 10**-5 and iiter <= iitermax:
        # ITERATIONS
        while np.linalg.norm(kriter) > dX_sum_stop and iiter <= iitermax:
        #while np.linalg.norm(kriter - kriterold) > 10**-4:
            print("------------ Iteration No" , iiter , "------------")
            iterdic = acls.get_iterdic_in_expdic(expdic,iiter)
        
            # Partie ASM
            if with_dzcst:
                dz_cst = PXPdZ_lis
            else:
                dz_cst = None
                
            ###### Determination of B : observation vector
            
#            # ligne adel
#            Xbato_lis2 = []
#            for xbatooo in Xbato_lis:
#                Xbato_lis2.append([e[0] for e in xbatooo])
#            Xbato_lis = Xbato_lis2
#            # XXXXXXXXXXXXXXXXXXXXXXXXXXX
            
            ObsASM , ModASM  , arglis = acls.vectorialize_ASM_multi(PXPapri_lis,
                                                           ObsASM_lis ,
                                                           Xbato_lis,Z,C,
                                                           nbprocs=nbproc,
                                                           dz_cst=dz_cst,
                                                           ObsASMBoolinp_lis=ObsASMgoodbool)
            
            
            # tambouille zone   
#            pos_fb = (array([  2.61612797e+02,  -1.17006982e+02,  -1.28876579e-01]), array([  2.62098541e+02,  -1.16621655e+02,   1.68159868e-01]))
#            PXP = [  -64.44905644 ,  727.96098426 , 1621.6416359 ]
#            
#            RTtup_lis3 = []
#            RTtup_lis4 = []
#
#            for i in range(len(Xbato_lis)):
#
#                for j in range(len(Xbato_lis[i])):
#                    
#                    if ObsASMgoodbool[i][j]:
#                        RTtup  = rt.raytrace_seek(Xbato_lis[i][j],PXPapri_lis[i], Z, C,
#                                                  thetaminin=0, thetamaxin=88,
#                                                  verbose=0,fulloutput=0,severe=False)
#                        RTtup_lis3.append(RTtup)
#                        
#                    RTtup4  = rt.raytrace_seek(Xbato_lis[i][j],PXPapri_lis[i], Z, C,
#                                              thetaminin=0, thetamaxin=88,
#                                              verbose=0,fulloutput=0,severe=False)
#                    RTtup_lis4.append(RTtup4)                        
#
#            RTtup_lis = []
#            for e in arglis:
#                RTtup_lis.append(rt.raytrace_seek(*e))
#
#            ObsASM[-1] , ModASM[-1]
#            
#            raise Exception
            
            
            Nb_nan = np.sum(np.isnan(ModASM))
            
            # case where NaN are found in the Model
            if Nb_nan != 0:
                print('WARN : ' , Nb_nan , 'NaN in ModASM (Sub Marine Acoustic LSQ Model)')
                print('those pings will be removed and a new Model will be determined')
                if iiter == 1 and False: # mode for the 1st iteration DISCONTINUED
                    L_ObsASMgoodbool =[len(e) for e in ObsASMgoodbool]
                    NaN_ModASM = np.isnan(ModASM)
                    NaN_ModASM_sublist = gf.sublistsIt(NaN_ModASM,L_ObsASMgoodbool,output_array=True)
                    ObsASMgoodbool_nan = []
                    
                    for Nanasm , Bool in zip(NaN_ModASM_sublist , ObsASMgoodbool):
                        ObsASMgoodbool_nan.append(np.logical_and(Bool,np.logical_not(Nanasm)))
                        
                    ObsASMgoodbool = ObsASMgoodbool_nan
                        
                else: # mode for the >1 iteration
                    ObsASMgoodbool_stacked = np.hstack(ObsASMgoodbool)
                    ObsASMgoodbool_stacked[ObsASMgoodbool_stacked == True] = np.logical_not(np.isnan(ModASM))
                    L_ObsASMgoodbool =[len(e) for e in ObsASMgoodbool]
                    ObsASMgoodbool = gf.sublistsIt(ObsASMgoodbool_stacked,L_ObsASMgoodbool,output_array=True)
                    
                if iiter > 1 and with_v_4_p_reinject:
                    V_4_P_reinject = np.array(V_4_P_reinject)[np.logical_not(np.isnan(ModASM))]
                 
                    
                
                # Recalc of the observation vector without NaN
                ObsASM , ModASM  , arglis = acls.vectorialize_ASM_multi(PXPapri_lis,
                                                               ObsASM_lis ,
                                                               Xbato_lis,Z,C,
                                                               nbprocs=nbproc,
                                                               dz_cst=dz_cst,
                                                               ObsASMBoolinp_lis=ObsASMgoodbool)
            
            
        
            
            # ============= TEST ZONE ================================
            #import raytrace 
            ##Xbato_lis[0][0][0]
            #raytrace.raytrace_seek(Xbato_lis[0][0],PXPapri_lis[0],Z,C,
            #                       verbose=0,fulloutput=0)
            #A = arglis[0]
            #raytrace.raytrace_seek(*A)
            # ============= TEST ZONE ================================
            
            B_ASM = ObsASM - ModASM
            
            if plot_mode:
                vertbar_pings = np.cumsum([np.sum(e) for e in ObsASMgoodbool])
                
                figObsMod , axObsMod = plt.subplots()
                figObsMod.suptitle('Acoustics Obs. and Model, iter. ' + str(iiter) + ' (red bars split each beacon data)' )

                axObsMod.plot(ObsASM,'bx',label='Observations')
                axObsMod.plot(ModASM,'g+',label='Model (Ray Tracing)')
                
                for vbar in vertbar_pings:
                    axObsMod.axvline(vbar,color='r')
                
                figObsMod.set_size_inches(11.69,8.27)

                axObsMod.set_xlabel('Pings')
                axObsMod.set_ylabel('Two way travel time (seconds)')
                axObsMod.legend()
                

                plt.savefig(protopath + 'iter' + str(iiter) + '_ObsMod.pdf')
                plt.savefig(protopath + 'iter' + str(iiter) + '_ObsMod.png')

                figB      , axB      = plt.subplots()
                figB.suptitle('Difference between Obs. & Model (Ray Tracing), iter. ' + str(iiter) + ' (red bars split each beacon data)' )

                axB.plot(B_ASM,'.',c='orange')
                for vbar in vertbar_pings:
                    axB.axvline(vbar,color='r')
                axB.set_xlabel('Pings', size=15)
                axB.set_ylabel('Residuals (seconds)', size=15)
                figB.set_size_inches(11.69,8.27)
                
                axB.tick_params(axis='both', which='major', labelsize=15)
                axB.tick_params(axis='both', which='minor', labelsize=11)

                plt.savefig(protopath + 'iter' + str(iiter) + '_B.pdf')
                plt.savefig(protopath + 'iter' + str(iiter) + '_B.png')
                plt.savefig(protopath + 'iter' + str(iiter) + '_B.svg')
                

                plt.close('all')

            
            ###### Determination of A : Jacobian matrix in differents cases
            if with_zmaster:
                JacobASM = acls.jacob_ASM((PXPref,dPXPapri_lis), ObsASM_lis,
                                          Xbato_lis,Z,C,h,
                                          nbprocs=nbproc,monoZ=with_mono_z,
                                          accur=1,dz_cst=dz_cst,
                                          ObsASMBoolinp_lis=ObsASMgoodbool)
            elif not with_barycenter:
                JacobASM = acls.jacob_ASM(PXPapri_lis,ObsASM_lis,Xbato_lis,Z,C,h,
                                              nbprocs=nbproc,monoZ=with_mono_z,
                                              accur=1,dz_cst=dz_cst,
                                              ObsASMBoolinp_lis=ObsASMgoodbool)
            else:
                JacobASM = acls.jacob_ASM((PXPbary,dPXPapri_lis),ObsASM_lis,Xbato_lis,Z,C,h,
                                              nbprocs=nbproc,monoZ=with_mono_z,
                                              accur=1,dz_cst=dz_cst,
                                              ObsASMBoolinp_lis=ObsASMgoodbool)  
            # debug bloc
            #rt.raytrace_diff_light(*JacobASM[1][0])
            #for argsjacobasm in JacobASM[1][0]:
            #    rt.raytrace_diff_light(*argsjacobasm)
            
            # les baselines sont injectées comme observables => mode contraint et non mode fixé
            #Partie BL
            if with_bl:
                if with_barycenter:
                    ObsBL,ModBL = acls.vectorialize_BL(BLnoised,dPXPapri_lis,dz_cst)
                else:
                    ObsBL,ModBL = acls.vectorialize_BL(BLnoised,PXPapri_lis,dz_cst)
                B_BL    = ObsBL - ModBL
                JacobBL = acls.jacob_BL(PXPapri_lis,with_mono_z,with_barycenter,dz_cst)
                
            if with_zmaster:
                # Vieux bloc foireux (mélange des PXP observés/theoriques p/r au pxp master/ barycentre)
                # On le garde en souvenir (161009)
                #ObsZ = PXPdZ_lis
                #ModZ = np.array([dpxp[-1] for dpxp in dPXPapri_lis])
                #ObsZ = np.delete(ObsZ,i_pxp_master)
                #ModZ = np.delete(ModZ,i_pxp_master)


                ObsZ = np.delete(PXPdZ_lis,i_pxp_master)
                ModZ = np.delete(dPXPapri[2::3] ,i_pxp_master) -  dPXPapri[2::3][i_pxp_master]
                
                B_Z  = ObsZ - ModZ
                
                nz = len(PXPdZ_lis)
                JacobZ = np.hstack((np.zeros((nz,JacobBL.shape[1] - nz)) , np.diag(np.ones((nz,nz))[:,0])))
                JacobZ = np.zeros(JacobZ.shape)
                for i,lin in enumerate(JacobZ):
                    lin[3+ i*3+2] = 1  
                    lin[3+ i_pxp_master*3+2] = -1  

                JacobZ = np.delete(JacobZ , i_pxp_master , 0 )
                nz -= 1 


# BLOC souvenir lorsque les obs de dZ ne
# se rapportent au bary et pas au pxp master
#                    B_Z  = ObsZ - ModZ
#                    
#                    nz = len(PXPdZ_lis)
#                    JacobZ = np.hstack((np.zeros((nz,JacobBL.shape[1] - nz)) , np.diag(np.ones((nz,nz))[:,0])))
#                    JacobZ = np.zeros(JacobZ.shape)
#                    for i,lin in enumerate(JacobZ):
#                        lin[3+ i*3+2] = 1  
#                        lin[2] = -1  

                    
            # Partie Commune

            print("INFO : End of Jacobian filling, start of redesign + weight management")
            A = JacobASM
            B = B_ASM
            

            #P = gf.diagonalize(0.001,np.max(A.shape))
            # Bloc des poids bien pas beau comme il faut 
            # réformé le 160603
            #if not with_v_4_p_reinject or not "V_4_P_reinject" in globals().keys():
            #    Ptmp = [1 / ((10**-3)**2)] * np.max(A.shape)
            #else:
            #    Ptmp = 1 / ((np.array(V_4_P_reinject) * 1) ** 2)
            #
            #if with_bl:
            #    #_,_,Ptmp = geok.weight_mat([10,1],[len(ObsASM),len(ObsBL)],
            #    #                        fuv,sparsediag=True)
            #    #Ptmp = gf.diagonalize(0.001,len(ObsBL))
            #    Ptmp = Ptmp + [0.001] * len(ObsBL)
            #
            #    #P = scipy.linalg.block_diag(P,Ptmp)
            #    A = np.vstack((A,JacobBL))
            #    B = np.hstack((B,B_BL))
            #    
            #if with_zmaster:
            #    #Ptmp = gf.diagonalize(10**0,nz)
            #    Ptmp = Ptmp + [10**0] *  nz
            #
            #    #P = scipy.linalg.block_diag(P,Ptmp)
            #    A = np.vstack((A,JacobZ))
            #    B = np.hstack((B,B_Z))
            #    
            #P = sprs.diags(Ptmp,0)
            #
            #del Ptmp

            sigma_4_P = []
            len_4_P   = []
            
            if not with_v_4_p_reinject or not "V_4_P_reinject" in list(globals().keys()): # V_4_P_reinject (2nd teste) désigne la liste et non le bouleen, pas touche
                sigma_4_P = sigma_4_P + [sigma_defo_asm]
                len_4_P.append(np.max(A.shape))
            else:
                sigma_4_P = sigma_4_P   + list(V_4_P_reinject)
                len_4_P   = len_4_P     +  len(V_4_P_reinject) * [1]

            if with_bl:
                #_,_,Ptmp = geok.weight_mat([10,1],[len(ObsASM),len(ObsBL)],
                #                        fuv,sparsediag=True)
                #Ptmp = gf.diagonalize(0.001,len(ObsBL))
                sigma_4_P = sigma_4_P + [sigma_defo_bl]
                len_4_P.append(len(ObsBL))
                A = np.vstack((A,JacobBL))
                B = np.hstack((B,B_BL))             

            if with_zmaster:
                #Ptmp = gf.diagonalize(10**0,nz)
                sigma_4_P = sigma_4_P + [sigma_defo_zmaster]
                len_4_P.append(nz)                   
                
            
                #P = scipy.linalg.block_diag(P,Ptmp)
                A = np.vstack((A,JacobZ))
                B = np.hstack((B,B_Z))

            _,_,P = geok.weight_mat(np.array(sigma_4_P),len_4_P,fuv,sparsediag=True)
            #_,_,P = geok.weight_mat(np.array(sigma_4_P),len_4_P,fuv,sparsediag=True)

            print("INFO : Weights end")

#                if not batch_mode:        
#                    A_stk.append(A)
            
            # la contrainte de barycentre est injectée suivant la methode fixée d'Helmert
            if with_barycenter:
                G     = acls.constraints_matrix(len(PXP_lis),with_barycenter,with_mono_z) #MODIF
                #G2    = np.zeros(np.max(G.shape)) #MODIF
                #G2[2] = 1               #MODIF
                #G = np.vstack((G,G2))   #MODIF
                n_lagrangian = 3         #MODIF


            # ==== A partir de là on fat les conv memap / sparse
#                A = np.matrix(A)
#                B = np.matrix(B).T
#                P = np.matrix(P)
             
            A = genefun.memmap_from_array(A)
            B = genefun.memmap_from_array(B)
            #P = genefun.memmap_from_array(P)
#                N = A.T * P * A
            
            print("INFO : Sparse matrix conversion")
            Asprs = sprs.csc_matrix(A)
            Bsprs = sprs.csc_matrix(B).T
            Psprs = sprs.csc_matrix(P)
            print("INFO : End of sparse conversion")
            
            Nsprs = (Asprs.T).dot(Psprs).dot(Asprs)
            N     = Nsprs.toarray()
        
#                N = np.dot(np.dot(A.T , P ) , A)
            
            # Liquidation des matrices inutiles
            JacobASM , ObsASM = [] , []

            # debut contraintes
            if with_barycenter:
                print("INFO : Creation of the constraint matrix for barycenter determination")
                Gsiz = G.shape[0]
                O = np.zeros((G.shape[0],G.shape[0]))
                N = np.vstack((N,G))    
                N = np.hstack((N,np.hstack((G,O)).T))
            # fin contraintes       
            Ninv = scipy.linalg.inv(N) 
            print("INFO : Normal matrix inversion done")

#                AtPB = A.T * P * B
#                AtPB = np.dot(np.dot(A.T , P ) , B)  

            print("INFO : Normal equation is about to be solved")
            AtPBsprs = (Asprs.T).dot(Psprs).dot(Bsprs)
            AtPB     = AtPBsprs.toarray()
            Corel = np.corrcoef(Ninv)
            
            if with_barycenter:
                AtPBbis     = np.vstack((AtPB,np.matrix(np.zeros(Gsiz)).T)) 
                #AtPBbis[-1] = -63.524075746381186 #acls.barycenter_calc(PXPtrue_arr[:,-1]) # LIGNE WRONG (160610)
                dX = np.dot( Ninv , AtPBbis )
            else:
                dX = np.dot( Ninv , AtPB )
            

        #    c,resid,rank,sigma = scipy.linalg.lstsq(A,B)
            # ==== fin de la section des matrix  
            
            dX = np.array(dX).squeeze()
            dXorig = dX
            print("INFO : Normal equation is solved")

                    
            if with_barycenter:
                if not with_mono_z:        
                    dX = dX[:-3]
                else:
                    dX = dX[:-2]
                    
            debugrkmode = 0 # Discontinué
            
            if with_zmaster and debugrkmode:
                # manip empirique de **** **** pour corriger la non nullité
                # de la Somme des dZ ... (160609)
                # On remarque empiriquement que sum dZ = une valeur dZcorrection
                # que l'on ajoute ensuite aux dZ individuel
                # Le debugrk mode est discontinué mais a grandement aidé 
                # au debugging, on le garde en souvenir (160610)
                PXPbarynew_debugrk = np.array(PXPbary) + dX[:3]
                PXPbary_debugrk    = np.array(PXPbarynew_debugrk) 
                
                dPXPnew_debugrk  = dPXPapri + dX[n_lagrangian:]
                bary_debugrk     = dPXPnew_debugrk.reshape(shape_arr)
                dZcorr_debugrk   = acls.barycenter_calc(bary_debugrk)[-1]
                  
                for i in range(nPXP):
                    dX[ 3+ i*3+2 ] = dX[ 3+ i*3+2 ] - dZcorr_debugrk

            print("INFO : Creating the new PXP coordinates")
      
            if with_mono_z:
                dX = np.insert(dX,np.arange(2,len(dX)-2,2),dX[-1]) # manip sioux pour ajouter les dZ communs à tous les PXP 
            if not with_barycenter:
                PXPold  = PXPapri
                PXPnew  = PXPapri + dX
                PXPapri = PXPnew
                PXPapri_lis_old = PXPapri_lis
                PXPapri_lis     = list(PXPnew.reshape(shape_arr))
            else:
                PXPbaryold = np.array(PXPbary)
                PXPbarynew = np.array(PXPbary) + dX[:3]
                PXPbary    = np.array(PXPbarynew) 

                dPXPnew      = dPXPapri + dX[n_lagrangian:]
                dPXPapri     = dPXPnew
                dPXPapri_lis = list(dPXPnew.reshape(shape_arr))    

                PXPbarynew_tmp    = np.array(PXPbarynew)
                #PXPbarynew_tmp[2] = PXPbaryold[2] # WHY CETTE LIGNE !!!!
                #OK c'était un relicat d'une tentative de resolution du bug
                #chronique de non nullité des dZ (160609)
                PXPnew = dPXPnew + np.tile(PXPbarynew_tmp,nPXP)


                PXPapri_lis_old = PXPapri_lis
                PXPapri_lis     = list(PXPnew.reshape(shape_arr))

#                ObsASM_V , ModASM4_V  = acls.vectorialize_ASM_multi(PXPapri_lis,
#                                                               ObsASM_lis,
#                                                               Xbato,Z,C,
#                                                               nbprocs=nbproc,
#                                                               dz_cst=dz_cst)
#                                                               
#                                                               
#                Vnaif = ObsASM_V - ModASM4_V
#                
#                if with_bl:
#                    if with_barycenter:
#                        ObsBL_V,ModBL_V = acls.vectorialize_BL(BL,dPXPapri_lis)
#                    else:
#                        ObsBL_V,ModBL_V = acls.vectorialize_BL(BL,PXPapri_lis)
#                        
#                    V_BL = ObsBL_V - ModBL_V
#                    Vnaif = np.concatenate((Vnaif,V_BL))

            Vrigour = B - A.dot(dXorig[:A.shape[1]])
            
            V = Vrigour
            
            
            if with_bl:
                V_for_cleaning = V[:-len(B_BL)]
            else:
                V_for_cleaning = V


            if cleaning_coef != 0 and not ObsASMgoodbool is None: # pour quoi cleaning_coef != 0, et pourquoi ici ? chépa (160922)
                print("INFO : Cleaning pings based on residuals")                
                lenObsASM = [np.sum(o) for o in ObsASMgoodbool]
                V_per_ObsASM = gf.sublistsIt(V_for_cleaning,lenObsASM,True)

                ObsASMgoodbool_tmp = []
                Vcleaned = []
                if with_v_4_p_reinject:
                    V_4_P_reinject = []
                if with_cleaning: 
                    for iii  , (boolping  , Voasm) in enumerate(zip(ObsASMgoodbool , V_per_ObsASM )):
                        print("INFO : Valid pings BEFORE cleaning",np.sum(boolping),'for PXP',iii)
                        indices_true  = np.where(boolping)[0]
                        indices_false = np.where(np.logical_not(boolping))[0]
                                                
                        actualised_bool = boolping[indices_true] * (np.abs(np.array(Voasm)) < cleaning_coef * np.std(Voasm))
                        boolping_new = np.array(boolping)
                        boolping_new[indices_true] = actualised_bool
                            
                        ObsASMgoodbool_tmp.append(boolping_new)
                        
                        Vcleaned = Vcleaned + list(Voasm[actualised_bool])

                        if with_v_4_p_reinject:
                            V_4_P_reinject = V_4_P_reinject + list(Voasm[actualised_bool])
                        print("INFO : Valid pings AFTER cleaning",np.sum(boolping_new),'for PXP',iii)
                   
                    ObsASMgoodbool = ObsASMgoodbool_tmp

                if with_v_4_p_reinject:
                    V_4_P_reinject = np.array(V_4_P_reinject)
                
            oldfuv = fuv

            fuv = geok.fuv_calc(V,A,P)
            fuv_noP = geok.fuv_calc(V,A)
                
            with np.errstate(invalid='ignore'):
                sigmas = np.sqrt(np.diag(Ninv) * fuv)
                sigmas2 = np.sqrt(np.diag(Ninv) *100000000000)

            if with_barycenter:
                if not with_mono_z:        
                    sigmas = sigmas[:-3]
                else:
                    sigmas = sigmas[:-2] 
                    
            iterdic['varcovar'] = Ninv
            iterdic['V']        = V
            iterdic['B_ASM']    = B_ASM
            iterdic['B']        = B
            
            print("")
            print("----- Results of the iteration " , iiter)
            print("")
            print("how to read 'dX' & 'sigmas' fields : ")
            if with_barycenter:
                print('1st row = results for Barycenter')
                print('2nd row = results for PXP1')
                print('3rd row = results for PXP2')
                print('etc ...')
                print('Thus: New PXP Coords = Apri PXP Coords + dXbary (1st row) + dX_PXP')
                
            else:
                print('1st row = results for PXP1')
                print('2nd row = results for PXP2')
                print('etc ...')
            print("")
            
            if with_barycenter:
                sigmas_reshape = np.reshape(sigmas,(nPXP + 1,3))
                dX_reshape     = np.reshape(dX,(nPXP + 1,3))

            else:
                sigmas_reshape = np.reshape(sigmas,(nPXP,3))
                dX_reshape     = np.reshape(dX,(nPXP,3))

                
            acls.print_n_dicwrit('f.u.v.' , fuv , iterdic , 1)
            acls.print_n_dicwrit('f.u.v. before inversion' , oldfuv , iterdic , 1)
            acls.print_n_dicwrit('coords. apriori' , np.array(PXPapri_lis_old) ,
                                 iterdic )

            acls.print_n_dicwrit('barycenter apriori (raw calc)' , 
                                 acls.barycenter_calc(np.array(PXPapri_lis_old)),            
                                 iterdic )

            acls.print_n_dicwrit('dX' , dX_reshape , iterdic )
            acls.print_n_dicwrit('sum dX' , np.sum(dX) , iterdic , 1)
            acls.print_n_dicwrit('sum abs dX' , np.sum(np.abs(dX)) , iterdic , 1)
             
            PXPnew_arr = PXPnew.reshape(shape_arr)
           
            acls.print_n_dicwrit( 'new coords.' , PXPnew_arr , iterdic )
            if not batch_mode:
                PXPnew_stk.append(PXPnew)

            acls.print_n_dicwrit( 'sigmas' , sigmas_reshape , iterdic )
            
            barybrut = acls.barycenter_calc(PXPnew.reshape(shape_arr))            
            barybrut_old = acls.barycenter_calc(PXPold.reshape(shape_arr))   
            
            kriterold = kriter
            kriter = np.linalg.norm(barybrut_old - barybrut)
            kriter = np.sum(np.abs(dX))
            
            acls.print_n_dicwrit( "critera/critera old" , [ kriter , kriterold ] ,
                                 iterdic)
            
            if with_print_simu_infos:
                acls.print_n_dicwrit( "delta 3D to the true position in distance (m) for each PXP" , 
                                np.linalg.norm( PXPnew_arr - PXPtrue_arr,axis=1) ,
                                iterdic)
                acls.print_n_dicwrit( "delta 2D to the true position in distance (m) for each PXP" , 
                                np.linalg.norm((PXPnew_arr - PXPtrue_arr)[:,:2],axis=1) ,
                                iterdic)
    
                acls.print_n_dicwrit( "delta of coordinates to the true position (m) for each PXP"  , 
                                PXPnew_arr - PXPtrue_arr , iterdic ) 
            
            if with_print_simu_infos and with_bl:
                acls.print_n_dicwrit( "true baselines, from no-noised coords. (only meanful for simu)" , 
                                     acls.print_BL(PXPtrue_arr) , iterdic)
            
            if with_bl:   
                acls.print_n_dicwrit( "observed baselines (as in B-file)" ,
                                     acls.print_BL(ObsBL,0) , iterdic)                

            acls.print_n_dicwrit("baselines of first apriori coordinates ('Mod')",            
                                 acls.print_BL(PXPapri0_arr),iterdic)
            
            acls.print_n_dicwrit("baselines of current iteration apriori coordinates ('Mod')",            
                                 acls.print_BL(np.vstack(np.vstack(PXPapri_lis_old))),iterdic)

            acls.print_n_dicwrit("new baselines after inversion",
                                 acls.print_BL(PXPnew_arr),iterdic)
            
            acls.print_n_dicwrit("raw barycentrer : Sum Xpxp / Npxp",
                                 barybrut , iterdic)

            baryvrai = acls.barycenter_calc(PXPtrue_arr)      

            if with_print_simu_infos:            
                acls.print_n_dicwrit("delta 3D bary raw/true in distance" , 
                np.linalg.norm(barybrut - baryvrai) , iterdic)

                acls.print_n_dicwrit("delta 2D bary raw/true in distance" ,
                np.linalg.norm(barybrut[0:2] - baryvrai[0:2]) , iterdic ) 
            
                acls.print_n_dicwrit("delta bary raw/true in coords.",
                barybrut - baryvrai , iterdic ) 

            if with_barycenter:
                acls.print_n_dicwrit("barycentrer estimated in LSQ",
                                     PXPbary , iterdic )

                acls.print_n_dicwrit("barycentrer dXs in LSQ",                
                acls.barycenter_calc(dPXPnew.reshape(shape_arr)) , iterdic)
                
                if np.sum(np.abs(acls.barycenter_calc(dPXPnew.reshape(shape_arr)))) > 10**-5:
                    print("WARN : sum Barycentrer dXs in LSQ > 10**-5 !!!!")

            if with_barycenter and with_print_simu_infos:                         
                acls.print_n_dicwrit("delta 3D bary LSQ/true in distance", 
                np.linalg.norm(PXPbary - baryvrai) , iterdic)
                
                acls.print_n_dicwrit("delta 2D bary LSQ/true in distance", 
                np.linalg.norm(PXPbary[0:2] - baryvrai[0:2]) , iterdic)
            
                acls.print_n_dicwrit("delta bary LSQ/true in coords.",
                PXPbary - baryvrai , iterdic )
            
            if with_zmaster:
                acls.print_n_dicwrit("delta to reference beacon",
                                     PXPnew_arr[:,2] - PXPnew_arr[i_pxp_master,2],
                                     iterdic)
                acls.print_n_dicwrit("delta in input to reference beacon",
                                     PXPdZ_lis,
                                     iterdic)                                         
                acls.print_n_dicwrit("delta true to reference beacon",
                                     PXPdZ_lis_true,
                                     iterdic)   

            if plot_mode:
                ####################### HISTO V NORMALIZED #######################
                
                Vnorma      = np.array(V / P.diagonal())
                # OLD CLEANING OBSOLETE ! (170626)
                Vnormclean  = Vnorma[np.abs(Vnorma) < 3 * np.std(Vnorma)]
                
                plt.clf()
                plt.subplots_adjust(left=0.2, right=0.9, top=0.9, bottom=0.2)
    
                n,bins,patches = plt.hist(Vnormclean ,100,normed=1)
                gauss = scipy.stats.norm.pdf(bins,np.mean(Vnormclean),np.std(Vnormclean))
                plt.plot(bins,gauss)
    
                plt.xlabel('TWTT residuals (s)')
    
                plt.savefig(protopath + 'iter' + str(iiter) + '_histVnorma.png')
                plt.savefig(protopath + 'iter' + str(iiter) + '_histVnorma.pdf')
                
                
                ####################### preliminar for TEMPORAL V ##########

                if len(Vcleaned) != 0:
                    Vhisto = np.array(Vcleaned)
                else:
                    Vhisto = V
                
                ####################### HISTO V REGULAR #######################
                
                plt.clf()
                plt.subplots_adjust(left=0.2, right=0.9, top=0.9, bottom=0.2)
    
                n,bins,patches = plt.hist(Vhisto ,100,normed=1)
                gauss = scipy.stats.norm.pdf(bins,np.mean(Vhisto),np.std(Vhisto))
                plt.plot(bins,gauss)
                
                plt.xlabel('TWTT residuals (s)')
                plt.ylabel('occurence')
                
                plt.savefig(protopath + 'iter' + str(iiter) + '_histV.png')
                plt.savefig(protopath + 'iter' + str(iiter) + '_histV.pdf')
          
                ####################### HISTO V REGULAR 3sig centered ##########
                
                plt.clf()
                plt.subplots_adjust(left=0.2, right=0.9, top=0.9, bottom=0.2)
    
                n,bins,patches = plt.hist(Vhisto ,100,normed=1)
                gauss = scipy.stats.norm.pdf(bins,np.mean(Vhisto),np.std(Vhisto))
                plt.plot(bins,gauss)
                plt.xlim((- np.std(Vhisto) * 3 , np.std(Vhisto) * 3))
                
                plt.xlabel('TWTT residuals (s)')
                plt.ylabel('occurence')
                
                plt.savefig(protopath + 'iter' + str(iiter) + '_histV_3sig_centered.png')
                plt.savefig(protopath + 'iter' + str(iiter) + '_histV_3sig_centered.pdf')
                
                
                
                
                ####################### preliminar for TEMPORAL V ##########
                lenObsASM_4_plot = [np.sum(o) for o in ObsASMgoodbool]
                if cleaning_coef != 0:
                    V_per_ObsASM_4_plot = gf.sublistsIt(Vcleaned,lenObsASM_4_plot,True)
                else:
                    V_per_ObsASM_4_plot = gf.sublistsIt(V,lenObsASM_4_plot,True)
                
                ####################### TEMPORAL V mk 1 ##########
                
                plt.clf()
                plt.subplots_adjust(left=0.2, right=0.9, top=0.9, bottom=0.2)
    
                for iidplot,(Arrbool , Ttt , Vv ) in enumerate(zip(ObsASMgoodbool,TTT_lis_4_resid,V_per_ObsASM_4_plot )):
                    plt.plot(geok.posix2dt(Ttt,1)[Arrbool],Vv,'.',label = 'ID' + str(idPXP_lis[iidplot]))
                
                    plt.ylabel('TWTT residuals (s)')
                    plt.legend()
                    
                fig = plt.gcf()
                fig.autofmt_xdate()
    
                plt.savefig(protopath + 'iter' + str(iiter) + '_temporalV.png')
                plt.savefig(protopath + 'iter' + str(iiter) + '_temporalV.pdf')
                

                ####################### TEMPORAL V mk 2 ########## 
                
                plt.clf()
                plt.subplots_adjust(left=0.2, right=0.9, top=0.9, bottom=0.2)
    
                fig , AXES = plt.subplots(len(idPXP_lis),1)
                
                if not genefun.is_iterable(AXES):
                    AXES = [AXES]
                
                fig.set_size_inches(8.27,11.69)
                
                smartV_stk = []
                COL = ['r','g','y','b']
                for iidplot,(Arrbool , Ttt , Vv , ax , col) in enumerate(zip(ObsASMgoodbool,TTT_lis_4_resid,V_per_ObsASM_4_plot , AXES,COL)):
                    ax.plot(geok.posix2dt(Ttt,1)[Arrbool],Vv,'.',label = 'ID' + str(idPXP_lis[iidplot]),color=col)
                    smartV_stk.append(np.column_stack(([iidplot] * len(Vv),Ttt[Arrbool],Vv)))
                
                    ax.set_ylabel('TWTT residuals (s)')
                    ax.legend()
    
                fig.autofmt_xdate()
    
                plt.savefig(protopath + 'iter' + str(iiter) + '_temporalV2.png')
                plt.savefig(protopath + 'iter' + str(iiter) + '_temporalV2.pdf')

                #### Saving associated file
                np.savetxt(protopath + 'iter' + str(iiter) + '.smartV' ,
                           np.vstack(smartV_stk) )            
            
            np.savetxt(protopath  + 'iter' + str(iiter) + '.V' , V)
            
            print("")
            iiter=iiter+1   

        print("===================== FINAL =====================")
        end = genefun.get_timestamp(0)
        acls.print_n_dicwrit("end"   , str(end) , iterdic , 1 )
        acls.print_n_dicwrit("duration" , str(end - start) , iterdic , 1 )
        acls.print_n_dicwrit("name" , exp_name , iterdic , 1)
        acls.print_n_dicwrit( "path" , exp_path , iterdic , 1)
        acls.print_n_dicwrit( "plateform" , gf.get_computer_name() , iterdic , 1)
        acls.print_n_dicwrit( "version python" , sys.version , iterdic , 1)
        acls.print_n_dicwrit( "version numpy" , np.version.version , iterdic , 1)
        acls.print_n_dicwrit( "version scipy" , scipy.version.version, iterdic , 1)
        acls.print_n_dicwrit("batch_mode" , bool(batch_mode) , iterdic , 1)
        acls.print_n_dicwrit("force_mode" , bool(force_mode) , iterdic , 1)
        
        with_var_lis = sorted([e for e in [e for e in list(globals().keys()) if type(e) is str] if e.startswith('with_')])
        
        for withvar in with_var_lis:
            acls.print_n_dicwrit( withvar , bool(globals()[withvar]) , iterdic,1)

        acls.print_n_dicwrit("weight ASM",sigma_defo_asm     , iterdic )
        acls.print_n_dicwrit("weight BL", sigma_defo_bl      , iterdic ) 
        acls.print_n_dicwrit("weight dZ", sigma_defo_zmaster , iterdic ) 
            
        acls.print_n_dicwrit("noise PXP apriori"     , sigma_pxp_apri  , iterdic )
        acls.print_n_dicwrit("seed noise PXP apriori", kmeta_pxp_apri  , iterdic )
        acls.print_n_dicwrit("noise dZ"              ,  sigma_dz_apri  , iterdic )
        acls.print_n_dicwrit("seed noise dZ"         , k_dz_apri       , iterdic )
        acls.print_n_dicwrit("noise  Z"              , sigma_z_apri    , iterdic ) 
        acls.print_n_dicwrit("seed noise Z"          , k_z_apri        , iterdic ) 
        acls.print_n_dicwrit("noise BL"              , sigma_bl_apri   , iterdic )
        acls.print_n_dicwrit("seed noise BL"         , k_bl_apri       , iterdic )

        if with_time_window:
            acls.print_n_dicwrit( "time window start" , (ssss), iterdic,1)
            acls.print_n_dicwrit( "time window end  " , (eeee), iterdic,1)
        if batch_mode:
            acls.print_n_dicwrit( "variables params" , bool_kw_lis , iterdic,1)
            acls.print_n_dicwrit( "var. params of the exp." , booldic , iterdic,1)
            
        if with_jackknife:
            acls.print_n_dicwrit( "jackknife inverted" , with_invert_jk , iterdic,1)
            acls.print_n_dicwrit( "jackknife random seed" , jk_rand_seed , iterdic,1)
            acls.print_n_dicwrit( "keep_ratio" , keep_ratio , iterdic,1)

            
        acls.print_n_dicwrit("Nb pings     " , len(Xbato) , iterdic , 1)
        acls.print_n_dicwrit("Nb PXPs      " , nPXP , iterdic , 1)
        acls.print_n_dicwrit("Size Jacob." , A.shape , iterdic , 1)
        acls.print_n_dicwrit("Size Resid." , V.shape , iterdic , 1)
        acls.print_n_dicwrit("Nb iter.     " , iiter-1  , iterdic , 1)
        
        F.stop()
        
        expdic[-1] = expdic[max(expdic.keys())]
        
        genefun.save_obj_as_file(expdic , exp_path_output , exp_name, suffix=bool_4_exp_str)
                                    
        geok.chi2_test_lsq(V,A,P)
        
        #%%
        if plot_mode:
            plt.clf()
            if with_fb_mode:
                if with_time_window:
                    Xplot = np.array([e[0] for e in Xbato_lis[-1]])[:,0][ObsASMgoodbool[-1]]
                    Yplot = np.array([e[0] for e in Xbato_lis[-1]])[:,1][ObsASMgoodbool[-1]]
                else:
                    Xplot = np.vstack(Xbato_f)[:,0][ObsASMgoodbool[-1]]
                    Yplot = np.vstack(Xbato_f)[:,1][ObsASMgoodbool[-1]]
                plt.plot(Xplot,Yplot, 'b-')
                #plt.plot(np.vstack(Xbato_b)[:,0], np.vstack(Xbato_b)[:,1], 'rx')
            else:
                print("WARN : plot traject may be wrong")
                plt.plot(zip(*Xbato)[0] , zip(*Xbato)[1] , 'x')
            plt.scatter(PXPnew_arr[:,0] , PXPnew_arr[:,1],s=200)
            
            for (pxp,(i,j)) in zip(PXPnew_arr,itertools.combinations(list(range(nPXP)),2)):
                print(i,j)
                
            for ipxp,pxp in enumerate(PXPnew_arr):
                if with_barycenter:
                    kkk = 3
                else:
                    kkk = 0
                if 1: # le bloc de test pour reproduire les ellipses
                    # Les sigma sont normalisé par le fuv (B) ou non (A)
                    sigxA  = np.sqrt(Ninv[ kkk + ipxp * 3    , kkk + ipxp * 3   ]) 
                    sigyA  = np.sqrt(Ninv[ kkk + ipxp * 3 +1 , kkk + ipxp*3 +1  ]) 
                    sigxyA = Ninv[kkk + ipxp * 3,kkk + ipxp*3 +1] 
                    sigxB  = np.sqrt(Ninv[ kkk + ipxp * 3    , kkk + ipxp * 3   ] * fuv) 
                    sigyB  = np.sqrt(Ninv[ kkk + ipxp * 3 +1 , kkk + ipxp*3 +1  ] * fuv) 
                    sigxyB = Ninv[kkk + ipxp * 3,kkk + ipxp*3 +1] * fuv
                    
                    # C'est la fonction que l'on utilise jusqu'a présent
                    xe1,ye1,_,_ = geok.error_ellipse(pxp[0],pxp[1], sigxB , sigyB , sigxyB, scale= 10000) 
                    xe2,ye2,_,_ = geok.error_ellipse(pxp[0],pxp[1], sigyB , sigxB , sigxyB, scale= 10000)
                    
                    # c'est la fonction de Ghiliani, considérié comme étant celle fiable
                    # phi => angle à partir de X dans le sens trigo
                    # t   => angle à partir de Y dans le sens horaire
                    qxx = Ninv[ kkk + ipxp * 3    , kkk + ipxp * 3  ]
                    qyy = Ninv[ kkk + ipxp * 3 +1 , kkk + ipxp*3 +1 ]
                    qxy = Ninv[ kkk + ipxp * 3    , kkk + ipxp*3 +1 ] 
                    a,b,phi  = geok.error_ellipse_parameters( qxx , qyy , qxy , fuv , 0 )
                    a,b,t    = geok.error_ellipse_parameters( qxx , qyy , qxy , fuv , 1 )
                    #génération de l'ellipse de Ghiliani, on inverse le signe de phi car
                    # ellipse_get_coords compte dans le sens horaire
                    xee , yee = geok.ellipse_get_coords(a*10000,b*10000,pxp[0],pxp[1],phi)
        
                    # une autre fonction basée sur Stang and Borre
                    geok.error_ellipse_parameters_2(sigyB , sigxB , sigxyB)
                    geok.error_ellipse_parameters_2(sigxB , sigyB , sigxyB)
        
                    # Bref on illustre le fait que les ellipses de Stang and Borre
                    # sont vaseuses parce que la convention d'orientation n'est pas claire
                    # et on utilise à l'avenir les ellipses de Ghiliani
                    # DAUTANT que, il y avait en plus une modification perso 
                    # (signe - devant D et dX1) pour la fonction error_ellipse,
                    # à investiguer à l'occaz pourquoi S&B sont foireux 
        
                    #plt.plot(xee,yee,"o")
                    #plt.plot(xe1 ,ye1 ,"+")
                    #plt.plot(xe2 ,ye2 ,"x")

                    plt.plot(xee,yee,"-")
                    #plt.plot(xe1 ,ye1 ,"+")
                    #plt.plot(xe2 ,ye2 ,"x")
                else:
                    # le mode operationnel de Ghiliani
                    qxx = Ninv[ kkk + ipxp * 3    , kkk + ipxp * 3   ]
                    qyy = Ninv[ kkk + ipxp * 3 +1 , kkk + ipxp*3 +1  ]
                    qxy = Ninv[kkk + ipxp * 3,kkk + ipxp*3 +1] 
                    a,b,phi  = geok.error_ellipse_parameters( qxx , qyy , qxy , fuv )
                    print("qxx , qyy , qxy", qxx , qyy , qxy , phi)
                    print("a   ,b    , phi",a,b,phi)

            scale=1000
            xee , yee = geok.ellipse_get_coords(a*scale,b*scale,pxp[0],pxp[1],phi)
            plt.plot(xee,yee)
            plt.xlabel('East (m)')
            plt.ylabel('North (m)')
            

            plt.axis('equal')
            plt.savefig(protopath + "geo.png")
            plt.savefig(protopath + "geo.pdf")
        
        if plot_mode:            
            plt.close('all')

        print("INFO : STOP")
        print("INFO : END WITHOUT FAIL")

    except Exception as e:
        if force_mode:
            print("ERR : exception ... but we continue !")     
            print(e)
            F.stop()
            continue
        else:
            print("ERR : exception ... everything is stopped !")
            print(e)
            F.stop()
            raise
try:
    exp_stk = []
    for root, dirs, files in os.walk(exp_path):
        for name in files:
            if '.exp' in name:
                exp_stk.append(os.path.join(root, name))
                   
    out_sum_fil = acls.exp_summary(exp_stk,path_input=False,
                                   outdir=exp_path)
    print("INFO : SUMMARY FILE CREATED")
    print(out_sum_fil)

    #acls.plot_cntr_from_expdics(exp_path,exp_path,exp)
    #acls.plot_dist_from_expdics(exp_path,exp_path,exp)
except:
    print("WARN : fail of summary creation")
    pass
