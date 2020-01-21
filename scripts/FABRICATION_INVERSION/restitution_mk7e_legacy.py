# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 13:04:35 2015

@author: psakicki

v4 : en mode mono Z , on injecte le dZi comme observable
MK6v2 : implementation du sigma_batch_mode (160616)
MK6v3 : fusion avec GEODESEA
MK7 LEGACY : fork directly based on MK6v3 => NO CONFIG FILE
MK7b : for TPX1-GFZ computer
MK7d : implementation of the multi SSP option (with_epoch_SSP)
MK7e : encapsulation in a main function for clustering
"""

print("import start")
import acouclass as acls
import SSP as ssp
import matplotlib
#matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
from tempfile import mkdtemp
import numpy as np
import scipy
from scipy import interpolate
import scipy.optimize as optimize
import scipy.stats
import scipy.sparse as sprs
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
print("import middle")
import genefun as gf
import re
import copy
import sympy
import pandas as pd
import tabulate
import geodetik as geok
import genefun
import sys
import main_restit_fct 
import softs_runner
print("import end")
            
multiexp_mode      = 1
batch_mode         = 1
force_mode         = 1
force_mode_reading = 1
purge_mode         = 0
purge_prelim_mode  = 1
purge_sumfile_mode = 0
sigma_batch_mode   = 0
plot_mode          = 0
seed_batch_mode    = 0
ASM_from_files_mode = 0
GFZ_cluster_mode   = 0

if platform.node() in ('calipso','finlande','diamant','VirBox5HPHP'):
    gene_path         = "/home/psakicki/THESE/DATA/1506_GEODESEA/WORKING/OPERA_DATA"
    alternat_SSP_path = '/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/exemple/SSP/SSP_NOAA_dep5781_20050827000000'
    smart_SSP_Sdic_path = '/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/exemple/Sdic/Sdic_geodesea.pik'

    gene_path = '/home/psakicki/THESE/DATA/1506_GEODESEA/WORKING/OPERA2//'
    gene_path = '/home/psakicki/THESE/DATA/1506_GEODESEA/WORKING/OPERA3//'
    gene_path = '/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/working/OFFSETS_mk1'
    gene_path = '/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/working/TRAJECT_NOISE_mk1'
    gene_path = '/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/working/BUOYstyle_mk1'
    gene_path = '/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/working/REPRIZ'
    gene_path = '/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/working/'

elif platform.node() == 'psakicki-MS-16Y1':
    gene_path         = '/home/psakicki/Documents/CODES/acoustyx_toolbox_2/working'
    alternat_SSP_path = '/home/psakicki/Documents/CODES/acoustyx_toolbox_2/exemple/SSP/SSP_NOAA_dep5781_20050827000000'

elif platform.node() == 'diamant-old':
    gene_path = '/home/psakicki/GEODESEA_DATA/OPERA1'
    gene_path = '/home/psakicki/GEODESEA_DATA/OPERA2'
    gene_path = '/home/psakicki/GEODESEA_DATA/'

elif platform.node() == 'TPX1-GFZ':
    gene_path = "/home/psakicki/aaa_FOURBI/GNSSA_REBOOT1802_mini/"
    gene_path = "/home/psakicki/aaa_FOURBI/GNSSA_REBOOT1802/"
    gene_path = "/home/psakicki/CODES/acoustyx_toolbox_2_py3/working/"
    gene_path = "/home/psakicki/aaa_FOURBI/GNSSA_SIMU_regroupement/GNSSA_SSP_tempo_update_for_SFGM-IPGP18_beta/"
    gene_path = "/home/psakicki/aaa_FOURBI/GNSSA_SIMU_regroupement/GNSSA_TAMBOUILLE_IMPLEMENT_STYLE_SHOULD_BE_DEL/"
    gene_path = "/home/psakicki/aaa_FOURBI/GNSSA_SIMU_regroupement/SYNOPSIG/REBOOT/"
    gene_path = "/home/psakicki/aaa_FOURBI/GNSSA_SIMU_regroupement/SYNOPSIG/INTEGRAL_SSP/"

elif ("pgal" in platform.node()) or ("sigs" in platform.node()):
    gene_path = "/dsk/ggsp_pf/PLAYGROUND/psakicki/GNSSA/REBOOT1802_WORK/"
    gene_path = "/dsk/ggsp_pf/PLAYGROUND/psakicki/GNSSA/SSP_UPDATE_simple_mk1/"
    gene_path = "/dsk/ggsp_pf/PLAYGROUND/psakicki/GNSSA/SSP_UPDATE_simple_mk1b_cluster_implement_test/"

elif "kg3" in platform.node():
    gene_path = "/dsk/ggsp_pf/PLAYGROUND/psakicki/GNSSA/REBOOT1802_WORK_FAST/"
    gene_path = "/dsk/ggsp_pf/PLAYGROUND/psakicki/GNSSA/REBOOT1803_more_exps_WORK/"
    gene_path = "/dsk/ggsp_pf/PLAYGROUND/psakicki/scripts_PS/zyx_TOOLBOXs/acoustyx_toolbox_2_py3/working/"
    gene_path = "/dsk/ggsp_pf/PLAYGROUND/psakicki/GNSSA/Les_anciens_BATCH_RC_avec_les_modes_les_plus_recents/"
    gene_path = "/dsk/ggsp_pf/PLAYGROUND/psakicki/GNSSA/RB_EGU18/"

MNfile = 'P' # M pour les SIMU, N pour GEODESEA

expdirsuffix = ''
if platform.node() == 'calipso':
    exp  = 'IUEM_LAPTOP-3Beacontracking'
    exp  = 'test_repriz_adel_3x50_x2000_y2000_noisFalse-0__'
    exp  = "batc2_3x1000_x2000_y2000_nois1-1e-06_"
    exp  = 'batc2c_deriv_decalZ_1000_R50_noisTrue-1e-06__'
    exp  = 'batc2c_deriv_decalZ_1000_R50_noisTrue-1e-06'
    exp  = 'test_repriz_adel_5000_R50_noisTrue-1e-06__'
    exp  = 'test_repriz_adel_50_R50_noisTrue-1e-06__'
    exp  = 'test_repriz_adel_2500_R5000_nois1-1e-06_vit7.2_'
    exp  = 'test_repriz_FB1_2500_R5000_noisTrue-1e-06_vit7.2_'
    exp  = 'test_repriz_adel_FB0_5000_R50_noisFalse-0_vit500_'
    exp  = 'test_repriz_FB1_5000_R50_noisTrue-1e-06_vit7.2_'
    exp  = 'test_repriz_adel_FB0_mk3_500_R50_noisFalse-0_vit500_'
    exp  = 'test_repriz_FB1_mk6_5000_R50_noisTrue-1e-06_vit7.2_'
    exp  = 'test_repriz_FB1_mk6_5000_R50_noisTrue-1e-06_vit7.2_'
    exp  = 'IUEM_LAPTOP-3Beacontracking'
    exp  = 'test_repriz_adel_50_R50_noisTrue-1e-06__'
    exp  = 'test_repriz_adel_FB0_50_R50_noisFalse-0_vit500_'
    exp  = "test_repriz_adel_FB0_mk3_50_R50_noisFalse-0_vit500_"
    exp  = 'test_repriz_FB1_mk6_500_R50_noisTrue-1e-06_vit7.2_'
    exp  = 'test_repriz_FB1_mk6_ImpNoiz0_500_R50_nois1-1e-06__'
    exp  = 'test_repriz_FB1_mk6_ImpNoiz0_50_R50_nois1-0__'
    exp  = 'test_repriz_FB1_mk6_ImpNoiz0_500_R50_nois1-0__'
    exp  = 'test_repriz_FB1_mk6_ImpNoiz0_5000_R50_nois1-0__'
    exp  = 'ADEL22_50_R5000_nois1-0__'
    exp  = 'batc2c_deriv_decalZ_noping_notraj_noise_5000_R50_noisFalse-0__'
    exp  = "ADEL33b_5000_R50_noisFalse-0__"
    exp  = "ADEL66_5000_R50_nois1-1e-06___"
    exp  = 'TRAJECT_NOISE_mk1_1000_R50_nois1-0__traject_noise_0_0_0_'
    exp  = "ADEL999_1000_R50_nois1-0_vit1__"
    exp  = "ADEL999_5000_R50_nois1-0_vit1__"
    exp  =  "ADEL9994_5000_R1000_nois1-0_vit1__"
    exp  = "ADEL9993_5000_R50_nois1-0_vit1__"
    exp = "ADEL9994_500_R1000_nois1-0_vit1__"
    exp = 'ULTRAKLASSIK_mk4_20000_R1000_noisTrue-0___'
    exp = 'ULTRAKLASSIK_mk7_IMPROVEDNOISE_decimed_2000_R1000_noisTrue-1e-06___'
    expdirsuffix = ''
    exp  = "ADEL666_50_R1000_nois1-0_vit1__"
    exp  = "ADEL22_50_R1000_nois1-0_vit1__"
    exp  = 'batc2b_deriv_REBOOT_1500_R50_noisTrue-1e-06__'
    exp  = 'BATCH_RC_mk2_1000_R1000_nois1-1e-06___'
    expdirsuffix = 'restitMK6'
    exp  = 'BATCH_RC_mk2_1000_R1000_nois1-1e-06___'
    expdirsuffix = 'dZADEL'
    exp = 'ADOCdZ1_50_R100_nois1-1e-06_vit1__'
    exp = 'ADOCdZ1_3x100_x2000_y2000_nois1-1e-06_vit1__'
    expdirsuffix = ''
    exp='BATCH_RC_mk6_OPERASCENAR_tempor_5000_R50_nois1-0_vit1___'
    exp='BATCH_RC_mk5b_authentik_dekal_zonal_noise_wo_10koef_1000_R50_nois1-1e-06___'
    exp='BATCH_RC_mk6_OPERASCENAR_tempor_10_R50_nois1-0_vit1___'
    exp = 'IUEM_LAPTOP-3Beacontracking_Zcst'
    exp = 'IUEM_LAPTOP-3Beacontracking'

else:
    exp = 'IUEM_LAPTOP-3Beacontracking'
    exp = 'test_repriz_adel_5000_R50_noisTrue-1e-06__'
    exp = 'batc2b_deriv_REBOOT_50_R50_noisTrue-1e-06__seed_1126_v151014'
    exp = "FINALb_VariTempor_basic_1000_R50_nois1-0.0_vit1__date20070424_grand-carree_"


if not multiexp_mode:
    exp_lis = [exp]
else:
    # == en manuel
    exp_lis = ["batc2_3x100_x2000_y2000_nois1-1e-06_" ,
               "batc2_3x500_x2000_y2000_nois1-1e-06_" ,
               "batc2_3x100_x0_y0_nois1-1e-06_"       ,
               "batc2_3x500_x0_y0_nois1-1e-06_"       ,
               "batc2_3x1500_x0_y0_nois1-1e-06_"      ,
               "batc2_3x1500_x2000_y2000_nois1-1e-06_",
               "batc2_3x1000_x0_y0_nois1-1e-06_"      ,
               "batc2_3x1000_x2000_y2000_nois1-1e-06_"]

    exp_lis = ['batc2b_deriv_1500_R50_noisTrue-1e-06__',
               'batc2b_deriv_5000_R50_noisTrue-1e-06__']

    exp_lis = ['OPERA1/IUEM_LAPTOP-3Beacontracking',
               'OPERA2/IUEM_LAPTOP-3Beacontracking']

    exp_lis = ['IUEM_LAPTOP-3Beacontracking']

    # == avec une wildcard
    expprefix = 'tpo1'
    expprefix = 'batc2_*x0_y0*'
    expprefix = 'batc2_*'
    expprefix = 'batc2c_*'
    expprefix = 'TRAJECT_NOISE_*'
    expprefix = 'OFFSET*mk1b*'
    expprefix = 'TRAJECT_NOISE_*1b*'
    expprefix = 'BUOYstyle_mk*'
    expprefix = "STRAIGHT_TRAJ_OFFSET_mk1*"
    expprefix = "STRAIGHT_TRAJ_OFFSET_mk*"
    expprefix = "KLASSIK_DECAL_CENTER*"
    expprefix = 'ULTRAKLASSIK_mk3*'
    expprefix = "ADEL666_500_R50_nois1-0_vit1__*"
    expprefix = "ADEL777*3x1500*"
    expprefix = "ADEL777*3x3000*"
    expprefix = "ADEL777_3x100_x2000_y2000_nois1-0_vit1__*"
    expprefix = "ADEL777_50_R50_nois1-0_vit1__*"
    expprefix = "ADEL777_500_R50_nois1-0_vit1__*"
    expprefix = "ADEL777_3x100_x2000_y2000_nois1-0_vit1__*"
    expprefix = 'ADEL777_3x100_x200_y200_nois1-0_vit1__*'
    expprefix = "ADEL777_500_R5000_nois1-0_vit1__*"
    expprefix = "ADEL777_500_R50_nois1-0_vit1_*"
    expprefix = "ADEL777_5000_R50_nois1-0_vit1__*"
    expprefix = "ADEL888_*"
    expprefix = "batc2b_deriv_1500_R50_noisTrue-1e-06__*"
    expprefix = "ADEL999_3x100_x2000_y2000_nois1-0_vit1__*"
    expprefix = "ADEL999_3x100_x2000_y7000_nois1-0_vit1__*"
    expprefix = "ULTRAKLASSIK_mk4*"
    expprefix = "ULTRAKLASSIK_mk5_noIMPROVEDNOISE*"
    expprefix = "ULTRAKLASSIK_mk6_IMPROVEDNOISE_fullrandomized*"
    expprefix = 'ULTRAKLASSIK_mk7_IMPROVEDNOISE_decimed*10000*'
    expprefix = 'ULTRAKLASSIK_mk8*'
    expprefix = 'BATCH_RC_mk2*'
    expprefix = 'ULTRAKLASSIK_mk9*'
    expprefix = 'BATCH_RC_mk3*'
    expprefix = "BATCH_RC_mk4*"
    expprefix = "BATCH_RC_mk5b_authentik_dekal_zonal_noise_wo_10koef_5000_R500_nois1-1e-06___*"
    expprefix = "BATCH_RC_mk5*"
    expprefix = "BATCH_RC_mk*"
    expprefix = "*OPERASCENAR*pseudograd*10p6*"
    expprefix = "GRAD_mk2d_50000_R50_nois1-0____*"
    expprefix = 'SSFssrealik_1*'
    expprefix = "SSFssrealik_1_1000_R50_nois1-0____*"
    expprefix = "SSFrealik_2*1000*_"
    expprefix = "SSFrealik_3_10p5_3500*3x3333*"
    expprefix = "*smallarray*"
    expprefix = "*proxyMOVE*"
    expprefix = "*METEOR_pseudoSSF_1_3x3333_x500_y500_ang0_nois1-0____*"
    expprefix = "*ATAL_pseudoSSF_1*"
    expprefix = "*500_pseudoSSF_1*"
    expprefix = "*test_triangle*"
    expprefix = "*batc2c*"
    expprefix = "*FINALb*"
    expprefix = '*FINALb_VariTempor_zonal1045*'
    expprefix = "*FINALb_BUOYstyle*mk6_*"
    expprefix = "*FINALb*basic*"
    expprefix = "*FINAL*"
    expprefix = "*FINALb_BUOYstyle*mk6*"
    expprefix = "*FINAL*"
    expprefix = "*VariTemporCORIGEDb*"
    expprefix = "1709C*"

    expprefix = "*BATCH_RC_mk3*R10_*"
    expprefix = "*BATCH_RC_mk3_authentik_dekal_3x333_x2000_y2000_nois1-1e-06___*"
    expprefix = "*BATCH_RC_mk3*"

    expprefix = "*FINAL*"
    expprefix = "*FINALb_VariTemporCORIGEDb_basic_REBOOT1803*"

    expprefix = "RB_EGU18*R50*"
    expprefix = "RB_EGU18*x33*"
    expprefix = "RB_EGU18*201002*"
    expprefix = "RB_EGU18*_R1_*"
    expprefix = "RB_EGU18_TEMPO-D*"
    expprefix = "RB_EGU18_TEMPO-E*_3x333*"

    expprefix = "*RB_EGU18*0_R*"
    expprefix = "*RB_EGU18*3x33*"

    expprefix = "*RB_EGU18_TEMPO-A_500_R50_nois1-0.0_vit1__date20100225_grand-carree_Kdrift10_*"

    expprefix = "*BATCH_RC_mk3*"

    expprefix = "*RB_EGU18_TEMPO-A*0_R*"
    expprefix = "*RB_EGU18_TEMPO-A*3x16*Kdrift10_*"
    
    expprefix = "*RB_EGU18_TEMPO-A*"


    exp_lis   = glob.glob(os.path.join(gene_path,expprefix + '*'))
    exp_lis   = list(reversed(sorted([os.path.basename(e) for e in exp_lis])))



    purge_ext_list = ("/*.log","/*.exp","/*.pdf","/*.png","/*.V")

    if purge_sumfile_mode:
        print("WARN : sum file(s) will be purged too !")
        purge_ext_list = purge_ext_list + ("/*.sum",)

    if purge_prelim_mode:
        print("WARN : purge PRELIMINARY mode !!!")
        print(gene_path)
        print("10sec for cancel")
        time.sleep(10)

        for exp in exp_lis:
            for ext in purge_ext_list:
                rmlis = glob.glob(os.path.join(gene_path,exp) + ext)
                if len(rmlis) == 0:
                    print('WARN : purge prelim. : remove list empty for ' + exp + ext)
                for f in rmlis:
                    print('INFO : remove : ' , f)
                    os.remove(f)

for iexp , exp in enumerate(exp_lis):

    print("##### EXPERIENCE : " , iexp + 1 , '/' , len(exp_lis))

    expini   = exp
    exp_path = os.path.join(gene_path,expini + expdirsuffix)
    exp      = os.path.basename(expini)
    try:
        bigdico  = acls.give_me_the_path(exp_path,exp)#,[1,3,5])
    except Exception as e:
        if force_mode_reading:
            print(e)
            print("WARN : error while reading but force_mode_reading is ON")
            continue
        else:
            raise e

    if 'OPERA' in exp_path and 'GEODESEA' in exp_path:
        MNfile = 'O'
    elif 0:
        MNfile = 'M'
        MNfile = 'P'
    else:
        MNfile = acls.auto_find_extension_MNOP_files(exp_path)

    if purge_mode and not purge_prelim_mode:
        print("WARN : purge mode !!!")
        print("5sec for cancel")
        time.sleep(8)

        for ext in purge_ext_list:
            rmlis = glob.glob(exp_path + ext)
            if len(rmlis) == 0:
                print('WARN : purge prelim. : remove list empty for ' + exp + ext)
            for f in rmlis:
                print('INFO : remove : ' , f)
                os.remove(f)

    if batch_mode:
        bool_kw_lis = [ 'with_alternat_SSP' ,
                        'with_monoZ'        ,
                        'with_barycenter'   ]

        # bool_kw_lis if you want to apply a Jack Knife
        if 0:
            bool_kw_lis = [ 'with_V_4_P_reinject'  ,
                            'with_jackknife'       ,
                            'with_invert_jk'       ]

            booldic_var = genefun.boolean_dict(bool_kw_lis)
            booldic_lis = geok.kwargs_for_jacobian({'with_ssp_bilin'      : 0,
                                                    'with_decimate_SSP'   : 0,
                                                    'with_alternat_SSP'   : 0,
                                                    'with_barycenter'     : 0,
                                                    'with_monoZ'          : 0,
                                                    'with_BL'             : 0,
                                                    'with_dzcst'          : 0,
                                                    'with_zmaster'        : 0,
                                                    'with_time_window'    : 0,
                                                    'with_specific_timwin_bool' : 0},
                                                    booldic_var)
        if 0:
            # bool_kw_lis if you want to apply a Time WIndows
            bool_kw_lis = [ 'with_V_4_P_reinject'  ,
            'with_time_window' ,
            'with_specific_timwin_bool' ]

            booldic_var = genefun.boolean_dict(bool_kw_lis)
            booldic_lis = geok.kwargs_for_jacobian({'with_ssp_bilin'    : 0,
                                                    'with_decimate_SSP' : 0,
                                                    'with_alternat_SSP' : 0,
                                                    'with_barycenter'   : 0,
                                                    'with_monoZ'        : 0,
                                                    'with_BL'           : 0,
                                                    'with_dzcst'        : 0,
                                                    'with_zmaster'      : 0,
                                                    'with_jackknife'    : 0,
                                                    'with_invert_jk'    : 0},
                                                    booldic_var)

        if 0:
            bool_kw_lis = ['with_alternat_SSP' ,
            'with_monoZ',
            'with_barycenter',
            'with_BL' ]

            booldic_var = genefun.boolean_dict(bool_kw_lis)

            booldic_lis = geok.kwargs_for_jacobian({'with_ssp_bilin'    : 0,
                                                   'with_decimate_SSP'  : 0,
                                                   'with_dzcst'         : 0,
                                                   'with_zmaster'       : 0,
                                                   'with_jackknife'     : 0,
                                                   'with_invert_jk'     : 0,
                                                   'with_FB_mode'       : 1,
                                                   'with_time_window'   : 0,
                                                   'with_V_4_P_reinject' : 0},
                                                   booldic_var)

        if 0:
            # BLOC pour retrouver a quoi correspond monoZ , Z master , dz_cst
            bool_kw_lis = ['with_alternat_SSP' ,
            'with_monoZ',
            'with_barycenter',
            'with_BL'  ,
            'with_dzcst',
            'with_zmaster']

            booldic_var = genefun.boolean_dict(bool_kw_lis)

            booldic_lis = geok.kwargs_for_jacobian({'with_ssp_bilin'     : 0,
                                                   'with_decimate_SSP'   : 0,
                                                   'with_zmaster'        : 0,
                                                   'with_jackknife'      : 0,
                                                   'with_invert_jk'      : 0,
                                                   'with_FB_mode'        : 1,
                                                   'with_time_window'    : 0,
                                                   'with_V_4_P_reinject' : 0},
                                                   booldic_var)
        if 0:
            # BLOC du barycentre
            bool_kw_lis = ['with_barycenter']

            booldic_var = genefun.boolean_dict(bool_kw_lis)

            booldic_lis = geok.kwargs_for_jacobian({'with_ssp_bilin'     : 0,
                                                   'with_decimate_SSP'   : 0,
                                                   'with_zmaster'        : 0,
                                                   'with_jackknife'      : 0,
                                                   'with_invert_jk'      : 0,
                                                   'with_FB_mode'        : 1,
                                                   'with_time_window'    : 0,
                                                   'with_V_4_P_reinject' : 0,
                                                     'with_alternat_SSP' : 0,
                                                     'with_monoZ'        : 0,
                                                     'with_barycenter'   : 0,
                                                     'with_BL'           : 0,
                                                     'with_dzcst'        : 0,
                                                     'with_zmaster'      : 0},
                                                   booldic_var)


        if 0:
            # BLOC pour un BATCH STANDARD
            bool_kw_lis = ['with_monoZ',
            'with_barycenter',
            'with_BL'  ,
            'with_dzcst',
            'with_zmaster']

            booldic_var = genefun.boolean_dict(bool_kw_lis)

            booldic_lis = geok.kwargs_for_jacobian({'with_ssp_bilin'     : 0,
                                                   'with_decimate_SSP'   : 0,
                                                   'with_zmaster'        : 0,
                                                   'with_jackknife'      : 0,
                                                   'with_invert_jk'      : 0,
                                                   'with_FB_mode'        : 1,
                                                   'with_time_window'    : 0,
                                                   'with_V_4_P_reinject' : 0,
                                                   'with_alternat_SSP'   : 0},
                                                   booldic_var)

        if 0:
            # BLOC des BASELINES
            bool_kw_lis = ['with_BL']

            booldic_var = genefun.boolean_dict(bool_kw_lis)

            booldic_lis = geok.kwargs_for_jacobian({'with_ssp_bilin'     : 0,
                                                   'with_decimate_SSP'   : 0,
                                                   'with_zmaster'        : 0,
                                                   'with_jackknife'      : 0,
                                                   'with_invert_jk'      : 0,
                                                   'with_FB_mode'        : 1,
                                                   'with_time_window'    : 0,
                                                   'with_V_4_P_reinject' : 0,
                                                     'with_alternat_SSP' : 0,
                                                     'with_monoZ'        : 0,
                                                     'with_barycenter'   : 0,
                                                     'with_barycenter'   : 0,
                                                     'with_dzcst'        : 0,
                                                     'with_zmaster'      : 0},
                                                   booldic_var)



        if 0:
            # BLOC pour REBOOT1802
            # basé sur le BATCH STANDARD, NE PAS MODIFIER
            bool_kw_lis = ['with_monoZ',
            'with_barycenter',
            'with_dzcst',
            'with_zmaster']

            booldic_var = genefun.boolean_dict(bool_kw_lis)

            booldic_lis = geok.kwargs_for_jacobian({'with_ssp_bilin'     : 0,
                                                   'with_decimate_SSP'   : 0,
                                                   'with_zmaster'        : 0,
                                                   'with_jackknife'      : 0,
                                                   'with_invert_jk'      : 0,
                                                   'with_FB_mode'        : 1,
                                                   'with_time_window'    : 0,
                                                   'with_V_4_P_reinject' : 0,
                                                   'with_alternat_SSP'   : 0,
                                                   'with_BL'             : 1},
                                                   booldic_var)


        if 0:
            # BLOC pour REBOOT1802
            # MODIFIABLE POUR ESSAIS
            # Ce bloc est pas mal, on le garde
            bool_kw_lis = ['with_barycenter',
                           'with_BL',
                           'with_zmaster']

            booldic_var = genefun.boolean_dict(bool_kw_lis)

            booldic_lis = geok.kwargs_for_jacobian({'with_ssp_bilin'     : 0,
                                                   'with_decimate_SSP'   : 1,
                                                   'with_jackknife'      : 0,
                                                   'with_invert_jk'      : 0,
                                                   'with_FB_mode'        : 1,
                                                   'with_time_window'    : 0,
                                                   'with_V_4_P_reinject' : 0,
                                                   'with_alternat_SSP'   : 0,
                                                   'with_monoZ'          : 0,
                                                   'with_dzcst'          : 0},
                                                   booldic_var)

        if 1:
            # BLOC pour REBOOT1802
            # MODIFIABLE POUR ESSAIS EQUIVALENT A DU NON BATCH

            booldic_var = genefun.boolean_dict(bool_kw_lis)

            booldic_lis = geok.kwargs_for_jacobian({'with_ssp_bilin'     : 0,
                                                   'with_decimate_SSP'   : 1,
                                                   'with_jackknife'      : 0,
                                                   'with_invert_jk'      : 0,
                                                   'with_FB_mode'        : 1,
                                                   'with_time_window'    : 0,
                                                   'with_V_4_P_reinject' : 0,
                                                   'with_alternat_SSP'   : 0,
                                                   'with_monoZ'          : 0,
                                                   'with_dzcst'          : 0,
                                                   'with_barycenter'     : 1,
                                                   'with_BL'             : 1,
                                                   'with_zmaster'        : 1,
                                                   "with_epoch_SSP"      : 0,
                                                   "with_smart_SSP"      : 0,
                                                   "with_integrate_SSP"  : 1,
                                                   },
                                                   dict())



        if 0:
            # on inaugure le BLOC SEMI BATCH

            basis_dic = {  'with_ssp_bilin'      : 0,
                           'with_decimate_SSP'   : 0,
                           'with_zmaster'        : 0,
                           'with_jackknife'      : 0,
                           'with_invert_jk'      : 0,
                           'with_FB_mode'        : 1,
                           'with_time_window'    : 0,
                           'with_V_4_P_reinject' : 0,
                           'with_alternat_SSP'   : 0}

            bool_kw_lis = ['with_monoZ',
            'with_barycenter',
            'with_BL'   ,
            'with_dzcst',
            'with_zmaster']

            vari_dic_lis =[{'with_monoZ'     : 0,
                           'with_barycenter' : 1,
                           'with_BL'         : 1,
                           'with_dzcst'      : 0,
                           'with_zmaster'    : 1 } ,
                          {'with_monoZ'      : 1,
                           'with_barycenter' : 1,
                           'with_BL'         : 1,
                           'with_dzcst'      : 1,
                           'with_zmaster'    : 0 } ,
                          {'with_monoZ'      : 1,
                           'with_barycenter' : 1,
                           'with_BL'         : 0,
                           'with_dzcst'      : 1,
                           'with_zmaster'    : 0 } ,
                          {'with_monoZ'      : 0,
                           'with_barycenter' : 0,
                           'with_BL'         : 0,
                           'with_dzcst'      : 0,
                           'with_zmaster'    : 0 } ]
            # Basé sur
            #BATCH_RC_mk3_authentik_dekal_1000_R10_nois1-1e-06____01101_20160613_015352.exp 0.00268858107561
            #BATCH_RC_mk3_authentik_dekal_1000_R10_nois1-1e-06____11110_20160613_011205.exp 0.00131287724058
            #BATCH_RC_mk3_authentik_dekal_1000_R10_nois1-1e-06____11010_20160613_011724.exp 0.458497986848
            #BATCH_RC_mk3_authentik_dekal_1000_R10_nois1-1e-06____00000_20160613_020902.exp 0.86706021259

            #proto_vari_dic     = {'with_monoZ'     : [0,1,1,0],
            #                     'with_barycenter' : [1,1,1,0],
            #                     'with_BL'         : [1,1,0,0],
            #                     'with_dzcst'      : [0,1,1,0],
            #                     'with_zmaster'    : [1,0,0,0]}
            #
            #vari_dic_lis = []
            #for iii  in range(len(proto_vari_dic[proto_vari_dic.keys()[0]])):
            #    vari_dic_lis.append(collections.OrderedDict())
            #    for k in proto_vari_dic.keys():
            #        vari_dic_lis[-1][k] = proto_vari_dic[k][iii]



            #for k,iitt in proto_vari_dic_lis.items():
            #    vari_dic_lis.append(dict())
            #    for it in iitt:
            #        vari_dic_lis[-1][k] = it

            booldic_lis = [ gf.dicts_merge(basis_dic,d) for d in vari_dic_lis ]


        if 0:
            # BLOC SEMI BATCH pour L'experience principale (zmaster/BLs)
            # EQUIVALENT A DU NON BATCH ...

            basis_dic = {  'with_ssp_bilin'      : 0,
                           'with_decimate_SSP'   : 0,
                           'with_zmaster'        : 0,
                           'with_jackknife'      : 0,
                           'with_invert_jk'      : 0,
                           'with_FB_mode'        : 1,
                           'with_time_window'    : 0,
                           'with_V_4_P_reinject' : 0,
                           'with_alternat_SSP'   : 0}

            bool_kw_lis = ['with_monoZ',
            'with_barycenter',
            'with_BL'   ,
            'with_dzcst',
            'with_zmaster']

            vari_dic_lis =[{'with_monoZ'     : 0,
                           'with_barycenter' : 1,
                           'with_BL'         : 1,
                           'with_dzcst'      : 0,
                           'with_zmaster'    : 1 }]

            booldic_lis = [ gf.dicts_merge(basis_dic,d) for d in vari_dic_lis ]

        if 0:
            # BLOC SEMI BATCH pour L'experience annexe basée sur les dz as constant (zmaster/BLs)
            # Que l'on considère sauf mention ctraire comme le meilleur mode au 08 sept 2016
            # EQUIVALENT A DU NON BATCH ...

            basis_dic = {  'with_ssp_bilin'      : 0,
                           'with_decimate_SSP'   : 1,
                           'with_zmaster'        : 0,
                           'with_jackknife'      : 0,
                           'with_invert_jk'      : 0,
                           'with_FB_mode'        : 1,
                           'with_time_window'    : 0,
                           'with_V_4_P_reinject' : 0,
                           'with_alternat_SSP'   : 0,
                           'with_specific_timwin_bool' : 0}

            bool_kw_lis = ['with_monoZ',
            'with_barycenter',
            'with_BL'   ,
            'with_dzcst',
            'with_zmaster']

            vari_dic_lis =[{'with_monoZ'     : 1,
                           'with_barycenter' : 0,
                           'with_BL'         : 1,
                           'with_dzcst'      : 1,
                           'with_zmaster'    : 0 }]

            booldic_lis = [ gf.dicts_merge(basis_dic,d) for d in vari_dic_lis ]

        if 0:
            # BLOC SEMI BATCH pour le reboot du 1709 on veut juste une restitution sans prise en cmpt du Z
            # EQUIVALENT A DU NON BATCH ...

            basis_dic = {  'with_ssp_bilin'      : 0,
                           'with_decimate_SSP'   : 1,
                           'with_zmaster'        : 0,
                           'with_jackknife'      : 0,
                           'with_invert_jk'      : 0,
                           'with_FB_mode'        : 1,
                           'with_time_window'    : 0,
                           'with_V_4_P_reinject' : 0,
                           'with_alternat_SSP'   : 0,
                           'with_specific_timwin_bool' : 0}


            vari_dic_lis =[{'with_monoZ'     : 0,
                           'with_barycenter' : 1,
                           'with_BL'         : 0,
                           'with_dzcst'      : 0,
                           'with_zmaster'    : 0 }]

            booldic_lis = [ gf.dicts_merge(basis_dic,d) for d in vari_dic_lis ]

        if 0:
            # BLOC SEMI BATCH pour le REBOOT1802
            # basée sur :
            # BLOC SEMI BATCH pour L'experience annexe basée sur les dz as constant (zmaster/BLs)
            # Que l'on considère sauf mention ctraire comme le meilleur mode au 08 sept 2016

            basis_dic = {  'with_ssp_bilin'      : 0,
                           'with_decimate_SSP'   : 1,
                           'with_zmaster'        : 0,
                           'with_jackknife'      : 0,
                           'with_invert_jk'      : 0,
                           'with_FB_mode'        : 1,
                           'with_time_window'    : 0,
                           'with_V_4_P_reinject' : 0,
                           'with_alternat_SSP'   : 0}

            bool_kw_lis = ['with_monoZ',
            'with_barycenter',
            'with_BL'   ,
            'with_dzcst',
            'with_zmaster']

            vari_dic_lis =[{'with_monoZ'     : 0,
                           'with_barycenter' : 1,
                           'with_BL'         : 1,
                           'with_dzcst'      : 0,
                           'with_zmaster'    : 1 }]

            booldic_lis = [ gf.dicts_merge(basis_dic,d) for d in vari_dic_lis ]


        if 0:
            # BLOC SEMI BATCH pour le REBOOT1802
            # basée sur :
            # BLOC SEMI BATCH pour L'experience annexe basée sur les dz as constant (zmaster/BLs)
            # Experimental on fait ce que l'on veux avec les options

            basis_dic = {  'with_ssp_bilin'      : 0,
                           'with_decimate_SSP'   : 1,
                           'with_zmaster'        : 0,
                           'with_jackknife'      : 0,
                           'with_invert_jk'      : 0,
                           'with_FB_mode'        : 1,
                           'with_time_window'    : 0,
                           'with_V_4_P_reinject' : 0,
                           'with_alternat_SSP'   : 0}

            bool_kw_lis = ['with_monoZ',
            'with_barycenter',
            'with_BL'   ,
            'with_dzcst',
            'with_zmaster']

            vari_dic_lis =[{'with_monoZ'     : 0,
                           'with_barycenter' : 0,
                           'with_BL'         : 1,
                           'with_dzcst'      : 0,
                           'with_zmaster'    : 1}]

            booldic_lis = [ gf.dicts_merge(basis_dic,d) for d in vari_dic_lis ]

    else:
        booldic_lis = [dict()]


    # ================================================
    # et vint le sigma_batch_mode 160616
    if sigma_batch_mode:

        sigma_defo_ASM     = 10**-5
        sigma_defo_BL      = 10**-3
        sigma_defo_zmaster = 10**-8

        sigmadic_var = {'sigma_pxp_apri'      : [10,5,1],
                        'sigma_dZ_apri'       : [0,10**-3,10**-1,1],
                        'sigma_BL_apri'       : [0,10**-3,10**-1,1],
                         'sigma_defo_ASM'     : 10**-5,
                         'sigma_defo_BL'      : 10**-3,
                         'sigma_defo_zmaster' : 10**-8}

        sigmadic_var = {'sigma_pxp_apri'     : [10]    ,
                        'sigma_dZ_apri'      : [10**-2],
                        'sigma_BL_apri'      : [10**-2],
                        'sigma_defo_ASM'     : np.power(10.,[-6,-5,-4,-3,-2,-1,0]),
                        'sigma_defo_BL'      : np.power(10.,[-6,-5,-4,-3,-2,-1,0]),
                        'sigma_defo_zmaster' : np.power(10.,[-6,-5,-4,-3,-2,-1,0])}

        # selection meilleur ED pour 1000ping R1000
        sigmadic_var = {'sigma_pxp_apri'     : [10]    ,
                        'sigma_dZ_apri'      : [10**-2],
                        'sigma_BL_apri'      : [10**-2],
                        'sigma_defo_ASM'     : [1.000000],
                        'sigma_defo_BL'      : [0.000100],
                        'sigma_defo_zmaster' : [0.010000]}
        # selection des 4 meilleurs E2D
        sigmadic_var = {'sigma_pxp_apri'     : [10]    ,
                        'sigma_dZ_apri'      : [10**-2],
                        'sigma_BL_apri'      : [10**-2],
                        'sigma_defo_ASM'     : [1.0, 1e-06 , 0.001, 0.01, 0.1],
                        'sigma_defo_BL'      : [1.0, 1e-06 , 1.0e-05, 0.0001, 0.001],
                        'sigma_defo_zmaster' : [1.0e-05, 0.001, 0.0001, 0.1, 0.01]}

        # selection meilleur ED pour 1000ping R10
        sigmadic_var = {'sigma_pxp_apri'     : [10]    ,
                        'sigma_dZ_apri'      : [10**-2],
                        'sigma_BL_apri'      : [10**-2],
                        'sigma_defo_ASM'     : [0.100000],
                        'sigma_defo_BL'      : [0.000100],
                        'sigma_defo_zmaster' : [0.010000]}

        #  meilleur ED pour 1000ping R10 MAIS avec sig z fortement contraint
        sigmadic_var = {'sigma_pxp_apri'     : [10]    ,
                        'sigma_dZ_apri'      : [10**-2],
                        'sigma_BL_apri'      : [10**-2],
                        'sigma_defo_ASM'     : [0.100000],
                        'sigma_defo_BL'      : [0.000100],
                        'sigma_defo_zmaster' : [10**-4,10**-5]}


        # LEGACY bruits et poids lorsque la sum dX bary est null
        sigmadic_var = {'sigma_pxp_apri'     : [1]     ,
                        'sigma_dZ_apri'      : [10**-3],
                        'sigma_BL_apri'      : [10**-2],
                        'sigma_defo_ASM'     : [10**-5],
                        'sigma_defo_BL'      : [10**-3],
                        'sigma_defo_zmaster' : [10**-3]}


        # LEGACY FOLLOWING ON, on cherche le paramètre qui fait que
        # la sum dX bary est pas null
        ### jusque là ca marche
        sigmadic_var = {'sigma_pxp_apri'     : [10]    ,
                        'sigma_dZ_apri'      : [10**-2],
                        'sigma_BL_apri'      : [10**-2],
                        'sigma_defo_ASM'     : [10**-5],
                        'sigma_defo_BL'      : [10**-3],
                        'sigma_defo_zmaster' : [10**-5]}

        # conclusion si le poids des BLs est trop fort, alors on perd
        # la condtion de nullité
        # ces paramètres empiriques vont bien
        # En fait non, les dZ sont ultra préponderant ...
        # et donne le même résultat même avec 10pings
        sigmadic_var = {'sigma_pxp_apri'     : [10]    ,
                        'sigma_dZ_apri'      : [10**-2],
                        'sigma_BL_apri'      : [10**-2],
                        'sigma_defo_ASM'     : [10**-2],
                        'sigma_defo_BL'      : [10**-3],
                        'sigma_defo_zmaster' : [10**-5]}

        # On teste ceux là pour la forme (on chage le poids ASM)
        # ca ca marche bien mais c'est un peu bof bof avec le
        # dZ en terme de précision
        sigmadic_var = {'sigma_pxp_apri'     : [10]    ,
                        'sigma_dZ_apri'      : [10**-0],
                        'sigma_BL_apri'      : [10**-2],
                        'sigma_defo_ASM'     : [10**-5],
                        'sigma_defo_BL'      : [10**-3],
                        'sigma_defo_zmaster' : [10**-0]}


        # On teste ceux là pour la forme (on chage le poids ASM)
        sigmadic_var = {'sigma_pxp_apri'     : [10]    ,
                        'sigma_dZ_apri'      : [10**-2],
                        'sigma_BL_apri'      : [10**-2],
                        'sigma_defo_ASM'     : [10**-0,10**-9],
                        'sigma_defo_BL'      : [10**-3],
                        'sigma_defo_zmaster' : [10**-3]}

        sigmadic_var = {'sigma_pxp_apri'     : [10]    ,
                        'sigma_dZ_apri'      : [10**-12],
                        'sigma_BL_apri'      : [10**-12],
                        'sigma_defo_ASM'     : [10**-6],
                        'sigma_defo_BL'      : [10**-3],
                        'sigma_defo_zmaster' : [10**-3]}

        # pour un nieme test sur les sigmas pour les sigma finaux (161004)
        sigmadic_var = {'sigma_pxp_apri'     : [10]    ,
                        'sigma_dZ_apri'      : [10**-0],
                        'sigma_BL_apri'      : [10**-0],
                        'sigma_defo_ASM'     : [10**-16],
                        'sigma_defo_BL'      : [10**-3],
                        'sigma_defo_zmaster' : [10**-3]}


        sigmadic_var = {'sigma_pxp_apri'     : [10]    ,
                        'sigma_dZ_apri'      : [10**-12],
                        'sigma_BL_apri'      : [10**-12],
                        'sigma_defo_ASM'     : [10**-6],
                        'sigma_defo_BL'      : [10**-3],
                        'sigma_defo_zmaster' : [10**-3]}

        # Trouver des sigma cohérents ...
        sigmadic_var = {'sigma_pxp_apri'     : [10]    ,
                        'sigma_dZ_apri'      : [10**-12],
                        'sigma_BL_apri'      : [10**-12],
                        'sigma_defo_ASM'     : [.5],
                        'sigma_defo_BL'      : [10**-3],
                        'sigma_defo_zmaster' : [10**-3]}

        # REBOOT1802 : basé sur les logs des tests BATCH_RC_mk3
        # correspondants aux tests de combinaisions du manuscrit
        sigmadic_var = {'sigma_pxp_apri'     : [1]     ,
                        'sigma_dZ_apri'      : [10**-3],
                        'sigma_BL_apri'      : [10**-2],
                        'sigma_defo_ASM'     : [10**-5],
                        'sigma_defo_BL'      : [10**-3],
                        'sigma_defo_zmaster' : [10**-2]}

        # REBOOT1802 : basé sur les logs des tests BATCH_RC_mk3
        # ici on blinde le sigma des dZ
        sigmadic_var = {'sigma_pxp_apri'     : [1]     ,
                        'sigma_dZ_apri'      : [10**-3],
                        'sigma_BL_apri'      : [10**-2],
                        'sigma_defo_ASM'     : [10**-5],
                        'sigma_defo_BL'      : [10**-3],
                        'sigma_defo_zmaster' : [10**-3,10**-2,5*10**-3]}

        # REBOOT1802 : basé sur les logs des tests BATCH_RC_mk3
        # avec sigma_defo_zmaster' : 10**-3 parce que ca semble etre le mieux
        sigmadic_var = {'sigma_pxp_apri'     : [1]      ,
                        'sigma_dZ_apri'      : [10**-2] , ### -2 TESTé et validé dans le REBOTT1802 Run9
                        #'sigma_BL_apri'      : [10**-2] , ### -2 TESTé et validé dans le REBOTT1802 Run8
                        'sigma_BL_apri'      : [10**-3] , ### -2 TESTé et validé dans le REBOTT1802 Run8
                        'sigma_defo_ASM'     : [10**-5] , # -5
                        #'sigma_defo_BL'      : [10**-3,.5*10**-4,10**-4,10**-5] , ### -3 initial
                        #'sigma_defo_BL'      : [10**-3] , ### -3 TESTé et validé dans le REBOTT1802 Run8
                        'sigma_defo_BL'      : [10**-4],
                        'sigma_defo_zmaster' : [10**-3]}  ### -3 TESTé et validé dans le REBOTT1802 Run8

                        # le bruit sur l'observation joue peu, c'est le poids qui joue, ccl REBOTT1802 Run8&9

        sigmadic_lis = geok.kwargs_for_jacobian(dict(),sigmadic_var)

    else:
        sigmadic_lis = [dict()]

    # End of sigma batch mode

    if seed_batch_mode:
        #### Les seeds sur la trajectoire
        seeddic_var = { 'k_z_apri'     : list(np.arange(10,12))  ,
                        'k_dZ_apri'    : list(np.arange(20,22))  ,
                        'k_BL_apri'    : list(np.arange(30,32))  }
        seeddic_lis = geok.kwargs_for_jacobian(dict(),seeddic_var)
        

        #### Les seeds sur la trajectoire en mode period 
        seeddic_var = { 'seed_epoch_SSP_period_parameter'     : [3600/2,3600,7200] , #Time in second
                        "seed_epoch_SSP_mode" : ["period_first",
                                                 "period_last",
                                                 "period_mean",
                                                 "period_middle"]}

        #### Les seeds sur la trajectoire en mode chunck 
        seeddic_var = { 'seed_epoch_SSP_chunk_parameter'     : [1,5,10] ,
                        "seed_epoch_SSP_mode" : ["chunk_first",
                                                 "chunk_last",
                                                 "chunk_mean",
                                                 "chunk_middle",
                                                 "none"]}
                        
                        
        seeddic_lis = geok.kwargs_for_jacobian(dict(),seeddic_var)

    else:
        seeddic_lis = [dict()]

    ####################################################
    ###############  BIG LOOP ##########################
    ####################################################  
    
    
    kommand_stk = []
    for booldic , sigmadic ,seeddic    in itertools.product(booldic_lis,
                                                            sigmadic_lis,
                                                            seeddic_lis):
        
        
        ##### PREPARE A PICKLE WITH THE VALUES
        global_dic = copy.copy(dict(globals()))
        outdir = "/home/psakicki/tmp_pickles/"
        outname_pickle_tmp = "GNSSA_param_tmp"
        global_dic_wo_fcts = dict()        
        for k,v in global_dic.items():
            try:
                genefun.pickle_saver(v,outdir,outname_pickle_tmp,
                                     timestamp=False)
                global_dic_wo_fcts[k] = v
            except:
                continue
            
        GENERAL_dic = global_dic
            
        if not GFZ_cluster_mode:
            main_restit_fct.main_GNSSA_restit_fct(GENERAL_dic)
            
        ##### HERE ARE LOADED THE COMMAND IN A LIST
        else:
            print("############### Cluster Mode")
            outname_pickle = "GNSSA_param_" + str(np.random.randint(10000)).zfill(6)
            pickle_path = genefun.pickle_saver(global_dic_wo_fcts,outdir,
                                               outname_pickle,
                                               timestamp=True)
            print("###### Pickle of variables :")
            print(pickle_path)    
            print(platform.node())
            #if ("sigs" in platform.node()) or ("pgal1" == platform.node()):
            if ("sigs" in platform.node()) or ("pgal" in platform.node()):
                kommand_prefix = "/dsk/igs2/soft_wrk/psakicki/SOFT_EPOS8_BIN_TOOLS/SCRIPTS/cluster_job.pl -c 'pyps" 
                function_path  = "/dsk/ggsp_pf/PLAYGROUND/psakicki/scripts_PS/zyx_TOOLBOXs/acoustyx_toolbox_2_py3/lib/main_restit_fct.py"
                kommand_sufix  = "'"
            else:
                kommand_prefix = "python3"
                function_path  = "/home/psakicki/CODES/acoustyx_toolbox_2_py3/lib/main_restit_fct.py"
                kommand_sufix  = ""
                
            kommand = " ".join((kommand_prefix,function_path,pickle_path,kommand_sufix))
            kommand_stk.append(kommand)



            Ptest = genefun.pickle_loader(pickle_path)
    
    
    
    
    if GFZ_cluster_mode:
        softs_runner.cluster_GFZ_run(kommand_stk)
        
    
    
    try:
        out_sum_fil = acls.exp_summary(exp_path,french_keys=True)
#        acls.plot_cntr_from_expdics(exp_path,exp_path,exp)
#        acls.plot_dist_from_expdics(exp_path,exp_path,exp)
    except:
        print("WARN : fail of summary creation")
        pass


print("restit MK6")

