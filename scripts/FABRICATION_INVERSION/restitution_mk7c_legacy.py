# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 13:04:35 2015

@author: psakicki

v4 : en mode mono Z , on injecte le dZi comme observable
MK6v2 : implementation du sigma_batch_mode (160616)
MK6v3 : fusion avec GEODESEA
MK7 LEGACY : fork directly based on MK6v3 => NO CONFIG FILE
MK7b : for TPX1-GFZ computer
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
import genefun as gf
import re
import copy
import sympy
import pandas as pd
import tabulate
import geodetik as geok
import genefun
from megalib import *
import sys

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
    gene_path = "/home/psakicki/aaa_FOURBI/GNSSA_SIMU_regroupement/GNSSA_TAMBOUILLE_IMPLEMENT_STYLE/"


elif "pgal" in platform.node():
    gene_path = "/dsk/ggsp_pf/PLAYGROUND/psakicki/GNSSA/REBOOT1802_WORK/"

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


        if 1:
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
                                                   'with_zmaster'        : 0},
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

        seeddic_var = { 'k_z_apri'     : list(np.arange(10,12))  ,
                        'k_dZ_apri'    : list(np.arange(20,22))  ,
                        'k_BL_apri'    : list(np.arange(30,32))  }

        seeddic_lis = geok.kwargs_for_jacobian(dict(),seeddic_var)

    else:
        seeddic_lis = [dict()]




    for booldic , sigmadic ,seeddic    in itertools.product(booldic_lis,
                                                             sigmadic_lis,
                                                             seeddic_lis):



        try:
            print("")
            print("===============    BATCH OPTION SUMMARY    ===============")
            print("------- ALL boolean options  (booldic_lis) -------")
            print(booldic_lis)
            print("------- ALL sigma options   (sigmadic_lis) -------")
            print(sigmadic_lis)
            print("------- ALL seed options     (seeddic_lis) -------")
            print(seeddic_lis)
            print("------- CURRENT boolean options  (booldic) -------")
            print(booldic)
            print("------- CURRENT sigma options   (sigmadic) -------")
            print(sigmadic)
            print("------- CURRENT seed options     (seeddic) -------")
            print(seeddic)
            print("=============== END of BATCH OPTION SUMMARY ===============")
            print("")

            if batch_mode:
                print("")
                print('*************** EVALUATE BATCH OPTIONS ***************')
                genefun.eval_a_dict(booldic,globals())
                if with_zmaster:
                    print("WARN : with_zmaster is True and with_barycenter is False => with_barycenter forced to True")
                    print("      (mode is not implemented)")
                    with_barycenter = 1

                bool_4_exp = [str(int(eval(e))) for e in bool_kw_lis]
                bool_4_exp_str = ''.join(bool_4_exp)
                print('******************************************************')
                print("")
            else:
                bool_4_exp_str = ''

            if 'V_4_P_reinject' in list(globals().keys()):
                del V_4_P_reinject

            F = genefun.Tee_frontend(exp_path , exp , bool_4_exp_str )

            # ==============================================
            # INITIALISATION DES PARAMETRES
            # ==============================================

            nPXP_proto = len(bigdico[MNfile])
            idPXP_lis = list(bigdico[MNfile].keys())


            # PARAMETRE BRUIT PXP / Z
            sigma_pxp_apri = 10
            kmeta_pxp_apri = 1 # le K normal c'est l'ID du PXP, le meta c'est pour tous
                               # le Kmeta sera un multiple de 10 par convention

            sigma_z_apri  = 1 # Ne sert qu'en z mono, pour l'apriori du zmono ...
            k_z_apri      = 42 + 10 - 10

            k_dZ_apri     = 1789 + 10 - 10
            sigma_dZ_apri = 10**-6

            k_BL_apri     = 852369 + 10 - 10
            sigma_BL_apri = 10**-2

            ## PARAMETRE SIGMA
            sigma_defo_ASM     = 10**-5
            sigma_defo_BL      = 10**-3
            sigma_defo_zmaster = 10**-8

            ## PARAMETRES ANNEXES
            kriter    = 10**7
            kriterold = 10**10
            iiter = 0 # pour le expdic
            h = 0
            # for jackknife
            keep_ratio   = 0.5
            keep_ratio   = 0.1
            jk_rand_seed = 1452 #Le jk_seed a été ramméné avant l'interprétation du dico



            if sigma_batch_mode:
                print('*************************************')
                print('INFO : sigma_batch_mode AVANT eval : ')
                for k in list(sigmadic.keys()):
                    print(k , eval(k))

                genefun.eval_a_dict(sigmadic,globals())

                print('INFO : sigma_batch_mode APRES eval : ')
                for k in list(sigmadic.keys()):
                    print(k , eval(k))
                print('*************************************')

            if not seed_batch_mode:
                k_z_apri      *= np.abs(np.round(np.log10(sigma_z_apri))  + 11)
                k_dZ_apri     *= np.abs(np.round(np.log10(sigma_dZ_apri)) + 22)
                k_BL_apri     *= np.abs(np.round(np.log10(sigma_BL_apri)) + 33)
                k_z_apri      = int(k_z_apri)
                k_dZ_apri     = int(k_dZ_apri)
                k_BL_apri     = int(k_BL_apri)
            else:

                print('*************************************')
                print('INFO : seed_batch_mode AVANT eval : ')
                for k in list(seeddic.keys()):
                    print(k , eval(k))

                genefun.eval_a_dict(seeddic,globals())

                print('INFO : seed_batch_mode APRES eval : ')
                for k in list(seeddic.keys()):
                    print(k , eval(k))
                print('*************************************')

                k_z_apri      = int(k_z_apri)
                k_dZ_apri     = int(k_dZ_apri)
                k_BL_apri     = int(k_BL_apri)

            if platform.node() in ('calipso','finlande','diamant','VirBox5HPHP'):
                nbproc   = 6   #proco
                iitermax = 3   #itero
            elif platform.node() == 'psakicki-MS-16Y1':
                nbproc   = 3
                iitermax = 1
            elif platform.node() == 'diamant-old':
                nbproc   = 1
                iitermax = 3
            elif platform.node() == 'TPX1-GFZ':
                nbproc   = 3
                iitermax = 4
            elif platform.node() == 'pgal2':
                nbproc   = 4
                iitermax = 3

            elif platform.node() == 'pgal1':
                nbproc   = 7
                iitermax = 3

            elif "kg3" in platform.node():
                nbproc   = 36
                iitermax = 3


            PXPold = np.array([9999999,9999999,9999999] * nPXP_proto)
            PXPnew_stk = []
            ObsASM_lis = []
            Xbato_lis  = []
            TTT_lis    = []

            PXP_lis = []
            PXPapri_lis = []
            fuv = 1
            cleaning_coef = 0
            expdic = collections.OrderedDict()
            iterdic = acls.get_iterdic_in_expdic(expdic,iiter)

            A , N , Ninv , P , B = [],[],[],[],[]

            # PARAMETRES IMPORTANTS

            with_specific_timwin_bool = 0
            if not batch_mode:
                with_BL           = 0
                with_ssp_bilin    = 0
                with_monoZ        = 1
                with_alternat_SSP = 1
                with_barycenter   = 1
                with_decimate_SSP = 1

                with_ssp_bilin    = 0
                with_decimate_SSP = 1
                with_alternat_SSP = True
                with_barycenter   = False
                with_monoZ        = True
                with_BL           = False

                with_ssp_bilin    = 0
                with_decimate_SSP = 1
                with_alternat_SSP = 1
                with_barycenter   = 1
                with_monoZ        = 1
                with_BL           = 1
                with_dzcst        = 1
                with_zmaster      = 0

                with_ssp_bilin    = 0
                with_decimate_SSP = 1
                with_alternat_SSP = 0
                with_barycenter   = 0
                with_monoZ        = 0
                with_BL           = 0
                with_dzcst        = 0
                with_zmaster      = 0

                with_jackknife    = 0
                with_invert_jk    = 1

                with_V_4_P_reinject = 0
                with_time_window    = 0
                with_specific_timwin_bool = 0

                with_specific_timwin_bool =  0
                with_ssp_bilin      = 0
                with_invert_jk      = False
                with_time_window    = 0
                with_V_4_P_reinject = 0
                with_decimate_SSP   = 1
                with_alternat_SSP   = 0
                with_zmaster        = 0
                with_barycenter     = 0
                with_monoZ          = 0
                with_jackknife      = 0
                with_BL             = 0
                with_dzcst          = 0
                with_cleaning       = 0

                with_BL           = 1
                with_ssp_bilin    = 0
                with_monoZ        = 0
                with_alternat_SSP = 0
                with_barycenter   = 0

                with_monoZ   = 0
                with_zmaster = 0

                if MNfile != 'P':
                    with_FB_mode         = 0 # True
                    with_old_style       = 1 # select the F or the B as the time of ping, else TWTT/2
                    with_Backward_select = 0 # WORKS ONLY IF with_FB_mode == False and_with_old_style == True !!!!
                else:
                    with_FB_mode         = 1 # True
                    with_old_style       = 0 # select the F or the B as the time of ping, else TWTT/2
                    with_Backward_select = 0 # WORKS ONLY IF with_FB_mode == False and_with_old_style == True !!!!


            if with_zmaster:
                print("WARN : with_zmaster is True and with_barycenter is False => with_barycenter forced to True")
                print("      (mode is not implemented)")

                with_barycenter = 1


            if with_BL:
                BL = bigdico['B']['d']

            if with_ssp_bilin:
                Z = bigdico['2lin_Z']['d']
                C = bigdico['2lin_C']['d']
            else:
                Z = bigdico['Z']['d']
                C = bigdico['C']['d']


            if with_alternat_SSP:
                alternat_SSP = np.loadtxt(alternat_SSP_path)
                Z = alternat_SSP[:,0]
                C = alternat_SSP[:,1]

            if 1:
                Zraw , Craw = Z , C
                Z,C = ssp.SSP_extrapolate(Z,C,3000,1)

            if with_decimate_SSP:
                print("decimating SSP")
                Zfull , Cfull = Z , C
                Z , C = ssp.SSP_light(Z,C)


            if with_monoZ and with_zmaster:
                print("ERR : with_monoZ and with_zmaster together, impossible combination, skipping ...")
                raise Exception

            if with_zmaster and not with_BL:
                print("ERR : with_zmaster and not with_BL, impossible combination, skipping ...")
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

            for ipxp in sorted(bigdico[MNfile].keys()):

                Mpxp = bigdico[MNfile][ipxp]

                Mdata = Mpxp['d']
                Mcomm = Mpxp['c']
                if MNfile == 'P':
                    Mtab    = Mpxp['t']

                if MNfile == 'N':
                    Xbato = list(Mdata[:,1:4])
                elif MNfile == 'M':
                    Xbato = list(Mdata[:,:3])
                elif MNfile == 'O':
                    Xbato_f = list(Mdata[:,1:4]  )
                    Xbato_b = list(Mdata[:,12:15])
                    Xbato   = list(zip(Xbato_f , Xbato_b))
                elif MNfile == 'P':
                    if with_FB_mode:
                        Xbato_f = list(np.column_stack((Mtab['X_emi_noise'],Mtab['Y_emi_noise'],Mtab['Z_emi_noise'])))
                        Xbato_b = list(np.column_stack((Mtab['X_rec_noise'],Mtab['Y_rec_noise'],Mtab['Z_rec_noise'])))
                        Xbato   = list(zip(Xbato_f , Xbato_b))
                    else:
                        if with_Backward_select:
                            Xbato = list(np.column_stack((Mtab['X_rec_noise'],Mtab['Y_rec_noise'],Mtab['Z_rec_noise'])))
                        else:
                            Xbato = list(np.column_stack((Mtab['X_emi_noise'],Mtab['Y_emi_noise'],Mtab['Z_emi_noise'])))

                Xbato_lis.append(Xbato)

                if MNfile == 'N':
                    ObsASM_load = Mdata[:,-1]
                elif MNfile == 'M':
                    ObsASM_load = Mdata[:,3]
                elif MNfile == 'O':
                    ObsASM_load = Mdata[:,-1]
                elif MNfile == 'P':
                    if with_FB_mode:
                        ObsASM_load = np.array(Mtab['TWTT_noise'])
                    elif with_old_style:
                        if with_Backward_select:
                            ObsASM_load = np.array(Mtab['T_rec_noise'])
                        else:
                            ObsASM_load = np.array(Mtab['T_emi_noise'])
                    else:
                        ObsASM_load = np.array(Mtab['TWTT_noise'] * .5)

                ObsASM_lis.append(ObsASM_load)
                if MNfile == 'N':
                    TTT_lis.append(Mdata[:,0])
                elif MNfile == 'O':
                    TTT_lis.append((Mdata[:,0] , Mdata[:,11]))
                elif MNfile == 'P':
                    if with_FB_mode:
                        TTT_lis.append((Mtab['E_emi_noise'] , Mtab['E_rec_noise']))
                    else:
                        if with_Backward_select:
                            TTT_lis.append(Mtab['E_rec_noise'])
                        else:
                            TTT_lis.append(Mtab['E_emi_noise'])

                PXP = acls.pxp_string_2_array(Mcomm['pxp_coords'])
                PXP_lis.append(PXP)
                nPXP = len(PXP_lis)

                k_pxp_apri = ipxp + kmeta_pxp_apri
                R_pxp_apri = np.random.RandomState(k_pxp_apri)
                if with_monoZ:
                    PXPapri_mono = PXP + np.concatenate((R_pxp_apri.randn(3)[0:2] * sigma_pxp_apri , err_z))
                else:
                    if MNfile in ('N','O'):
                        PXPapri_mono = PXP
                    elif MNfile == 'M':
                        PXPapri_mono = PXP + R_pxp_apri.randn(3) * sigma_pxp_apri
                    elif MNfile == 'P':
                        print('This par must be splited ...')
                        PXPapri_mono = PXP + R_pxp_apri.randn(3) * sigma_pxp_apri

                print("DEBUG ipxp, PXPapri_mono")
                print(ipxp,PXPapri_mono)
                PXPapri_lis.append(PXPapri_mono)

            if with_monoZ:
                for i,pxp in enumerate(PXPapri_lis):
                    if i_pxp_master != i:
                        pxp[2] = PXPapri_lis[i_pxp_master][2]

            PXPZ_lis       = np.array(PXP_lis)[:,-1]
            PXPdZ_lis_true = PXPZ_lis - PXPZ_lis[i_pxp_master]

            noise4PXPdZ = np.random.RandomState(k_dZ_apri).randn(nPXP) * sigma_dZ_apri
            PXPdZ_lis   = noise4PXPdZ +  PXPdZ_lis_true  # <= par defaut c'est le PXP1 qui est maitre
            PXPdZ_lis[i_pxp_master] = 0
            PXPtrue_arr = np.array(PXP_lis)
            shape_arr   = PXPtrue_arr.shape
            PXPapri0_arr = np.array(PXPapri_lis)
            PXPapri0_vct = genefun.vectorialize(PXPapri0_arr)
            PXPapri = PXPapri0_vct

            if with_BL:
                 noise4BLs = np.random.RandomState(k_BL_apri).randn(nPXP , nPXP) * sigma_BL_apri
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

                print("**************DEBUG with_barycenter ********************+++")
                print("PXPbary0     ",PXPbary0      )
                print("PXPbary      ",PXPbary       )
                print("dPXPapri0_arr",dPXPapri0_arr )
                print("dPXPapri0_vct",dPXPapri0_vct )
                print("dPXPapri     ",dPXPapri      )
                print("dPXPapri0    ",dPXPapri0     )
                print("dPXPapri_lis ",dPXPapri_lis  )


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

                print("**************DEBUG with_zmaster ********************+++")
                print("PXPref       ",PXPref        )
                print("dPXPapri0_arr",dPXPapri0_arr )
                print("dPXPapri0_vct",dPXPapri0_vct )
                print("dPXPapri0    ",dPXPapri0     )
                print("dPXPapri     ",dPXPapri      )
                print("dPXPapri_lis ",dPXPapri_lis  )







            if with_time_window:
                ObsASM_lis_orig = list(ObsASM_lis)
                Xbato_lis_orig  = list(Xbato_lis)
                ObsASM_lis  = []
                Xbato_lis   = []
                TTT_lis_4_resid = []


                # Dates en UTC+2
                ssss = dt.datetime(2015,6,22,3)
                eeee = dt.datetime(2015,6,22,7)

                ssss = dt.datetime(2015,6,21,22)
                eeee = dt.datetime(2015,6,22,3)

                # Dates en Pacific Time
                ssss = dt.datetime(2015,6,21,16)
                eeee = dt.datetime(2015,6,21,22)

                ssss = dt.datetime(2015,6,21,12)
                eeee = dt.datetime(2015,6,21,16)

                # Dates en Zulu :)
                ssss = dt.datetime(2015,6,22,1)
                eeee = dt.datetime(2015,6,22,2)

                if with_specific_timwin_bool == True: # CAS TRUE ou 1
                    ssss = dt.datetime(2015,6,21,20)
                    eeee = dt.datetime(2015,6,22,1)
                # ====== CAS INTERMEDIAIRES =======
                elif with_specific_timwin_bool == 21:
                    ssss = dt.datetime(2015,6,21,20)
                    eeee = dt.datetime(2015,6,21,22)
                elif with_specific_timwin_bool == 22:
                    ssss = dt.datetime(2015,6,21,22)
                    eeee = dt.datetime(2015,6,22,0)
                elif with_specific_timwin_bool == 23:
                    ssss = dt.datetime(2015,6,22,0)
                    eeee = dt.datetime(2015,6,22,2)
                elif with_specific_timwin_bool == 24:
                    ssss = dt.datetime(2015,6,22,2)
                    eeee = dt.datetime(2015,6,22,4)

                elif with_specific_timwin_bool == 11:
                    ssss = dt.datetime(2015,6,21,20)
                    eeee = dt.datetime(2015,6,21,21)
                elif with_specific_timwin_bool == 12:
                    ssss = dt.datetime(2015,6,21,21)
                    eeee = dt.datetime(2015,6,21,22)
                elif with_specific_timwin_bool == 13:
                    ssss = dt.datetime(2015,6,21,22)
                    eeee = dt.datetime(2015,6,21,23)
                elif with_specific_timwin_bool == 14:
                    ssss = dt.datetime(2015,6,21,23)
                    eeee = dt.datetime(2015,6,22,0)
                elif with_specific_timwin_bool == 14:
                    ssss = dt.datetime(2015,6,21,23)
                    eeee = dt.datetime(2015,6,22,0)
                elif with_specific_timwin_bool == 15:
                    ssss = dt.datetime(2015,6,22,0)
                    eeee = dt.datetime(2015,6,22,1)
                elif with_specific_timwin_bool == 15:
                    ssss = dt.datetime(2015,6,22,0)
                    eeee = dt.datetime(2015,6,22,1)
                elif with_specific_timwin_bool == 16:
                    ssss = dt.datetime(2015,6,22,1)
                    eeee = dt.datetime(2015,6,22,2)
                elif with_specific_timwin_bool == 17:
                    ssss = dt.datetime(2015,6,22,2)
                    eeee = dt.datetime(2015,6,22,3)
                elif with_specific_timwin_bool == 18:
                    ssss = dt.datetime(2015,6,22,3)
                    eeee = dt.datetime(2015,6,22,4)
                elif with_specific_timwin_bool == 19:
                    ssss = dt.datetime(2015,6,22,4)
                    eeee = dt.datetime(2015,6,22,5)
                elif with_specific_timwin_bool == 99:
                    ssss = dt.datetime(2015,6,22,4)
                    eeee = dt.datetime(2015,6,22,5)
                else:            # CAS FALSE ou 1
                    ssss = dt.datetime(2015,6,22,1)
                    eeee = dt.datetime(2015,6,22,5)

                ssssp = geok.dt2posix(ssss)
                eeeep = geok.dt2posix(eeee)

                if with_smart_SSP:
                    Sdic = genefun.pickle_loader(smart_SSP_Sdic_path)
                    ref_time = ssss + ((eeee - ssss)/2)
                    print('INFO : smart SSP time' , ref_time)
                    tnear , inear = gf.find_nearest(list(Sdic.keys()),ref_time)
                    Z , C = Sdic[tnear][:,0] , Sdic[tnear][:,1]
                    Z,C = ssp.SSP_extrapolate(Z,C,3000,50)


                for ttt,xbato in zip(TTT_lis,Xbato_lis_orig):
                    if MNfile != 'O':
                        _ , outxbato  = geok.time_win_basic(ssssp,eeeep,ttt,xbato)
                    else:
                        _ , outxbato  = geok.time_win_basic(ssssp,eeeep,ttt[-1],xbato)
                        outxbato = [tuple(e) for e in outxbato]
                    Xbato_lis.append(outxbato)
                for iiiddd , (ttt,obsasm) in enumerate(zip(TTT_lis,ObsASM_lis_orig)):
                    if MNfile != 'O':
                        epocoutobsasm , outobsasm = geok.time_win_basic(ssssp,eeeep,ttt,obsasm)
                    else:
                        epocoutobsasm , outobsasm = geok.time_win_basic(ssssp,eeeep,ttt[-1],obsasm)

                    print("avant/apres time win. " , len(obsasm) , len(outobsasm))

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
            print("===================== DEBUT =====================")
            start = genefun.get_timestamp(0)
            acls.print_n_dicwrit( "debut" , str(start) , iterdic , 1 )
            acls.print_n_dicwrit( "nom" , exp , iterdic , 1)
            acls.print_n_dicwrit( "path" , exp_path , iterdic , 1)
            acls.print_n_dicwrit( "plateforme" , gf.get_computer_name() , iterdic , 1)
            acls.print_n_dicwrit( "version python" , sys.version , iterdic , 1)
            acls.print_n_dicwrit( "version numpy" , np.version.version , iterdic , 1)
            acls.print_n_dicwrit( "version scipy" , scipy.version.version, iterdic , 1)

            acls.print_n_dicwrit( "batch_mode" , bool(batch_mode) , iterdic,1)
            acls.print_n_dicwrit( "force_mode" , bool(force_mode) ,  iterdic,1)
            acls.print_n_dicwrit( "with_BL" , bool(with_BL) , iterdic,1)
            acls.print_n_dicwrit( "with_ssp_bilin" , bool(with_ssp_bilin) , iterdic,1)
            acls.print_n_dicwrit( "with_monoZ" , bool(with_monoZ), iterdic,1)
            acls.print_n_dicwrit( "with_alternat_SSP" , bool(with_alternat_SSP), iterdic,1)
            acls.print_n_dicwrit( "with_barycenter" , bool(with_barycenter), iterdic,1)
            acls.print_n_dicwrit( "with_decimate_SSP" , bool(with_decimate_SSP), iterdic,1)
            acls.print_n_dicwrit( "with_dzcst" , bool(with_dzcst), iterdic,1)
            acls.print_n_dicwrit( "with_zmaster" , bool(with_zmaster), iterdic,1)
            acls.print_n_dicwrit( "with_time_window" , bool(with_time_window), iterdic,1)
            acls.print_n_dicwrit( "with_V_4_P_reinject" , bool(with_V_4_P_reinject), iterdic,1)
            acls.print_n_dicwrit( "with_jackknife" , bool(with_jackknife), iterdic,1)

            if with_time_window:
                acls.print_n_dicwrit( "time window start" , (ssss), iterdic,1)
                acls.print_n_dicwrit( "time window end  " , (eeee), iterdic,1)


            if batch_mode:
                acls.print_n_dicwrit( "params variables" , bool_kw_lis , iterdic,1)
                acls.print_n_dicwrit( "params var. de l'exp." , booldic , iterdic,1)

            if with_jackknife:
                acls.print_n_dicwrit( "jackknife inverted" , with_invert_jk , iterdic,1)
                acls.print_n_dicwrit( "jackknife random seed" , jk_rand_seed , iterdic,1)
                acls.print_n_dicwrit( "keep_ratio" , keep_ratio , iterdic,1)

            print("")

            iiter = 1

#            while np.linalg.norm(PXPapri - PXPold) > 5 * 10**-5 and iiter <= iitermax:
            while np.linalg.norm(kriter) > 1 * 10**-5 and iiter <= iitermax:
            #while np.linalg.norm(kriter - kriterold) > 10**-4:
                print("------------ iter No" , iiter , "------------")
                iterdic = acls.get_iterdic_in_expdic(expdic,iiter)

                # preliminaire a la parte BL & ASM
                if with_dzcst:
                    dz_cst = PXPdZ_lis
                else:
                    dz_cst = None


                # Partie ASM



                if not ASM_from_files_mode:
                    ObsASM , ModASM , args_lis_ASM = acls.vectorialize_ASM_multi(PXPapri_lis,
                                                                   ObsASM_lis,
                                                                   Xbato_lis,Z,C,
                                                                   nbprocs=nbproc,
                                                                   dz_cst=dz_cst,
                                                                   ObsASMBoolinp_lis=ObsASMgoodbool)

                    np.savetxt(os.path.join(exp_path,exp + "_" + "ObsASM_1.mat"),ObsASM)
                    np.savetxt(os.path.join(exp_path,exp + "_" + "ModASM_1.mat"),ModASM)

                else:
                    print("WARN : ObsASM & ModASM loaded from a file !!!!")
                    ObsASM = np.loadtxt(os.path.join(exp_path,exp + "_" + "ObsASM_1.mat"))
                    ModASM = np.loadtxt(os.path.join(exp_path,exp + "_" + "ModASM_1.mat"))








                B_ASM = ObsASM - ModASM





                if not ASM_from_files_mode:
                    if with_zmaster:
                        JacobASM = acls.jacob_ASM((PXPref,dPXPapri_lis), ObsASM_lis,
                                                  Xbato_lis,Z,C,h,
                                                  nbprocs=nbproc,monoZ=with_monoZ,
                                                  accur=1,dz_cst=dz_cst,
                                                  ObsASMBoolinp_lis=ObsASMgoodbool)
                    elif not with_barycenter:
                        JacobASM = acls.jacob_ASM(PXPapri_lis,ObsASM_lis,Xbato_lis,Z,C,h,
                                                      nbprocs=nbproc,monoZ=with_monoZ,
                                                      accur=1,dz_cst=dz_cst,
                                                      ObsASMBoolinp_lis=ObsASMgoodbool)
                    else:
                        JacobASM = acls.jacob_ASM((PXPbary,dPXPapri_lis),ObsASM_lis,Xbato_lis,Z,C,h,
                                                      nbprocs=nbproc,monoZ=with_monoZ,
                                                      accur=1,dz_cst=dz_cst,
                                                      ObsASMBoolinp_lis=ObsASMgoodbool)

                    np.savetxt(os.path.join(exp_path,exp + "_" + "JacobASM_1.mat"),JacobASM)

                else:
                    print("WARN : JacobASM loaded from a file !!!!")
                    JacobASM = np.loadtxt(os.path.join(exp_path,exp + "_" + "JacobASM_1.mat"))




                # les baselines sont injectées comme observables => mode contraint et non mode fixé
                #Partie BL
                if with_BL:
                    if with_barycenter:
                        ObsBL,ModBL = acls.vectorialize_BL(BLnoised,dPXPapri_lis,dz_cst)
                    else:
                        ObsBL,ModBL = acls.vectorialize_BL(BLnoised,PXPapri_lis,dz_cst)
                    B_BL    = ObsBL - ModBL

                    print("")
                    print("DEBUG BL")
                    print("**********************************+++")
                    print("B_BL,ObsBL,ModBL")
                    print(B_BL,ObsBL,ModBL)
                    print("**********************************+++")
                    print("")

                    #JacobBL = acls.jacob_BL(PXPapri_lis,0,1,dz_cst)

                    JacobBL = acls.jacob_BL(PXPapri_lis,with_monoZ,with_barycenter,dz_cst)


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

                print("fin du remplissage, debut du redesign + gestion des poids")
                A = JacobASM
                B = B_ASM



                #P = gf.diagonalize(0.001,np.max(A.shape))
                # Bloc des poids bien pas beau comme il faut
                # réformé le 160603
                #if not with_V_4_P_reinject or not "V_4_P_reinject" in globals().keys():
                #    Ptmp = [1 / ((10**-3)**2)] * np.max(A.shape)
                #else:
                #    Ptmp = 1 / ((np.array(V_4_P_reinject) * 1) ** 2)
                #
                #if with_BL:
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

                if not with_V_4_P_reinject or not "V_4_P_reinject" in list(globals().keys()): # V_4_P_reinject (2nd teste) désigne la liste et non le bouleen, pas touche
                    sigma_4_P = sigma_4_P + [sigma_defo_ASM]
                    len_4_P.append(np.max(A.shape))
                else:
                    sigma_4_P = sigma_4_P + list(V_4_P_reinject)
                    len_4_P   = len_4_P     + len(V_4_P_reinject)  * [1]

                if with_BL:
                    #_,_,Ptmp = geok.weight_mat([10,1],[len(ObsASM),len(ObsBL)],
                    #                        fuv,sparsediag=True)
                    #Ptmp = gf.diagonalize(0.001,len(ObsBL))
                    sigma_4_P = sigma_4_P + [sigma_defo_BL]
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

                _,_,P = geok.weight_mat(sigma_4_P,len_4_P,fuv,sparsediag=1)
                #P[-7,-7] = 2.370362140050853e+29


                print("fin des poids ")

#                if not batch_mode:
#                    A_stk.append(A)

                # la contrainte de barycentre est injectée suivant la methode fixée d'Helmert
                if with_barycenter:
                    G     = acls.constraints_matrix(len(PXP_lis),with_barycenter,with_monoZ) #MODIF
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

                print("conversion en sparse ")
                Asprs = sprs.csc_matrix(A)
                Bsprs = sprs.csc_matrix(B).T
                Psprs = sprs.csc_matrix(P)
                print("fin de la conversion en sparse ")

                Nsprs = (Asprs.T).dot(Psprs).dot(Asprs)
                N     = Nsprs.toarray()

#                N = np.dot(np.dot(A.T , P ) , A)

                # Liquidation des matrices inutiles
                JacobASM , ObsASM = [] , []

                # debut contraintes
                if with_barycenter:
                    Gsiz = G.shape[0]
                    O = np.zeros((G.shape[0],G.shape[0]))
                    N = np.vstack((N,G))
                    N = np.hstack((N,np.hstack((G,O)).T))
                # fin contraintes
                Ninv = scipy.linalg.inv(N)
#                AtPB = A.T * P * B
#                AtPB = np.dot(np.dot(A.T , P ) , B)
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

                if with_barycenter:
                    if not with_monoZ:
                        dX = dX[:-3]
                    else:
                        dX = dX[:-2]

                gruikmode = 0 # Discontinué

                if with_zmaster and gruikmode:
                    # manip empirique de gros porc pour corriger la non nullité
                    # de la Somme des dZ ... (160609)
                    # On remarque empiriquement que sum dZ = une valeur dZcorrection
                    # que l'on ajoute ensuite aux dZ individuel
                    # Le gruik mode est discontinué mais a grandement aidé
                    # au debugging, on le garde en souvenir (160610)
                    PXPbarynew_gruik = np.array(PXPbary) + dX[:3]
                    PXPbary_gruik    = np.array(PXPbarynew_gruik)

                    dPXPnew_gruik  = dPXPapri + dX[n_lagrangian:]
                    bary_gruik = dPXPnew_gruik.reshape(shape_arr)
                    dZcorr_gruik   = acls.barycenter_calc(bary_gruik)[-1]

                    for i in range(nPXP):
                        dX[ 3+ i*3+2 ] = dX[ 3+ i*3+2 ] - dZcorr_gruik

                if with_monoZ:
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
#                if with_BL:
#                    if with_barycenter:
#                        ObsBL_V,ModBL_V = acls.vectorialize_BL(BL,dPXPapri_lis)
#                    else:
#                        ObsBL_V,ModBL_V = acls.vectorialize_BL(BL,PXPapri_lis)
#
#                    V_BL = ObsBL_V - ModBL_V
#                    Vnaif = np.concatenate((Vnaif,V_BL))


                Vrigour = B - A.dot(dXorig[:A.shape[1]])

                V = Vrigour

                if cleaning_coef != 0 and not ObsASMgoodbool is None: # pour quoi cleaning_coef != 0, et pourquoi ici ? chépa (160922)

                    lenObsASM = [np.sum(o) for o in ObsASMgoodbool]
                    V_per_ObsASM = gf.sublistsIt(V,lenObsASM,True)

                    ObsASMgoodbool_tmp = []
                    Vcleaned = []
                    if with_V_4_P_reinject:
                        V_4_P_reinject = []
                    if with_cleaning:
                        for iii  , (boolping  , Voasm) in enumerate(zip(ObsASMgoodbool , V_per_ObsASM )):
                            print("ping valides AVANT nettoyage",np.sum(boolping),'for PXP',iii)
                            indices_true  = np.where(boolping)[0]
                            indices_false = np.where(np.logical_not(boolping))[0]

                            actualised_bool = boolping[indices_true] * (np.abs(np.array(Voasm)) < cleaning_coef * np.std(Voasm))
                            boolping_new = np.array(boolping)
                            boolping_new[indices_true] = actualised_bool

                            ObsASMgoodbool_tmp.append(boolping_new)

                            Vcleaned = Vcleaned + list(Voasm[actualised_bool])

                            if with_V_4_P_reinject:
                                V_4_P_reinject = V_4_P_reinject + list(Voasm[actualised_bool])
                            print("ping valides APRES nettoyage",np.sum(boolping_new),'for PXP',iii)

                        ObsASMgoodbool = ObsASMgoodbool_tmp

                    if with_V_4_P_reinject:
                        V_4_P_reinject = np.array(V_4_P_reinject)

                if not with_BL:
                    fuv = geok.fuv_calc(V,A,P)
                else:
                    fuv = geok.fuv_calc(V,A,P)

                acls.print_n_dicwrit('f.u.v.' , fuv , iterdic , 1)
                acls.print_n_dicwrit('coords. apriori' , np.array(PXPapri_lis_old) ,
                                     iterdic )

                acls.print_n_dicwrit('dX' , dX , iterdic )
                acls.print_n_dicwrit('som. dX' , np.sum(dX) , iterdic , 1)
                acls.print_n_dicwrit('som. abs dX' , np.sum(np.abs(dX)) , iterdic , 1)

                PXPnew_arr = PXPnew.reshape(shape_arr)

                acls.print_n_dicwrit( 'nouvelles coords.' , PXPnew_arr , iterdic )
                if not batch_mode:
                    PXPnew_stk.append(PXPnew)
                with np.errstate(invalid='ignore'):
                    sigmas = np.sqrt(np.diag(Ninv) * fuv)
                acls.print_n_dicwrit( 'sigmas' , sigmas , iterdic )

                barybrut = acls.barycenter_calc(PXPnew.reshape(shape_arr))
                barybrut_old = acls.barycenter_calc(PXPold.reshape(shape_arr))

                kriterold = kriter
                kriter = np.linalg.norm(barybrut_old - barybrut)
                kriter = np.sum(np.abs(dX))

                iterdic['varcovar'] = Ninv
                acls.print_n_dicwrit( "kriter/kriterold" , [ kriter , kriterold ] , iterdic)
                acls.print_n_dicwrit( "ecart 3D a la postion vraie en distance" ,
                                np.linalg.norm( PXPnew_arr - PXPtrue_arr,axis=1) , iterdic)
                acls.print_n_dicwrit( "ecart 2D a la postion vraie en distance" ,
                                np.linalg.norm((PXPnew_arr - PXPtrue_arr)[:,:2],axis=1) , iterdic)

                acls.print_n_dicwrit( "ecart a la postion vraie en coordonnees"  ,
                                PXPnew_arr - PXPtrue_arr , iterdic )

                acls.print_n_dicwrit( "BLs Vraies " , acls.print_BL(PXPtrue_arr) , iterdic)

                if with_BL:
                    acls.print_n_dicwrit( "BLs Observees" , acls.print_BL(ObsBL,0) , iterdic)

                acls.print_n_dicwrit("BL des coords a prioris('Mod')",
                                     acls.print_BL(PXPapri0_arr),iterdic)

                acls.print_n_dicwrit("BLs News à l'issue de l'inversion",
                                     acls.print_BL(PXPnew_arr),iterdic)

                acls.print_n_dicwrit("Barycentre 'brut' : Sum Xpxp / Npxp",
                                     barybrut , iterdic)

                baryvrai = acls.barycenter_calc(PXPtrue_arr)
                acls.print_n_dicwrit("ecart 3D bary brut/vrai en distance" ,
                np.linalg.norm(barybrut - baryvrai) , iterdic)

                acls.print_n_dicwrit("ecart 2D bary brut/vrai en distance" ,
                np.linalg.norm(barybrut[0:2] - baryvrai[0:2]) , iterdic )

                acls.print_n_dicwrit("ecart au bary brut/vrai en coords.",
                barybrut - baryvrai , iterdic )

                if with_barycenter:
                    acls.print_n_dicwrit("Barycentre estime dans le MC",
                                         PXPbary , iterdic )

                    acls.print_n_dicwrit("Barycentre des dXs dans le MC",
                    acls.barycenter_calc(dPXPnew.reshape(shape_arr)) , iterdic)

                    if np.sum(np.abs(acls.barycenter_calc(dPXPnew.reshape(shape_arr)))) > 10**-5:
                        print("WARN : som barycentre des dXs dans le MC > 10**-5 !!!!")



                    acls.print_n_dicwrit("ecart 3D bary MC/vrai en distance",
                    np.linalg.norm(PXPbary - baryvrai) , iterdic)

                    acls.print_n_dicwrit("ecart 2D bary MC/vrai en distance",
                    np.linalg.norm(PXPbary[0:2] - baryvrai[0:2]) , iterdic)

                    acls.print_n_dicwrit("ecart au bary MC/bary vrai en coords.",
                    PXPbary - baryvrai , iterdic )

                if with_zmaster:
                    acls.print_n_dicwrit("ecart à la balise de réference",
                                         PXPnew_arr[:,2] - PXPnew_arr[i_pxp_master,2],
                                         iterdic)
                    acls.print_n_dicwrit("ecart en input à la balise de réference",
                                         PXPdZ_lis,
                                         iterdic)
                    acls.print_n_dicwrit("ecart vrai à la balise de réference",
                                         PXPdZ_lis_true,
                                         iterdic)

                iiter=iiter+1

            print("===================== FINAL =====================")
            end = genefun.get_timestamp(0)
            acls.print_n_dicwrit("fin"   , str(end) , iterdic , 1 )
            acls.print_n_dicwrit("duree" , str(end - start) , iterdic , 1 )
            acls.print_n_dicwrit("nom" , exp , iterdic , 1)
            acls.print_n_dicwrit( "path" , exp_path , iterdic , 1)
            acls.print_n_dicwrit( "plateforme" , gf.get_computer_name() , iterdic , 1)
            acls.print_n_dicwrit( "version python" , sys.version , iterdic , 1)
            acls.print_n_dicwrit( "version numpy" , np.version.version , iterdic , 1)
            acls.print_n_dicwrit( "version scipy" , scipy.version.version, iterdic , 1)
            acls.print_n_dicwrit("batch_mode" , bool(batch_mode) , iterdic , 1)
            acls.print_n_dicwrit("force_mode" , bool(force_mode) , iterdic , 1)
            acls.print_n_dicwrit("with_BL" , bool(with_BL) , iterdic , 1)
            acls.print_n_dicwrit("with_ssp_bilin" , bool(with_ssp_bilin) , iterdic , 1)
            acls.print_n_dicwrit("with_monoZ" , bool(with_monoZ) , iterdic , 1)
            acls.print_n_dicwrit("with_alternat_SSP" , bool(with_alternat_SSP) , iterdic , 1)
            acls.print_n_dicwrit( "with_barycenter" , bool(with_barycenter), iterdic,1)
            acls.print_n_dicwrit( "with_decimate_SSP" , bool(with_decimate_SSP), iterdic,1)
            acls.print_n_dicwrit( "with_dzcst" , bool(with_dzcst), iterdic,1)
            acls.print_n_dicwrit( "with_zmaster" , bool(with_zmaster), iterdic,1)
            acls.print_n_dicwrit( "with_time_window" , bool(with_time_window), iterdic,1)
            acls.print_n_dicwrit( "with_V_4_P_reinject" , bool(with_V_4_P_reinject), iterdic,1)
            acls.print_n_dicwrit( "with_jackknife" , bool(with_jackknife), iterdic,1)
            acls.print_n_dicwrit("poids ASM",sigma_defo_ASM     , iterdic )
            acls.print_n_dicwrit("poids BL", sigma_defo_BL      , iterdic )
            acls.print_n_dicwrit("poids dZ", sigma_defo_zmaster , iterdic )

            acls.print_n_dicwrit("bruit PXP apriori"     , sigma_pxp_apri  , iterdic )
            acls.print_n_dicwrit("seed bruit PXP apriori", kmeta_pxp_apri  , iterdic )
            acls.print_n_dicwrit("bruit dZ"              ,  sigma_dZ_apri  , iterdic )
            acls.print_n_dicwrit("seed bruit dZ"         , k_dZ_apri       , iterdic )
            acls.print_n_dicwrit("bruit  Z"              , sigma_z_apri    , iterdic )
            acls.print_n_dicwrit("seed bruit Z"          , k_z_apri        , iterdic )
            acls.print_n_dicwrit("bruit BL"              , sigma_BL_apri   , iterdic )
            acls.print_n_dicwrit("seed bruit BL"         , k_BL_apri       , iterdic )

            if with_time_window:
                acls.print_n_dicwrit( "time window start" , (ssss), iterdic,1)
                acls.print_n_dicwrit( "time window end  " , (eeee), iterdic,1)
            if batch_mode:
                acls.print_n_dicwrit( "params variables" , bool_kw_lis , iterdic,1)
                acls.print_n_dicwrit( "params var. de l'exp." , booldic , iterdic,1)

            if with_jackknife:
                acls.print_n_dicwrit( "jackknife inverted" , with_invert_jk , iterdic,1)
                acls.print_n_dicwrit( "jackknife random seed" , jk_rand_seed , iterdic,1)
                acls.print_n_dicwrit( "keep_ratio" , keep_ratio , iterdic,1)


            acls.print_n_dicwrit("Nb pings     " , len(Xbato) , iterdic , 1)
            acls.print_n_dicwrit("Nb PXPs      " , nPXP , iterdic , 1)
            acls.print_n_dicwrit("Taille Jacob." , A.shape , iterdic , 1)
            acls.print_n_dicwrit("Taille Resid." , V.shape , iterdic , 1)
            acls.print_n_dicwrit("Nb iter.     " , iiter  , iterdic , 1)

            F.stop()

            expdic[-1] = expdic[max(expdic.keys())]

            # AD HOC pickle saver
            # save_obj_as_file(objin,pathin,prefix,ext='.exp',suffix='')
            genefun.save_obj_as_file(expdic , exp_path , exp, suffix=bool_4_exp_str)
            # GENERIC pickle saver not compatible with acoustic inversion nomenclature
            # pickle_saver(datain , outdir=None , outname=None , ext='.pik' , timestamp=False,full_path=None)
            # genefun.pickle_saver(expdic,exp_path, "_".join((exp,bool_4_exp_str)) , ext='.pik',timestamp=True)

            Vnorma      = np.array(V / P.diagonal())
            Vnormclean  = Vnorma[np.abs(Vnorma) < 3 * np.std(Vnorma)]

            protopath = os.path.join( exp_path , exp + '_' + bool_4_exp_str + '_' + gf.get_timestamp() )

            if plot_mode:
                plt.clf()
                plt.subplots_adjust(left=0.2, right=0.9, top=0.9, bottom=0.2)

                n,bins,patches = plt.hist(Vnormclean ,100,normed=1)
                gauss = scipy.stats.norm.pdf(bins,np.mean(Vnormclean),np.std(Vnormclean))
                plt.plot(bins,gauss)

                plt.xlabel('TWTT residuals (s)')

                plt.savefig(protopath + 'histVnorma.png')
                plt.savefig(protopath + 'histVnorma.pdf')


                plt.clf()
                plt.subplots_adjust(left=0.2, right=0.9, top=0.9, bottom=0.2)

                n,bins,patches = plt.hist(V ,100,normed=1)
                gauss = scipy.stats.norm.pdf(bins,np.mean(V),np.std(V))
                plt.plot(bins,gauss)

                plt.xlabel('TWTT residuals (s)')
                plt.ylabel('occurence')

                plt.savefig(protopath + 'histV.png')
                plt.savefig(protopath + 'histV.pdf')


                plt.clf()
                plt.subplots_adjust(left=0.2, right=0.9, top=0.9, bottom=0.2)

                n,bins,patches = plt.hist(V ,100,normed=1)
                gauss = scipy.stats.norm.pdf(bins,np.mean(V),np.std(V))
                plt.plot(bins,gauss)
                plt.xlim((- np.std(V) * 3 , np.std(V) * 3))

                plt.xlabel('TWTT residuals (s)')
                plt.ylabel('occurence')

                plt.savefig(protopath + 'histV.png')
                plt.savefig(protopath + 'histV.pdf')


                plt.clf()
                plt.subplots_adjust(left=0.2, right=0.9, top=0.9, bottom=0.2)

                lenObsASM_4_plot = [np.sum(o) for o in ObsASMgoodbool]
                V_per_ObsASM_4_plot = gf.sublistsIt(Vcleaned,lenObsASM_4_plot,True)
                for iidplot,(Arrbool , Ttt , Vv ) in enumerate(zip(ObsASMgoodbool,TTT_lis_4_resid,V_per_ObsASM_4_plot )):
                    plt.plot(geok.posix2dt(Ttt,1)[Arrbool],Vv,'.',label = 'ID' + str(idPXP_lis[iidplot]))

                    plt.ylabel('TWTT residuals (s)')
                    plt.legend()

                fig = plt.gcf()
                fig.autofmt_xdate()

                plt.savefig(protopath + 'temporalV.png')
                plt.savefig(protopath + 'temporalV.pdf')

                plt.clf()
                plt.subplots_adjust(left=0.2, right=0.9, top=0.9, bottom=0.2)

                fig , AXES = plt.subplots(len(idPXP_lis),1)
                fig.set_size_inches(8.27,11.69)
                lenObsASM_4_plot = [np.sum(o) for o in ObsASMgoodbool]
                V_per_ObsASM_4_plot = gf.sublistsIt(Vcleaned,lenObsASM_4_plot,True)

                smartV_stk = []
                COL = ['r','g','y','b']
                for iidplot,(Arrbool , Ttt , Vv , ax , col) in enumerate(zip(ObsASMgoodbool,TTT_lis_4_resid,V_per_ObsASM_4_plot , AXES,COL)):
                    ax.plot(geok.posix2dt(Ttt,1)[Arrbool],Vv,'.',label = 'ID' + str(idPXP_lis[iidplot]),color=col)
                    smartV_stk.append(np.column_stack(([iidplot] * len(Vv),Ttt[Arrbool],Vv)))

                    ax.set_ylabel('TWTT residuals (s)')
                    ax.legend()

                np.savetxt(protopath + '.smartV' , np.vstack(smartV_stk) )


                fig.autofmt_xdate()

                plt.savefig(protopath + 'temporalV2.png')
                plt.savefig(protopath + 'temporalV2.pdf')


            np.savetxt(protopath + '.V' , V)

            geok.chi2_test_lsq(V,A,P)

            #%%
            if plot_mode:
                plt.clf()
                if with_FB_mode:
                    plt.plot(np.vstack(Xbato_f)[:,0], np.vstack(Xbato_f)[:,1], 'b-')
                    #plt.plot(np.vstack(Xbato_b)[:,0], np.vstack(Xbato_b)[:,1], 'rx')
                else:
                    print("WARN : plot traject potentiellement foireux")
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
                        reload(acls)
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

                plt.axis('equal')
                plt.savefig(protopath + "_geo.png")
                plt.savefig(protopath + "_geo.pdf")

            if plot_mode:
                plt.close('all')

            print("STOP")
            print("FIN SANS ECHEC")

        except Exception as e:
            if force_mode:
                print("ERR : exception ... mais on continue !")
                print(e)
                F.stop()
                continue
            else:
                print("ERR : exception ... on arrête tout !")
                print(e)
                F.stop()
                raise
    try:
        out_sum_fil = acls.exp_summary(exp_path,french_keys=True)
#        acls.plot_cntr_from_expdics(exp_path,exp_path,exp)
#        acls.plot_dist_from_expdics(exp_path,exp_path,exp)
    except:
        print("WARN : fail of summary creation")
        pass





print("restit MK6")
