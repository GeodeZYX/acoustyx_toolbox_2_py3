# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 18:22:29 2015

Fabrication des shoots acoustiques + la trajectoire

Utilsiation comme base des scripts 'move as SSP'
mais refonte du code pour raison de propreté

On appelle (à la manière de GAMIT)
les Mfiles : la grosse matrice par PXP
les Zfiles et Cfiles : le SSP
le BLfile : la matrice de baseline

160419 : une version MK3.2 histoire de faire un backup
         avec la volonté de coder l'erreur de datation

160626 : integration du pseudoSSF 

160627 : v6 plus de if dans le header pour faciliter le batch

@author: psakicki
"""

import glob
import acouclass as acls
import numpy as np
import raytrace as rt
from megalib import *
import matplotlib.pyplot as plt
import multiprocessing as mp

# ===============================
# PARAMETERS
# ===============================

if platform.node() == 'calipso':
    prm_ssp_file_path = '/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/exemple/SSP/SSP_NOAA_dep5645_20030608000000'
    path_gene         = '/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/working'
    prm_sspt_dic_path = "/home/psakicki/THESE/bd_CTD/NEW_DATA_2015/MOVE/MOVE_reform_1507/SIMPLE/sspdic_sensor.pik"
    prm_Gdic_path     = "/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/exemple/Gdic/Gdic.pik"
    prm_Gdic_path     = "/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/exemple/Gdic/Gdic10p6.pik"
    prm_Gdic_path     = "/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/exemple/Gdic/Gdic_10p5_3500_10.pik"
    prm_Gdic_path     = "/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/exemple/Gdic/Gdic_realik15.pik"
    prm_Gdic_path     = '/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/exemple/Gdic/Gdic_proxyMOVE_median2_10.pik'
    

    prm_Kfile_path    = "/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/exemple/Kfiles_kourents_gradients/kfile_simplest.K"
    prm_Kfile_path    = "/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/exemple/Kfiles_kourents_gradients/kfiletestpoly.K"
    prm_Kfile_path    = "/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/exemple/Kfiles_kourents_gradients/kfile_poly_proxy_MOVE.K"
    
    prm_Udic_path = "/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/exemple/Udic/Udictest.pik"

elif platform.node() == 'psakicki-MS-16Y1':
    prm_ssp_file_path = '/home/psakicki/Documents/CODES/acoustyx_toolbox_2/exemple/SSP/SSP_NOAA_dep5645_20030608000000'
    path_gene         = '/home/psakicki/Documents/CODES/acoustyx_toolbox_2/working'
    prm_sspt_dic_path = "/media/psakicki/D34F-24E3/plots/SSPs/sspdic"
    prm_sspt_dic_path = "/home/psakicki/Documents/CODES/acoustyx_toolbox_2/exemple/SSPdic/sspdic_sensor.pik"
    prm_Gdic_path     = "/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/exemple/Gdic/Gdic.pik"
    prm_Udic_path = "/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/exemple/Udic/Udictest.pik"


exp_prefix  = 'tpo1' # for temporal
exp_prefix  = 'tst1'
exp_prefix  = 'repris1'

exp_prefix  = 'batc2b_deriv'
exp_prefix  = 'batc2d_deriv_megadecalZ'
exp_prefix  = 'batc3b_2000_bw_deriv_decalZ'
exp_prefix  = 'batc2c_deriv_decalZ_notrajnoise_stdpingnoise'
exp_prefix  = 'batc2c_deriv_decalZ_nopingnoise_trajnoise5decim'
exp_prefix  = 'batc4_nadir'
exp_prefix  = 'ADEL1ASYNC0'
exp_prefix  = 'RC1_SSF'
exp_prefix  = 'GRAD_mk2d'

exp_prefix  = 'SSFssrealik_1'
exp_prefix  = 'ADELLL666'

exp_prefix  = 'Grid4Gdic_10p5_3500'
exp_prefix  = 'SSFrealik_3_10p5_3500_smallarray'
exp_prefix  = 'SSFrealik_2_smallarray'

exp_prefix  = 'Grid4Gdic_proxyMOVE_median2'
exp_prefix  = 'SSFrealik_4_proxyMOVE_median2'

exp_prefix  = 'SSFinterpo_1'
exp_prefix  = 'SSFinterpo_1_mode_SD_pour_kompar_avecSSP_issu_du_SSF'

exp_prefix  = 'test_de_rephasage'

exp_suffix  = ''

# PXP COORDINATES
PXP1 = np.array([-2500,-2500,4010])
PXP2 = np.array([2500,2500,4005])
PXP3 = np.array([-2500,2500,3098])
PXP4 = np.array([2500,-2500,3085])

PXP1 = np.array([-2500,-2500,4010])
PXP2 = np.array([2500,2500,4005])
PXP3 = np.array([-2500,2500,3998])
PXP4 = np.array([2500,-2500,3985])

PXP1 = np.array([-500,-500,4010])
PXP2 = np.array([ 500, 500,4005])
PXP3 = np.array([-500, 500,3998])
PXP4 = np.array([ 500,-500,3985])

ddd  = 1000 
PXP1 = np.array([-ddd,-ddd,4010])
PXP2 = np.array([ ddd, ddd,4005])
PXP3 = np.array([-ddd, ddd,3998])
PXP4 = np.array([ ddd,-ddd,3985])

zPXP = 5000
PXP1 = np.array([-2500,-2500,zPXP])
PXP2 = np.array([2500 ,2500 ,zPXP])
PXP3 = np.array([-2500,2500 ,zPXP])
PXP4 = np.array([2500 ,-2500,zPXP])


zPXP = 4000
PXP1 = np.array([-500,-500,zPXP+10])
PXP2 = np.array([ 500, 500,zPXP-20])
PXP3 = np.array([-500, 500,zPXP+30])
PXP4 = np.array([ 500,-500,zPXP-40])

zPXP = 4000
PXP1 = np.array([-2500,-2500,zPXP+10])
PXP2 = np.array([2500 , 2500,zPXP-20])
PXP3 = np.array([-2500, 2500,zPXP+30])
PXP4 = np.array([2500 ,-2500,zPXP-40])


#PXP1 = np.array([0,0,4000])
PXP_lis = [PXP1,PXP2,PXP3,PXP4]

#PXP_lis = [PXP1]
tat     = 0

#PXP_lis = [PXP1]

#ligne obligatoire pour tout déporter à 10000 afin d'éviter les effets de bord
PXP_lis = [pxp + np.array([10000,10000,0]) for pxp in PXP_lis]
PXP_arr = np.vstack(PXP_lis)

prm_ssp_decimate = True

# TEMPORAL
# L'UNITE DE REF EST LA SECONDE !!!
prm_with_temporal   = 0
prm_tempor_start    = dt.datetime(2009, 6, 23, 0, 00) #dt.datetime(2008, 1, 4, 19, 20)
prm_tempor_len      = 86400
prm_tempor_ctd_time = prm_tempor_start #dt.datetime(2008, 1, 5, 19, 20)
prm_vit             = 1   #1 * 3600. / 500.  # en m/s
prm_epoch_init      = 0 

# TIDE (only activ if temporal activ)
prm_with_tide      = True
prm_tide_amplitude = 0.057840665308627685

# NOISE
prm_sigma_x = 0.001 
prm_sigma_y = 0.001 
prm_sigma_z = 0.002
prm_sigma_t = 0

prm_sigma_x = 0.10
prm_sigma_y = 0.10
prm_sigma_z = 0.15

prm_sigma_x = 0.50
prm_sigma_y = 0.50
prm_sigma_z = 1.

prm_sigma_x = 0.03
prm_sigma_y = 0.03 
prm_sigma_z = 0.05

prm_sigma_x = 0 
prm_sigma_y = 0 
prm_sigma_z = 0
prm_sigma_t = 0

# OFFSET
prm_add_offset = False
prm_offset     = [0.1,0.1,0.1]

imporved_noising            = 1
prm_null_water_column_noise = 1
prm_coef_sigma_zones        = 0
 
prm_with_ssf           = 1
prm_bypass_ssf_rec     = 1 # La recherche du eigenray retour étant très longue, on peut la bypasser et prendre la valeur de l'émission
prm_ssf_realistik_grad = 0
prm_SSF_based_on_interpolation = 0


# if NOT realistic grad
prm_eiko_h          = 1000
prm_eiko_resotype   = 'rkck'
prm_x_grad          = 0 #10**-5 #2.82625440384e-06   #1 * 10**-4
prm_y_grad          = 0       
prm_z_mixed         = 3500 #500

# =======================

prm_with_pseudo_ssf = 0 # ne doit être ON que si prm_with_ssf est OFF
                        # On exploite alors le Gdic

with_async = 1

###if imporved_noising:
mono_epoc        = False
prm_zones_bound  = [100,500,1500,2500]
prm_sigma_zones  =          [10**-4,10**-4,10**-4,10**-4,10**-4] 
prm_sigma_zones  = np.array([10**-2,10**-2,10**-3,10**-3,10**-4])

prm_zones_bound  = [85.0, 128.0, 202.0, 336.0, 521.5, 737.0, 977.5, 1233.0, 1498.0, 1773.5, 2054.0, 2354.0, 2674.5, 2995.0, 3314.5, 3634.5, 3954.5, 4249.5, 4520.0, 4775.5]
# écart type GLOBAL en terme de vitesse (m/s) maxi pour chaque sensors
prm_sigma_zones  = [10.242326276256273, 9.0194988496170989, 9.037603915864123, 6.5552606451821811, 5.3675675636850562, 2.6245083826908431, 1.995944311352654, 0.85310239595884385, 0.59606180974035217, 0.68827807182587697, 0.69401959561788751, 0.561718589619987, 0.57023253817564767, 0.7093758989173311, 0.49801892775203543, 0.74872331306260986, 0.59637578194047847, 0.3754758173676534, 0.2714521430911136, 0.14603678520159574, 0.22089892006374223]
# écart type GLOBAL en terme de vitesse (m/s) du jour 2009, 10, 2 PIRE JOUR
prm_sigma_zones  = [1.1294490768056602, 2.3271346245141826, 3.8989585830945623, 2.7872329075498423, 1.0992552152802941, 0.59283006896481016, 0.8441119312422859, 0.26105981242779575, 0.2719900156906187, 0.29139880761415776, 0.21710974820564075, 0.25541689969306347, 0.18210410560979798, 0.19192284887362085, 0.16249427940069811, 0.13108661713705372, 0.12538376620953084, 0.11622005621190985, 0.094064705767912038, 0.048114545814723304, 0.021960003178055625]
# écart type des RESIDUS en terme de vitesse (m/s) du jour 2009, 10, 2 PIRE JOUR
prm_sigma_zones  = [0.2056013911230388, 0.41526282258268971, 0.48261023939831421, 0.38557432287789739, 0.3530380782748811, 0.24376327806241502, 0.23364579787961265, 0.122484875302312, 0.033055723731374116, 0.045232104608395371, 0.028361756873242223, 0.015833337921709083, 0.033274359263870901, 0.03373283408434026, 0.0095287555473885684, 0.016874797279464185, 0.0093775434859970502, 0.0080870122548441174, 0.0057296494925734576, 0.0078149294749356602, 0.0062546386727679948]

prm_sigma_t_hdwr = 0        # erreur "hardware" résiduel propre à chaque PXP
prm_sigma_t_hdwr = 1*10**-6 # erreur "hardware" résiduel propre à chaque PXP

prm_zones_bound  = [500]
prm_sigma_zones  = [10,10] 

# écart type des RESIDUS en terme de vitesse (m/s) du jour 2007-04-24 JOUR MEDIAN
prm_sigma_zones =  np.array(list(reversed([0.0084860975169150348, 0.0077746395258356284, 0.0056111516219907011, 0.010427777463161744, 0.016723956769679072, 0.020665159598152891, 0.02244101699021828, 0.026310974740559833, 0.053630364030446659, 0.018654381782600555, 0.0229119090855158, 0.034312443410280405, 0.067845114424280942, 0.052028408496052253, 0.27195619026242834, 0.3390498221574656, 0.31754658978250788, 0.21331715011416538, 0.70105577322691237, 0.36711343216037745, 0.17608383536233305])))
# écart type des RESIDUS en terme de vitesse (m/s) du jour 2010, 2, 25, JOUR MEILLEUR
prm_sigma_zones =  np.array(list(reversed([0.010616052278649434, 0.0051059136080509188, 0.0029715997188615646, 0.005880264264937993, 0.0078627834091491988, 0.012549859648335077, 0.014427741950775038, 0.024639823545407852, 0.029128444426642083, 0.036350431577136674, 0.028188024108818176, 0.02725764678877789, 0.042119526127577585, 0.047136849776481846, 0.13346235979606466, 0.14713275558898559, 0.18252139375909865, 0.1467923524657169, 0.18407545304758163, 0.098006631286506524, 0.033948356716065296])))
# écart type des RESIDUS en terme de vitesse (m/s) du jour 2009, 6, 23, JOUR WORST    
prm_sigma_zones =  np.array(list(reversed([0.0062546386727679948, 0.0078149294749356602, 0.0057296494925734576, 0.0080870122548441174, 0.0093775434859970502, 0.016874797279464185, 0.0095287555473885684, 0.03373283408434026, 0.033274359263870901, 0.015833337921709083, 0.028361756873242223, 0.045232104608395371, 0.033055723731374116, 0.122484875302312, 0.23364579787961265, 0.24376327806241502, 0.3530380782748811, 0.38557432287789739, 0.48261023939831421, 0.41526282258268971, 0.2056013911230388])))

#prm_sigma_zones_OK = list(prm_sigma_zones)
prm_zones_bound = list(reversed([4902.5, 4654.3, 4384.2, 4099.2, 3777.9, 3457.3, 3131.0, 2810.4, 2489.0, 2171.0, 1910.4, 1614.6, 1359.6, 1099.0, 838.3, 611.5, 411.1, 240.8, 131.6, 81.5]))

#prm_sigma_zones  =  len(prm_sigma_zones) * [0]

# il faut créer des input pour ces deux paramètres aller savoir pourquoi ...
prm_zones_bound_input = list(prm_zones_bound)
prm_sigma_zones_input = list(prm_sigma_zones)
    
### else:
prm_sigma_t = 2 * 10**-4
prm_sigma_t = 0

# K is "Mersenne Twister pseudo-random number generator"
prm_K_xyz = 1111

prm_K_t_hdwr = 5000
prm_K_t_zone = 9993

prm_K_t = 4000
    
prm_traject = 'cross'
prm_traject = 'derive'
prm_traject = 'droite'

#DROITE
prm_drt_xsize   = 4100
prm_drt_ysize   = 4100
prm_drt_xcenter = 10000
prm_drt_ycenter = 10000
prm_drt_zcenter = 0
prm_drt_nb_pass = 10
prm_drt_nb_obs  = 10
prm_drt_angle   = 90

prm_drt_xsize   = 500
prm_drt_ysize   = 500
prm_drt_xcenter = 10000
prm_drt_ycenter = 10000
prm_drt_zcenter = 0
prm_drt_nb_pass = 3
prm_drt_nb_obs  = 33
prm_drt_angle   = 0


# CROSS
prm_crs_xsize   = 0
prm_crs_ysize   = 0
prm_crs_xsize   = 2000
prm_crs_ysize   = 2000
prm_crs_xcenter = 10000
prm_crs_ycenter = 10000
prm_crs_zcenter = 0
prm_crs_nb_pass = 3
prm_crs_nb_obs  = 1000 
prm_crs_angle   = 20
    
prm_drv_x0        = 10000
prm_drv_y0        = 10000 # start of the point 
prm_drv_xcenter   = 10000
prm_drv_ycenter   = 10000 # center of the circle
prm_drv_nb_obs    = 10000
prm_drv_R         = 50
prm_drv_step_size = 10

prm_drv_param     = 1
prm_drv_K_derive  = 1400

plot        = 0
plot_in_fct = 0
procs       = 8 #procos
timestamp_in_filename = 0

prm_noise_on_datation_std  = 10**-6
prm_noise_on_datation_seed = 54321
# NB : doit encore être implémenté ...
    
prm_forwd_bkwrd_mode = True
    
with_export = 1

ping_export_file_type = 'P'
with_stock_SSF_n_C = False

# ===============================
# FILE MANAGEMENT
# ===============================


#midfix = ''
print(vari_dico)

for k,v in vari_dico.items():
    #locals()[k] = v
    globals()[k] = v
    
for k in vari_dico.keys():
    print(k)
    exec('print ' + k)

#### gestion des bruits
# il faut créer des input pour ces deux paramètres aller savoir pourquoi ...
prm_sigma_zones = prm_coef_sigma_zones * np.array(prm_sigma_zones_input)

if prm_null_water_column_noise: 
    prm_sigma_zones  =  len(prm_sigma_zones) * [0]
    prm_sigma_t_hdwr = 0        # erreur "hardware" résiduel propre à chaque PXP    

prm_sigma_zones = list(prm_sigma_zones)
prm_zones_bound = list(prm_zones_bound_input)

if prm_traject == 'droite':
    prm_nb_obs = prm_drt_nb_obs
elif prm_traject == 'derive': 
    prm_nb_obs = prm_drv_nb_obs
elif prm_traject == 'cross': 
    prm_nb_obs = prm_crs_nb_obs
    

if not with_export:
    print("WARN : export désactivé !!!!")
    print("activer variable with_export")

if imporved_noising:
    noise4exp = str(prm_sigma_t_hdwr)
else:
    noise4exp = str(prm_sigma_t)


if prm_with_temporal:
    vitstr = 'vit' + str(prm_vit)
else:
    vitstr = ''

if prm_with_ssf and not prm_ssf_realistik_grad:
    gradstr = 'grad' + str(prm_x_grad) + '_'  + str(prm_y_grad)
elif prm_with_ssf and prm_ssf_realistik_grad:
    gradstr = 'grad_realik'
else:
    gradstr = ""
    
if prm_traject == 'droite':
    exp = '_'.join((exp_prefix  , str(prm_drt_nb_pass) + 'x' +  str(prm_drt_nb_obs) ,
                    'x' + str(prm_drt_xsize) , 'y' + str(prm_drt_ysize) , 'ang' + str(prm_drt_angle) ,
                    "nois" + str(imporved_noising) + '-' + noise4exp , vitstr ,gradstr,
                     midfix , exp_suffix))
elif prm_traject == 'cross':
    exp = '_'.join((exp_prefix  , "cross" , str(prm_crs_nb_pass) + 'x' +  str(prm_crs_nb_obs) ,
                    'x' + str(prm_crs_xsize) , 'y' + str(prm_crs_ysize) , 
                    "nois" + str(imporved_noising) + '-' + noise4exp , vitstr ,gradstr,
                     midfix , exp_suffix))
elif prm_traject == 'derive':
    exp = '_'.join((exp_prefix  , str(prm_drv_nb_obs) ,
                    'R' + str(prm_drv_R) , 
                    "nois" + str(imporved_noising) + '-' + noise4exp , vitstr ,gradstr,
                    midfix , exp_suffix))    

if imporved_noising:
    print('NOISE : ' , prm_sigma_zones , prm_sigma_t_hdwr)

path_exp  = '/'.join((path_gene,exp))
prm_timestamp = genefun.get_timestamp()

if not os.path.exists(path_exp):
    os.makedirs(path_exp)

# ===============================
# IMPORT DU SSP 
# ===============================

# cas 1 :
# Import d'un SSP réaliste, matlab-like trouvé dans un coin 
SSP = np.loadtxt(prm_ssp_file_path)
Z = np.squeeze(SSP[:,0])
C = np.squeeze(SSP[:,1])

min_PXP_depth = np.min([pxp[-1] for pxp in PXP_lis])
if np.max(Z) < min_PXP_depth:
    Z,C = ssp.SSP_extrapolate(Z,C,min_PXP_depth + 500,1)

if prm_ssp_decimate:
    Z,C = ssp.SSP_light(Z,C)

# cas 2 : 
if prm_with_temporal:
    sspt_dic  = genefun.pickle_loader(prm_sspt_dic_path)
    print("Iterateur en cours de generation")
    I , Zt    = acls.sspdic2InterpoSSPT(sspt_dic,prm_tempor_start,prm_tempor_len)
    print("Iterateur généré")
    tempor_ctd_diff_tmp  = prm_tempor_ctd_time - prm_tempor_start 
    prm_tempor_ctd_epoch = tempor_ctd_diff_tmp.seconds + tempor_ctd_diff_tmp.days * 86400.
# cas 3 : le cas echeant : import du Gdic :
if prm_with_pseudo_ssf:
    Ggraddic = genefun.pickle_loader(prm_Gdic_path)
# ===============================
# FABRICATION DE LA TRAJECTOIRE
# ===============================

# Génération 
print("#### Generation trajectoire" , prm_traject)

traject_return_mode = 'interp'
    
if prm_traject == 'droite':
    XYZreturned , E , _ = acls.fabriq_traject_droite(prm_drt_xsize,prm_drt_ysize,prm_drt_xcenter,
                                             prm_drt_ycenter,prm_drt_zcenter,prm_drt_nb_pass,
                                             prm_drt_nb_obs,prm_drt_angle,prm_vit,
                                             epoch_init=prm_epoch_init,
                                             plot=plot_in_fct,
                                             noise_on_datation_seed = 0,
                                             noise_on_datation_std  = 0,
                                             return_mode=traject_return_mode)
elif prm_traject == 'derive':                          
    XYZreturned , E , circle = acls.fabriq_traject_derive(prm_drv_x0,prm_drv_y0,
                                         prm_drv_xcenter,prm_drv_ycenter,
                                         prm_drv_R,prm_drv_step_size,prm_drv_nb_obs,prm_vit,
                                         prm_epoch_init,prm_drv_K_derive,plot=plot_in_fct,
                                         noise_on_datation_seed = 0,
                                         noise_on_datation_std  = 0,
                                         return_mode=traject_return_mode)

elif prm_traject == 'cross':                          
    XYZreturned , E , _ = acls.fabriq_traject_cross(prm_crs_xsize,prm_crs_ysize,prm_crs_xcenter,
                                             prm_crs_ycenter,prm_crs_zcenter,prm_crs_nb_pass,
                                             prm_crs_nb_obs,prm_crs_angle,prm_vit,
                                             epoch_init=prm_epoch_init,
                                             plot=plot_in_fct,
                                             noise_on_datation_seed = 0,
                                             noise_on_datation_std  = 0,
                                             return_mode=traject_return_mode)
    




IXYZ          = XYZreturned
XYZ           = IXYZ(E).T
XYZ_emi       = np.array(XYZ)
XYZ_emi_clean = np.array(XYZ)

if prm_traject == 'derive':
    print(XYZ_emi_clean[0] , prm_drv_x0 , prm_drv_y0)


#N_noise_on_datation_seed  = np.random.RandomState(prm_noise_on_datation_seed)
#prm_noise_on_datation_std


# =======================================================
# FABRICATION DES PINGS (MODE CLASSIC UNIQUEMENT FORWARD)
# =======================================================                                         

Tlen     = XYZ.shape[0]
Tshape   = (Tlen,1)

if not prm_forwd_bkwrd_mode:
    None
#    # Noise commun pour chaque beacon   
#    if imporved_noising:  
#        N_t_zone = np.random.RandomState(prm_K_t_zone)
#        if mono_epoc:
#            Tnse_proto = np.tile(N_t_zone.randn(len(prm_sigma_zones)),Tlen).reshape((Tlen,len(prm_sigma_zones)))
#        else:
#            Tnse_proto = N_t_zone.randn(Tlen,len(prm_sigma_zones))
#        Tnse = np.multiply( Tnse_proto ,  prm_sigma_zones )
#    
#    Ctempo_stk_stk = []                  
#    for prm_ipxp,prm_pxp_coords in enumerate(PXP_lis): 
#        prm_ipxp = prm_ipxp + 1  # having a valid ID of PXP
#        print "pings for PXP", prm_ipxp
#        args_lis = []
#        pool = mp.Pool(processes=procs)
#        
#        # Generation
#        xyzlis,PXPlis,Zlis,Clis = [],[],[],[]
#        Ctempo_stk = []
#        
#        for i in range(XYZ.shape[0]):  
#            xyz = XYZ[i,:]
#            if prm_with_temporal:
#                Z = Zt
#                C = acls.SSPT_from_Interpo(I,Zt,E[i])
#                Ctempo_stk.append(C)
#            args_lis.append((xyz,prm_pxp_coords,Z,C,0,88,False,True))
#        
#        Ctempo_stk_stk.append(Ctempo_stk)
#    
#        rt.raytrace_seek(*args_lis[0])
#    
#        results = [pool.apply_async(rt.raytrace_seek, args=x) for x in args_lis] 
#        
#        # For debug
#        #for x in args_lis:
#        #    rt.raytrace_seek(*x)        
#    
#        A_stk       = [e.get()[0] for e in results]
#        R_stk       = [np.sum(e.get()[1]) for e in results]
#        T_stk       = [np.sum(e.get()[2]) for e in results]
#        R_full_stk  = [e.get()[1] for e in results]
#        T_full_stk  = [e.get()[2] for e in results]
#               
#        Dz_full_stk , S_full_stk = rt.find_dZ_for_fabrik_ping(Z,C,A_stk,XYZ,prm_pxp_coords[-1],True)
#    
#        T_clean  = np.expand_dims(T_stk,1)
#        A_clean  = np.expand_dims(A_stk,1)
#        R_clean  = np.expand_dims(R_stk,1)
#    
#        # Bruitage 
#        if imporved_noising:
#            # Hardware
#            prm_K_t_hdwr_pxp = prm_ipxp + prm_K_t_hdwr
#            N_t_hdwr  = np.random.RandomState(prm_K_t_hdwr_pxp)
#            Tnse_hdwr = N_t_hdwr.randn(Tlen) * prm_sigma_t_hdwr
#            Tnse_hdwr = np.expand_dims(Tnse_hdwr,1)
#                    
#            # Zone
#            Tnse_zone = rt.noise_in_zones(prm_zones_bound,Tnse,Dz_full_stk,T_full_stk).T
#            Tnse_zone_total = np.expand_dims(np.sum(Tnse_zone,axis=1),1)
#            
#            Tnse_total = Tnse_hdwr + Tnse_zone_total
#            T_noise  = T_clean + Tnse_total      
#        else:
#            prm_K_t_pxp = prm_ipxp + prm_K_t
#            N_t = np.random.RandomState(prm_K_t_pxp)
#            T_noise = T_clean + N_t.randn(*Tshape) * prm_sigma_t
            

# =======================================================
# FABRICATION DES PINGS MODE FORWARD - BACKWARD (160420)
# ======================================================= 

else:
    Ctempo_stk_stk   = []
    SSFtempo_stk_stk = []
    
    if imporved_noising:  
        N_t_zone = np.random.RandomState(prm_K_t_zone)
        if mono_epoc:
            Tnse_proto = np.tile(N_t_zone.randn(len(prm_sigma_zones)),Tlen).reshape((Tlen,len(prm_sigma_zones)))
        else:
            Tnse_proto = N_t_zone.randn(Tlen,len(prm_sigma_zones))
        Tnse = np.multiply( Tnse_proto ,  prm_sigma_zones )
        
    for prm_ipxp,prm_pxp_coords in enumerate(PXP_lis): 
        prm_ipxp = prm_ipxp + 1  # having a valid ID of PXP
        print("pings for PXP", prm_ipxp)
        args_lis = []
        pool = mp.Pool(processes=procs)
        
        #### Generation emission
        print("#### Generation emission, SD mode")
        xyzlis,PXPlis,Zlis,Clis = [],[],[],[]
        Ctempo_stk = []
        
        for i in range(XYZ.shape[0]):  
            xyz = XYZ[i,:]
            if prm_with_temporal:
                Z = Zt
                C = acls.SSPT_from_Interpo(I,Zt,E[i])

                if np.max(Z) < min_PXP_depth:
                    Z,C = ssp.SSP_extrapolate(Z,C,min_PXP_depth + 500,1)
                                        
                if np.any( np.logical_not( np.logical_and( 1400. < C  ,  C < 1600. ))):
                    print("ERR : NaN in the SSP")
                    raise Exception
                
                Ctempo_stk.append(C)
            args_lis.append((xyz,prm_pxp_coords,Z,C,0,88,False,True))
        
        if prm_nb_obs < 11000:
            print('INFO : prm_nb_obs < 11000' , i)
            Ctempo_stk_stk.append(Ctempo_stk)
        
        if 0: # BLOC DEBUGG
            rt.raytrace_seek(*args_lis[0])
        if not with_async:
            results = [pool.apply(rt.raytrace_seek, args=x) for x in args_lis] 
        else:
            results = [pool.apply_async(rt.raytrace_seek, args=x) for x in args_lis] 
            results = [e.get() for e in results]
        pool.close()

        T_emi_clean = np.array([np.sum(e[2]) for e in results])
        
        A_emi            = [e[0] for e in results]
        R_stk            = [np.sum(e[1]) for e in results]
        #R_full_stk      = [e.get()[1] for e in results]
        T_emi_clean_full = [e[2] for e in results]
        E_emi_clean      = E
        
        T_emi_clean_SD = T_emi_clean
        
        #### Generation reception
        print("#### Generation reception, SD mode")
        pool = mp.Pool(processes=procs)
        args_lis = []
        for iii in range(len(T_emi_clean)):                
            t_em   = T_emi_clean[iii]
            e_em   = E[iii] 
            if prm_with_temporal:
                Z = Zt
                C = Ctempo_stk[iii]

            args_lis.append(((E,XYZ) , prm_pxp_coords , e_em , t_em , Z , C))
        
        
                    
        if not with_async:
            results = [pool.apply(acls.t_rec_finder, args=x) for x in args_lis] 
        else:
            results = [pool.apply_async(acls.t_rec_finder, args=x) for x in args_lis] 
            results = [e.get() for e in results]
        

        pool.close()

        if 0: # BLOC DEBUGG
            t_rec , e_rec , xyz_rec , a_rec , T_rec_full = acls.t_rec_finder(IXYZ , prm_pxp_coords ,
                                                        e_em , t_em , Z , C )
                                                        

        A_rec            = [e[3] for e in results]
        T_rec_clean      = np.array([e[0] for e in results])
        T_rec_clean_SD   = T_rec_clean
        E_rec_clean      = np.array([e[1] for e in results])
        XYZ_rec_clean    = np.array([e[2] for e in results])
        XYZ_rec          = np.array(XYZ_rec_clean)
        T_rec_clean_full = [e[4] for e in results]    
        
        ############## SSF MODE ##################

        if not prm_SSF_based_on_interpolation:
            Znorm , Cnorm = ssp.munk_pseudo_coef_calc(Z,C,norma_step=.5)
            SSPobjt = acls.SSP(Znorm,Cnorm,e=0)
            ssf_real = acls.make_SSF3D_from_SSP(SSPobjt,
                                           xmin=np.min(PXP_arr[:,0]) - 1000,
                                           xmax=np.max(PXP_arr[:,0]) + 1000,
                                           ymin=np.min(PXP_arr[:,1]) - 1000,
                                           ymax=np.max(PXP_arr[:,1]) + 1000)
                   
            zm,cm      = ssp.munk(6000,.1)
            SSPmunk    = acls.SSP(zm,cm,e=0)
            ssf_munk   = acls.make_SSF3D_from_SSP(SSPmunk,
                                           xmin=np.min(PXP_arr[:,0]) - 1000,
                                           xmax=np.max(PXP_arr[:,0]) + 1000,
                                           ymin=np.min(PXP_arr[:,1]) - 1000,
                                           ymax=np.max(PXP_arr[:,1]) + 1000)

            SSF = ssf_real 

        else:
            
            Udic = gf.pickle_loader(prm_Udic_path)
            X4I,Z4I,CgridXZ4I = Udic['X'] , Udic['Z_X'] ,  Udic['CgridXZ']
            SSF_a =  acls.make_SSF3D_from_SSP_n_distance((X4I,Z4I,CgridXZ4I),'X',
                                           xmin=np.min(PXP_arr[:,0]) - 1000,
                                           xmax=np.max(PXP_arr[:,0]) + 1000,
                                           ymin=np.min(PXP_arr[:,1]) - 1000,
                                           ymax=np.max(PXP_arr[:,1]) + 1000)
                                           


            SSF_b =  acls.make_SSF3D_from_SSP_n_distance((X4I,Z4I,CgridXZ4I),'Y',
                                           xmin=np.min(PXP_arr[:,0]) - 1000,
                                           xmax=np.max(PXP_arr[:,0]) + 1000,
                                           ymin=np.min(PXP_arr[:,1]) - 1000,
                                           ymax=np.max(PXP_arr[:,1]) + 1000)
                                           


            SSF =  SSF_a
            
            
            
            
                                                                                          
            

        if prm_ssf_realistik_grad:
            Karr          = np.loadtxt(prm_Kfile_path)
            Tensor_grad_x = SSF.add_a_gradient_from_vectors(0,Karr)
            #Tensor_grad_y = SSF.add_a_gradient_from_vectors(1,Karr)
        else:
            if prm_x_grad != 0:
                Tensor_grad_x = SSF.add_a_gradient(0,prm_x_grad,z_max_smooth=prm_z_mixed)
            if prm_y_grad != 0:
                Tensor_grad_y = SSF.add_a_gradient(1,prm_y_grad,z_max_smooth=prm_z_mixed)
                                                             
        #### Generation emission
        print("#### Generation emission, SSF mode")
        xyzlis,PXPlis,Zlis,Clis = [],[],[],[]
        SSFtempo_stk = []
        Ctempo_stk   = []
        
        for i in range(XYZ.shape[0]):  
            xyz = XYZ[i,:]
            if prm_with_temporal:
                Z = Zt
                C = acls.SSPT_from_Interpo(I,Zt,E[i])

                if np.any(np.logical_not(np.logical_and( 1400. < C ,
                                                         C < 1600. ))):
                    print("ERR : NaN in the SSP")
                    raise Exception
                    
                Ctempo_stk.append(C)
                
                Znorm , Cnorm = ssp.munk_pseudo_coef_calc(Z,C,
                                                          norma_step=.5)
                SSPobjt = acls.SSP(Znorm,Cnorm,e=0)
                ssf_real = acls.make_SSF3D_from_SSP(SSPobjt,
                                               xmin=np.min(PXP_arr[:,0]) - 1000,
                                               xmax=np.max(PXP_arr[:,0]) + 1000,
                                               ymin=np.min(PXP_arr[:,1]) - 1000,
                                               ymax=np.max(PXP_arr[:,1]) + 1000)
                                               
                SSF = ssf_real   

                if with_stock_SSF_n_C:
                    SSFtempo_stk.append(SSF)
                    SSFtempo_stk = [SSF] * XYZ.shape[0]
            ang_s_apri = acls.canonical_shooting_angle(xyz,prm_pxp_coords)
            args_lis.append((xyz,prm_pxp_coords,SSF,ang_s_apri,
                             prm_eiko_h,prm_eiko_resotype))
        
        if with_stock_SSF_n_C:
            Ctempo_stk_stk.append(SSFtempo_stk)


        if 0: #DEBUG CANONICAL ANGLE
            for pxp in PXP_lis: #canonical_shooting_angle
                for xyz in XYZ:
                    CSA = acls.canonical_shooting_angle(xyz,pxp)
                    print(CSA[0] , CSA[1])
                    print(np.deg2rad(CSA[0]) , np.deg2rad(CSA[1]))
                
        if 0: # BLOC DEBUGG
            xyz = np.array([10000,10000,0])

            ang_s_apri     = acls.canonical_shooting_angle(xyz,PXP4)
            ang_s_apri     = acls.canonical_shooting_angle(xyz,prm_pxp_coords)
            results_eiko   = acls.raytrace_ODE_seek(xyz,prm_pxp_coords,
                                                  SSF,ang_s_apri,
                                                  prm_eiko_h,
                                                  prm_eiko_resotype)

        
        if not with_async:
            results = [pool.apply(acls.raytrace_ODE_seek,args=x) for x in args_lis] 
        else:
            results = [pool.apply_async(acls.raytrace_ODE_seek,args=x) for x in args_lis] 
            results = [e.get() for e in results]
        

        pool.close()
        
        pool = mp.Pool(processes=procs)
        args_lis_direct = []
        for i,e in enumerate(results):
            xyz = XYZ[i,:]
            args_lis_direct.append((SSF,xyz,e[0][0:2],e[-1][1],prm_eiko_h,prm_eiko_resotype))
        if not with_async:
            results_directs = [pool.apply(acls.raytrace_ODE_2or3d_new,args=x) for x in args_lis_direct]
        else:
            results_directs = [pool.apply_async(acls.raytrace_ODE_2or3d_new,args=x) for x in args_lis_direct]
            #results_directs = [e.get() for e in results_directs]
            results_directs_2 = []
            for e in results_directs:
                try:
                    results_directs_2.append(e.get())
                except:
                    print("WARN : something went wrong in results_directs, we use the last working result ...")
                    results_directs_2.append(results_directs_2[-1])
            results_directs = results_directs_2
                    
            
        pool.close()

        T_emi_clean =  np.array([e[-1][-1]       for e in results])
        A_emi       =  np.array([np.abs(e[0][0]) for e in results])
        R_stk       =  np.array([e[-1][1]        for e in results])
        #R_full_stk = [e.get()[1] for e in results]
        E_emi_clean      = E
        T_emi_clean_full = []
        for e in results_directs:
            c, s, z = e[2] , e[1] , e[0][:,2]
            T_emi_clean_full.append(acls.t_cumul_from_c_s_out(e[2],e[1]))   

        #### Generation reception
        if not prm_bypass_ssf_rec:
            print("#### Generation reception, SSF mode (not bypassed)")

            pool = mp.Pool(processes=procs)
            kwargs_lis     = []
            args_lis       = []
            SSFtempo_stk   = []
            for iii in range(len(T_emi_clean)):                
                t_em   = T_emi_clean[iii]
                e_em   = E[iii] 
                if prm_with_temporal:
                    Z = Zt
                    C = acls.SSPT_from_Interpo(I,Zt,E[i])

                    if np.any(np.logical_not(np.logical_and( 1400. < C ,
                                                             C < 1600. ))):
                        print("ERR : NaN in the SSP")
                        raise Exception
                    
                    Znorm , Cnorm = ssp.munk_pseudo_coef_calc(Z,C,norma_step=.1)
                    SSPobjt = acls.SSP(Znorm,Cnorm,e=0)
                    ssf_real = acls.make_SSF3D_from_SSP(SSPobjt,
                                                   xmin=np.min(PXP_arr[:,0]) - 1000,
                                                   xmax=np.max(PXP_arr[:,0]) + 1000,
                                                   ymin=np.min(PXP_arr[:,1]) - 1000,
                                                   ymax=np.max(PXP_arr[:,1]) + 1000)
                                                   
                    SSF = ssf_real
                    Tensor_grad_x = SSF.add_a_gradient(0,prm_x_grad,z_max_smooth=prm_z_mixed)
                    Tensor_grad_y = SSF.add_a_gradient(1,prm_y_grad,z_max_smooth=prm_z_mixed)

                    #AAASSFtempo_stk.append(SSF)
    
                kwargs_lis.append({'InterpXYZ' : (E,XYZ) ,
                                  'XYZpxp'    : prm_pxp_coords ,
                                  'E_em_inp'  : e_em  ,
                                  't_em_inp'  : t_em  ,
                                  'Z' : []  , 'C' : []  ,
                                  'SDmode'  : False , 'SSF' : SSF,
                                  'theta_apri' : - A_emi[iii] , 
                                  's_apri' : R_stk[iii] })
                                
                args_lis.append(((E,XYZ) , prm_pxp_coords , e_em  , t_em  ,
                                 [],  []  , 0 , 's' , False ,  SSF , - A_emi[iii] , None , R_stk[iii]))

            if not with_async:
                results = [pool.apply(acls.t_rec_finder, tuple() , x) for x in kwargs_lis] 
            else:
                results = [pool.apply_async(acls.t_rec_finder, tuple() , x) for x in kwargs_lis] 
                results = [e.get() for e in results]
            
            A_rec            = [e[3] for e in results]
            T_rec_clean      = np.array([e[0] for e in results])
            E_rec_clean      = np.array([e[1] for e in results])
            XYZ_rec_clean    = np.array([e[2] for e in results])
            XYZ_rec          = np.array(XYZ_rec_clean)
            T_rec_clean_full = [e[4] for e in results]       
            T_rec_clean_SSF  = T_rec_clean
        
        else:
            print("#### Generation reception, SSF mode, BYPASSED !")
            A_rec       = A_emi
            T_rec_clean = T_emi_clean
            E_rec_clean = E_emi_clean
            XYZ_rec_clean = XYZ_emi_clean
            XYZ_rec       = np.array(XYZ_rec_clean)                
            T_rec_clean_full = T_emi_clean_full 
            


        

        

print("name of exp.")
print(exp)
print("path of exp.")
print(path_exp)

print(" ==== FIN de la FABRICATION ==== ")

if plot:
    fig , ax = plt.subplots()
    if prm_traject == 'derive':
        ax.plot(circle[:,0],circle[:,1],'gx')
    for pxp in PXP_lis:
        ax.scatter(pxp[0],pxp[1],color='r',s=300)
    ax.plot(XYZ[:,0],XYZ[:,1],'b.')
    plt.axis('equal')
    

if 1:
    print("Si bloqué, verifier que le lancement de la fct n'est pas activé dans le script core")
    outup  = fabrik_fct(exp_prefix)
#    A, B = outup
#        
#BC = B.C - B.C[:,0,:]
#
#A.C + BC 
#
#
#plt.plot(A.C[1,1,:])