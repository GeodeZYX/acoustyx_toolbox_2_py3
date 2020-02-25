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

@author: psakicki
"""

import glob
import acouclass as acls
import numpy as np
import raytrace as rt
from megalib import *
import matplotlib.pyplot as plt

# ===============================
# PARAMETERS
# ===============================

if platform.node() == 'calipso':
    prm_ssp_file_path = '/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/exemple/SSP/SSP_NOAA_dep5645_20030608000000'
    path_gene = '/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/working'
    sspt_dic_path = "/home/psakicki/THESE/bd_CTD/NEW_DATA_2015/MOVE/MOVE_reform_1507/SIMPLE/sspdic_sensor.pik"
elif platform.node() == 'pierre-MS-16Y1':
    prm_ssp_file_path = '/home/pierre/Documents/CODES/acoustyx_toolbox_2/exemple/SSP/SSP_NOAA_dep5645_20030608000000'
    path_gene = '/home/pierre/Documents/CODES/acoustyx_toolbox_2/working'
    sspt_dic_path = "/media/pierre/D34F-24E3/plots/SSPs/sspdic"

exp_prefix  = 'tpo1' # for temporal
exp_prefix  = 'tst1'
exp_prefix  = 'repris1'

exp_prefix  = 'batc2b_deriv'
exp_prefix  = 'batc2d_deriv_megadecalZ'
exp_prefix  = 'batc3b_2000_bw_deriv_decalZ'
exp_prefix  = 'batc2c_deriv_decalZ_notrajnoise_stdpingnoise'
exp_prefix  = 'batc2c_deriv_decalZ_nopingnoise_trajnoise5decim'
exp_prefix  = 'batc4_nadir'

exp_prefix  = 'test_repriz_adel_FB0_mk3'
exp_prefix  = 'ADEL33b'

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

ddd  = 4000 
PXP1 = np.array([-ddd,-ddd,4010])
PXP2 = np.array([ ddd, ddd,4005])
PXP3 = np.array([-ddd, ddd,3998])
PXP4 = np.array([ ddd,-ddd,3985])

PXP1 = np.array([-2500,-2500,4000])
PXP2 = np.array([2500,2500,4000])
PXP3 = np.array([-2500,2500,4000])
PXP4 = np.array([2500,-2500,4000])


#PXP1 = np.array([0,0,4000])
PXP_lis = [PXP1]

PXP_lis = [PXP1,PXP2,PXP3,PXP4]


#PXP_lis = [PXP1]

PXP_lis = [pxp + np.array([10000,10000,0]) for pxp in PXP_lis]

# TEMPORAL
prm_with_temporal   = False
prm_tempor_start    = dt.datetime(2008, 1, 4, 19, 20)
prm_tempor_len      = 86400*7
prm_tempor_ctd_time = dt.datetime(2008, 1, 5, 19, 20)
prm_vit             = 500 # en m/h
prm_epoch_init      = 0 

# NOISE
prm_sigma_x = 0 
prm_sigma_y = 0 
prm_sigma_z = 0
prm_sigma_t = 0

prm_sigma_x = 0.001 
prm_sigma_y = 0.001 
prm_sigma_z = 0.002
prm_sigma_t = 0

prm_sigma_x =  0.10
prm_sigma_y =  0.10
prm_sigma_z =  0.15

prm_sigma_x = 0.03
prm_sigma_y = 0.03 
prm_sigma_z = 0.05


prm_sigma_x =  0.50
prm_sigma_y =  0.50
prm_sigma_z =  1.

prm_sigma_x =  0.0
prm_sigma_y =  0.0
prm_sigma_z =  0.0

imporved_noising = False

if imporved_noising:
    mono_epoc = False
    prm_zones_bound  = [100,500,1500,2500]
    prm_sigma_zones  = [10**-2,10**-2,10**-3,10**-3,10**-4] 
    prm_sigma_zones  = [10**-4,10**-4,10**-4,10**-4,10**-4] 
    prm_sigma_t_hdwr = 0 # erreur "hardware" résiduel propre à chaque PXP
    prm_sigma_t_hdwr = 1*10**-6 # erreur "hardware" résiduel propre à chaque PXP
else:
    prm_sigma_t = 2 * 10**-4
    prm_sigma_t = 0

# K is "Mersenne Twister pseudo-random number generator"
prm_K_xyz = 1111
if imporved_noising:
    prm_K_t_hdwr = 5000
    prm_K_t_zone = 5100
else:
    prm_K_t = 4000

prm_traject = 'droite'
prm_traject = 'derive'

if prm_traject == 'droite':
    prm_xsize   = 00
    prm_ysize   = 00
    prm_xsize   = 2000
    prm_ysize   = 2000
    prm_xcenter = 10000
    prm_ycenter = 10000
    prm_zcenter = 0
    prm_nb_pass = 3
    prm_nb_obs  = 50
    prm_angle   = 45
    
elif prm_traject == 'derive':
    prm_x0        = 10000
    prm_y0        = 10000
    prm_xcenter   = 10000
    prm_ycenter   = 10000
    prm_nb_obs    = 5000
    prm_R         = 50
    prm_step_size = 1
    prm_param     = 1
    prm_K_derive  = 1123
    prm_K_derive  = 1122

plot = 0
plot_in_fct = 0
procs = 7
timestamp_in_filename = 0

with_export = 1

# ===============================
# FILE MANAGEMENT
# ===============================

if not with_export:
    print("WARN : export désactivé !!!!")
    print("activer a True with_export")

if imporved_noising:
    noise4exp = str(prm_sigma_t_hdwr)
else:
    noise4exp = str(prm_sigma_t)

if prm_with_temporal:
    vitstr = 'vit' + str(prm_vit)
else:
    vitstr = ''
    
    
if prm_traject == 'droite':
    exp = '_'.join((exp_prefix  , str(prm_nb_pass) + 'x' +  str(prm_nb_obs) ,
                    'x' + str(prm_xsize) , 'y' + str(prm_ysize) , 
                    "nois" + str(imporved_noising) + '-' + noise4exp , vitstr ,
                    exp_suffix))
elif prm_traject == 'derive':
    exp = '_'.join((exp_prefix  , str(prm_nb_obs) ,
                    'R' + str(prm_R) , 
                    "nois" + str(imporved_noising) + '-' + noise4exp , vitstr ,
                    exp_suffix))    

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

# cas 2 : 
if prm_with_temporal:
    sspt_dic = genefun.pickle_loader(sspt_dic_path)
    I , Zt   = acls.sspdic2InterpoSSPT(sspt_dic,prm_tempor_start,prm_tempor_len)
    prm_tempor_ctd_epoch = (prm_tempor_ctd_time - prm_tempor_start).seconds / 3600.

# ===============================
# FABRICATION DE LA TRAJECTOIRE
# ===============================

# Génération 

if prm_traject == 'droite':
    XYZ , E , _ = acls.fabriq_traject_droite(prm_xsize,prm_ysize,prm_xcenter,
                                             prm_ycenter,prm_zcenter,prm_nb_pass,
                                             prm_nb_obs,prm_angle,prm_vit,
                                             epoch_init=prm_epoch_init,plot=plot_in_fct)

elif prm_traject == 'derive':                          
    XYZ , E , circle = acls.fabriq_traject_derive(prm_x0,prm_y0,prm_xcenter,prm_ycenter,
                                         prm_R,prm_step_size,prm_nb_obs,prm_vit,
                                         prm_epoch_init,plot=plot_in_fct)


# Bruitage 
N_xyz = np.random.RandomState(prm_K_xyz)
XYZ_noise = XYZ + N_xyz.randn(*XYZ.shape) * \
np.repeat([[prm_sigma_x,prm_sigma_y,prm_sigma_z]],XYZ.shape[0],axis=0)


# ===============================
# FABRICATION DES PINGS
# ===============================                                         
     
# Noise commun pour chaque beacon
Tlen     = np.max(XYZ.shape)
Tshape   = (Tlen,1)

if imporved_noising:  
    N_t_zone = np.random.RandomState(prm_K_t_zone)
    if mono_epoc:
        Tnse_proto = np.tile(N_t_zone.randn(len(prm_sigma_zones)),Tlen).reshape((Tlen,len(prm_sigma_zones)))
    else:
        Tnse_proto = N_t_zone.randn(Tlen,len(prm_sigma_zones))
    Tnse = np.multiply( Tnse_proto ,  prm_sigma_zones )
    
 
Ctempo_stk_stk = []                  
for prm_ipxp,prm_pxp_coords in enumerate(PXP_lis): 
    prm_ipxp = prm_ipxp + 1  # having a valid ID of PXP
    print("pings for PXP", prm_ipxp)
    args_lis = []
    pool = mp.Pool(processes=procs)
    
    # Generation
    xyzlis,PXPlis,Zlis,Clis = [],[],[],[]
    Ctempo_stk = []
    
    for i in range(XYZ.shape[0]):  
        xyz = XYZ[i,:]
        if prm_with_temporal:
            Z = Zt
            C = acls.SSPT_from_Interpo(I,Zt,E[i])
            Ctempo_stk.append(C)
        args_lis.append((xyz,prm_pxp_coords,Z,C,0,88,False,True))
    
    Ctempo_stk_stk.append(Ctempo_stk)

    rt.raytrace_seek(*args_lis[0])

    results = [pool.apply_async(rt.raytrace_seek, args=x) for x in args_lis] 
    
    # For debug
    #for x in args_lis:
    #    rt.raytrace_seek(*x)        


    A_stk = [e.get()[0] for e in results]
    R_stk = [np.sum(e.get()[1]) for e in results]
    T_stk = [np.sum(e.get()[2]) for e in results]
    R_full_stk = [e.get()[1] for e in results]
    T_full_stk = [e.get()[2] for e in results]
           
    Dz_full_stk = rt.find_dZ_for_fabrik_ping(Z,C,A_stk,XYZ,prm_pxp_coords[-1])

    T_clean  = np.expand_dims(T_stk,1)
    A_clean  = np.expand_dims(A_stk,1)
    R_clean  = np.expand_dims(R_stk,1)

    # Bruitage 
    if imporved_noising:
        # Hardware
        prm_K_t_hdwr_pxp = prm_ipxp + prm_K_t_hdwr
        N_t_hdwr  = np.random.RandomState(prm_K_t_hdwr_pxp)
        Tnse_hdwr = N_t_hdwr.randn(Tlen) * prm_sigma_t_hdwr
        Tnse_hdwr = np.expand_dims(Tnse_hdwr,1)
                
        # Zone
        Tnse_zone = rt.noise_in_zones(prm_zones_bound,Tnse,Dz_full_stk,T_full_stk).T
        Tnse_zone_total = np.expand_dims(np.sum(Tnse_zone,axis=1),1)
        
        Tnse_total = Tnse_hdwr + Tnse_zone_total
        T_noise  = T_clean + Tnse_total      
    else:
        prm_K_t_pxp = prm_ipxp + prm_K_t
        N_t = np.random.RandomState(prm_K_t_pxp)
        T_noise = T_clean + N_t.randn(*Tshape) * prm_sigma_t

    # recherche des "prm_" 
    loctemp = dict(locals())
    outlis = []
    for k,v in loctemp.items():
        if 'prm_' in k:
            outlis.append(str(k[4:]) + ' : ' + str(v))
    outlis.sort()
    header = '\n'.join(outlis)
            
    # Enregistrement du Mfile
    if len(PXP_lis) != 1 :
        idpxp = 'PXP' + str(prm_ipxp)
    else:
        idpxp = ''
    
    if timestamp_in_filename:    
        ts = prm_timestamp
    else:
        ts = ''
    if not imporved_noising:
        header = header + '\nfields : XYZ_noise,T_noise,XYZ,T_clean,R_clean,A_clean'
        Mout   = np.hstack((XYZ_noise,T_noise,XYZ,T_clean,R_clean,A_clean))
    else:        
        header = header + '\nfields : XYZ_noise,T_noise,XYZ,T_clean,R_clean,A_clean,Tnse_zone,Tnse_zone_total,Tnse_hdwr,Tnse_total'
        Mout   = np.hstack((XYZ_noise,T_noise,XYZ,T_clean,R_clean,A_clean,Tnse_zone,Tnse_zone_total,Tnse_hdwr,Tnse_total))

    if with_export:
        Mpath = os.path.join(path_exp,'_'.join((exp,ts))+'.'+idpxp+'.M.dat')
        np.savetxt(Mpath,Mout,header=header,fmt='%21.12f')
    
# ===============================
# FABRICATION & EXPORT DES BL
# ===============================
BL_clean = acls.BL_from_PXPlist(PXP_lis)
if with_export:
    Bpath = os.path.join(path_exp,'_'.join((exp,ts+'.B.dat')))
    np.savetxt(Bpath,BL_clean)

# ===============================
# EXPORT DU SSP
# ===============================
Zpath = os.path.join(path_exp,'_'.join((exp,ts+'.Z.dat')))
Cpath = os.path.join(path_exp,'_'.join((exp,ts+'.C.dat')))


if prm_with_temporal:
    Z = Zt
    C = acls.SSPT_from_Interpo(I,Zt,prm_tempor_ctd_epoch)
    header = 'ssp_original : ' + sspt_dic_path + ' @ ' + str(prm_tempor_ctd_time)
else:
    header = 'ssp_original : ' + prm_ssp_file_path
    
if with_export:
    np.savetxt(Zpath,Z,header=header)
    np.savetxt(Cpath,C,header=header)

print("name of exp.")
print(exp)

print(" ==== FIN de la FABRICATION ==== ")


if plot:
    fig , ax = plt.subplots()
    if prm_traject == 'derive':
        ax.plot(circle[:,0],circle[:,1],'gx')
    for pxp in PXP_lis:
        ax.scatter(pxp[0],pxp[1],color='r',s=300)
    ax.plot(XYZ[:,0],XYZ[:,1],'b.')
    plt.axis('equal')
