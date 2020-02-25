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
elif platform.node() == 'psakicki-MS-16Y1':
    prm_ssp_file_path = '/home/psakicki/Documents/CODES/acoustyx_toolbox_2/exemple/SSP/SSP_NOAA_dep5645_20030608000000'
    path_gene = '/home/psakicki/Documents/CODES/acoustyx_toolbox_2/working'
    sspt_dic_path = "/media/psakicki/D34F-24E3/plots/SSPs/sspdic"
    sspt_dic_path = "/home/psakicki/THESE/bd_CTD/NEW_DATA_2015/MOVE/MOVE_reform_1507/SIMPLE/sspdic_sensor.pik"

exp_prefix  = 'tpo1' # for temporal
exp_prefix  = 'tst1'
exp_prefix  = 'repris1'

exp_prefix  = 'batc2b_deriv'
exp_prefix  = 'batc2d_deriv_megadecalZ'
exp_prefix  = 'batc3b_2000_bw_deriv_decalZ'
exp_prefix  = 'batc2c_deriv_decalZ_notrajnoise_stdpingnoise'
exp_prefix  = 'batc2c_deriv_decalZ_nopingnoise_trajnoise5decim'
exp_prefix  = 'batc4_nadir'


exp_prefix  = 'test_repriz_adel'


exp_suffix  = ''

# PXP COORDINATES
PXP1 = np.array([-2500,-2500,4000])
PXP2 = np.array([2500,2500,4000])
PXP3 = np.array([-2500,2500,4000])
PXP4 = np.array([2500,-2500,4000])

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


PXP1 = np.array([0,0,4000])
PXP_lis = [PXP1]

PXP_lis = [PXP1,PXP2,PXP3,PXP4]


tat = 0

#PXP_lis = [PXP1]

PXP_lis = [pxp + np.array([10000,10000,0]) for pxp in PXP_lis]

# TEMPORAL
prm_with_temporal   = False
prm_tempor_start    = dt.datetime(2008, 1, 4, 19, 20)
prm_tempor_len      = 86400*7
prm_tempor_ctd_time = dt.datetime(2008, 1, 5, 19, 20)
prm_vit             = 500. / 3600. # en m/s
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

prm_sigma_x =  0.0
prm_sigma_y =  0.0
prm_sigma_z =  0.0

prm_sigma_x =  0.50
prm_sigma_y =  0.50
prm_sigma_z =  1.

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

prm_traject = 'derive'
prm_traject = 'droite'

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

prm_noise_on_datation_std  = 10**-6
prm_noise_on_datation_seed = 54321
# NB : doit encore être implémenté ...
    
prm_forwd_bkwrd_mode = True
    
with_export = 1

ping_export_file_type = 'P'



# ===============================
# FILE MANAGEMENT
# ===============================

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

traject_return_mode = 'interp'
    
if prm_traject == 'droite':
    XYZreturned , E , _ = acls.fabriq_traject_droite(prm_xsize,prm_ysize,prm_xcenter,
                                             prm_ycenter,prm_zcenter,prm_nb_pass,
                                             prm_nb_obs,prm_angle,prm_vit,
                                             epoch_init=prm_epoch_init,
                                             plot=plot_in_fct,
                                             noise_on_datation_seed = prm_noise_on_datation_seed,
                                             noise_on_datation_std  = prm_noise_on_datation_std,
                                             return_mode=traject_return_mode)
elif prm_traject == 'derive':                          
    XYZreturned , E , circle = acls.fabriq_traject_derive(prm_x0,prm_y0,prm_xcenter,prm_ycenter,
                                         prm_R,prm_step_size,prm_nb_obs,prm_vit,
                                         prm_epoch_init,prm_K_derive,plot=plot_in_fct,
                                         noise_on_datation_seed = prm_noise_on_datation_seed,
                                         noise_on_datation_std  = prm_noise_on_datation_std,
                                         return_mode=traject_return_mode)


IXYZ = XYZreturned
XYZ  = IXYZ(E).T
    

# Bruitage 
N_xyz = np.random.RandomState(prm_K_xyz)
XYZ_noise = XYZ + N_xyz.randn(*XYZ.shape) * \
np.repeat([[prm_sigma_x,prm_sigma_y,prm_sigma_z]],XYZ.shape[0],axis=0)
XYZ_emi_noise = XYZ_noise

# =======================================================
# FABRICATION DES PINGS (MODE CLASSIC UNIQUEMENT FORWARD)
# =======================================================                                         


Tlen     = np.max(XYZ.shape)
Tshape   = (Tlen,1)

if not prm_forwd_bkwrd_mode:
    # Noise commun pour chaque beacon   
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
    
        A_stk       = [e.get()[0] for e in results]
        R_stk       = [np.sum(e.get()[1]) for e in results]
        T_stk       = [np.sum(e.get()[2]) for e in results]
        R_full_stk  = [e.get()[1] for e in results]
        T_full_stk  = [e.get()[2] for e in results]
               
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
            

# =======================================================
# FABRICATION DES PINGS MODE FORWARD - BACKWARD (160420)
# ======================================================= 

else:
    Ctempo_stk_stk = []
    for prm_ipxp,prm_pxp_coords in enumerate(PXP_lis): 
        prm_ipxp = prm_ipxp + 1  # having a valid ID of PXP
        print("pings for PXP", prm_ipxp)
        args_lis = []
        pool = mp.Pool(processes=procs)
        
        # Generation emission
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
        
        if 0:
            rt.raytrace_seek(*args_lis[0])
        
        results = [pool.apply_async(rt.raytrace_seek, args=x) for x in args_lis] 
        pool.close()

        T_emi_clean = np.array([np.sum(e.get()[2]) for e in results])
        
        A_stk      = [e.get()[0] for e in results]
        R_stk      = [np.sum(e.get()[1]) for e in results]
        #R_full_stk = [e.get()[1] for e in results]
        #T_full_stk = [e.get()[2] for e in results]
        
        E_emi_clean   = E
        XYZ_emi_clean = XYZ
        
        pool = mp.Pool(processes=procs)
        # Generation reception
        args_lis = []
        for iii in range(len(T_emi_clean)):                
            t_em   = T_emi_clean[iii]
            e_em   = E[iii] 
            if prm_with_temporal:
                Z = Zt
                C = Ctempo_stk[iii]

            args_lis.append(((E,XYZ) , prm_pxp_coords , e_em , t_em , Z , C))

        results = [pool.apply_async(acls.t_rec_finder, args=x) for x in args_lis] 
        pool.close()

        if 1:
            t_rec , e_rec , xyz_rec = acls.t_rec_finder(IXYZ , prm_pxp_coords ,
                                                        e_em , t_em , Z , C )
                                                        
        T_rec_clean   = np.array([e.get()[0] for e in results])
        E_rec_clean   = np.array([e.get()[1] for e in results])
        XYZ_rec_clean = np.array([e.get()[2] for e in results])
        
        # Bruitage
        
        XYZ_rec_noise = XYZ_rec_clean + N_xyz.randn(*XYZ_rec_clean.shape) * \
        np.repeat([[prm_sigma_x,prm_sigma_y,prm_sigma_z]],XYZ_rec_clean.shape[0],axis=0)
        
        Ttmp_stk = []

        if imporved_noising:
            prm_K_t_hdwr_pxp = prm_ipxp + prm_K_t_hdwr
            N_t_hdwr         = np.random.RandomState(prm_K_t_hdwr_pxp)
        
        prm_K_t_pxp = prm_ipxp + prm_K_t
        N_t         = np.random.RandomState(prm_K_t_pxp)

        for T_clean in (T_emi_clean , T_rec_clean):
            if imporved_noising:                
                # Hardware
                Tnse_hdwr = N_t_hdwr.randn(Tlen) * prm_sigma_t_hdwr
                Tnse_hdwr = np.expand_dims(Tnse_hdwr,1)
                        
                # Zone
                Tnse_zone = rt.noise_in_zones(prm_zones_bound,Tnse,Dz_full_stk,T_full_stk).T
                Tnse_zone_total = np.expand_dims(np.sum(Tnse_zone,axis=1),1)
                
                Tnse_total = Tnse_hdwr + Tnse_zone_total
                Ttmp_noise  = T_clean + Tnse_total      
            else:
                Ttmp_noise  = T_clean + np.squeeze(N_t.randn(*Tshape)) * prm_sigma_t 

            Ttmp_stk.append(Ttmp_noise)
        
        TWTT_clean  = T_emi_clean + T_rec_clean + tat
        T_emi_noise = Ttmp_stk[0]
        T_rec_noise = Ttmp_stk[1]
        TWTT_noise  = T_emi_noise + T_rec_noise + tat
        TAT = np.array([tat] * len(TWTT_noise))
        
        # recherche des "prm_" 
        loctemp = dict(locals())
        outlis = []
        for k,v in loctemp.items():
            if 'prm_' in k:
                outlis.append(str(k[4:]) + ' : ' + str(v))
        outlis.sort()
        header = '\n'.join(outlis)
    
        # EXPORT DES PINGs
        # Enregistrement du Mfile
        if 'M' in ping_export_file_type:
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
    
        # le N file est discontinué car incomplet
    
        if 'O' in ping_export_file_type:
            if len(PXP_lis) != 1 :
                idpxp = 'PXP' + str(prm_ipxp)
            else:
                idpxp = ''
            
            if timestamp_in_filename:    
                ts = prm_timestamp
            else:
                ts = ''
    
    
            header = header + '\nfields : Epoch_emi, XYZ_emi ,ymdhms_emi, Epoch_rec , XYZ_rec , Ping_Raw , TAT , Ping_Ope_(minus_TAT_2_divided) '
            Mout   = np.hstack((XYZ_noise,T_noise,XYZ,T_clean,R_clean,A_clean,Tnse_zone,Tnse_zone_total,Tnse_hdwr,Tnse_total))
        
            if with_export:
                Mpath = os.path.join(path_exp,'_'.join((exp,ts))+'.'+idpxp+'.M.dat')
                np.savetxt(Mpath,Mout,header=header,fmt='%21.12f')
        
        
        # Le Pfile est le nouveau fichier, voulu flexible et définitif :
        # les titres des colonnes sont TOUTES explicitement nommées 
        # pour être ensuite interpétées par Pandas
        
        if 'P' in ping_export_file_type:
            if len(PXP_lis) != 1 :
                idpxp = 'PXP' + str(prm_ipxp)
            else:
                idpxp = ''
            
            if timestamp_in_filename:    
                ts = prm_timestamp
            else:
                ts = ''
                
                
                
            Datalis  = [E_emi_clean   ,
                        XYZ_emi_noise ,
                        T_emi_noise   ,     
                        E_rec_clean   ,
                        XYZ_emi_noise ,
                        T_rec_noise   ,
                        TAT           ,
                        TWTT_noise    ,
                        E_emi_clean   ,
                        XYZ_emi_clean ,
                        T_emi_clean   ,
                        E_rec_clean   ,
                        XYZ_rec_clean ,
                        T_rec_clean   ,
                        TWTT_clean    ]
            Fieldlis = ['E_emi_noise',
                        'X_emi_noise',
                        'Y_emi_noise',
                        'Z_emi_noise',
                        'T_emi_noise',
                        'E_rec_noise',
                        'X_rec_noise',
                        'Y_rec_noise',
                        'Z_rec_noise',
                        'T_rec_noise',
                        'TAT',
                        'TWTT_noise',
                        'E_emi_clean',
                        'X_emi_clean',
                        'Y_emi_clean',
                        'Z_emi_clean',
                        'T_emi_clean',
                        'E_rec_clean',
                        'X_rec_clean',
                        'Y_rec_clean',
                        'Z_rec_clean',
                        'T_rec_clean',
                        'TWTT_clean']     
            
            header = header + '\nfields : ' + ' '.join(Fieldlis)
            Mout   = np.column_stack(Datalis)
            if with_export:
                Mpath = os.path.join(path_exp,'_'.join((exp,ts))+'.'+idpxp+'.P.dat')
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
