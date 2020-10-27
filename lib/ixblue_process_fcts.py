#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  1 11:11:19 2020

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
import acouclass as acls

import numpy as np
import pandas as pd

from scipy.spatial.transform import Rotation


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
    
    diff_c = -(np.sqrt((xb - xe) ** 2 + (yb - ye) ** 2 + (zb - ze) ** 2) + np.sqrt((xb - xr) ** 2 + (yb - yr) ** 2 + (zb - zr) ** 2)) / c ** 2

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



def fct_obs_dir_vect(xbea_x_apri_in,
                     xbea_y_apri_in,
                     xbea_z_apri_in,
                     xrec_in):
        
    xbea_apri_in =  np.array([xbea_x_apri_in,
                              xbea_y_apri_in,
                              xbea_z_apri_in])
    
    dist = conv.dist(xbea_apri_in,xrec_in)
    
    vx = (xbea_x_apri_in - xrec_in[0]) / dist
    vy = (xbea_y_apri_in - xrec_in[1]) / dist
    vz = (xbea_z_apri_in - xrec_in[2]) / dist
    
    return vx , vy , vz
    
    
def fct_obs_dir_vect_deriv_ana(xbea_x_apri_in,
                               xbea_y_apri_in,
                               xbea_z_apri_in,
                               xrec_in):
    xbea = xbea_x_apri_in
    ybea = xbea_y_apri_in
    zbea = xbea_z_apri_in
    xrec,yrec,zrec = xrec_in
    
    vx_diff_x = ((xbea - xrec) ** 2 + (ybea - yrec) ** 2 + (zbea - zrec) ** 2) ** (-0.3e1 / 0.2e1) * (ybea ** 2 - 2 * yrec * ybea + yrec ** 2 + (zbea - zrec) ** 2)
    vx_diff_y = -(xbea - xrec) * (ybea - yrec) * ((xbea - xrec) ** 2 + (ybea - yrec) ** 2 + (zbea - zrec) ** 2) ** (-0.3e1 / 0.2e1)
    vx_diff_z = -(xbea - xrec) * (zbea - zrec) * ((xbea - xrec) ** 2 + (ybea - yrec) ** 2 + (zbea - zrec) ** 2) ** (-0.3e1 / 0.2e1)
    
    vy_diff_x = -(xbea - xrec) * (ybea - yrec) * ((xbea - xrec) ** 2 + (ybea - yrec) ** 2 + (zbea - zrec) ** 2) ** (-0.3e1 / 0.2e1)
    vy_diff_y = ((xbea - xrec) ** 2 + (ybea - yrec) ** 2 + (zbea - zrec) ** 2) ** (-0.3e1 / 0.2e1) * (xbea ** 2 - 2 * xbea * xrec + xrec ** 2 + (zbea - zrec) ** 2)
    vy_diff_z = -(ybea - yrec) * (zbea - zrec) * ((xbea - xrec) ** 2 + (ybea - yrec) ** 2 + (zbea - zrec) ** 2) ** (-0.3e1 / 0.2e1)
    
    
    vz_diff_x = -(xbea - xrec) * (zbea - zrec) * ((xbea - xrec) ** 2 + (ybea - yrec) ** 2 + (zbea - zrec) ** 2) ** (-0.3e1 / 0.2e1)
    vz_diff_y = -(zbea - zrec) * (ybea - yrec) * ((xbea - xrec) ** 2 + (ybea - yrec) ** 2 + (zbea - zrec) ** 2) ** (-0.3e1 / 0.2e1)
    vz_diff_z = ((xbea - xrec) ** 2 + (ybea - yrec) ** 2 + (zbea - zrec) ** 2) ** (-0.3e1 / 0.2e1) * (xbea ** 2 - 2 * xbea * xrec + xrec ** 2 + (ybea - yrec) ** 2)

    out_tup = (vx_diff_x,vx_diff_y,vx_diff_z,
               vy_diff_x,vy_diff_y,vy_diff_z,
               vz_diff_x,vz_diff_y,vz_diff_z)
    
    return out_tup
    
    
    
    
    
    
    
    

def fct_obs_dir_vect_deriv_symbo(xbea_x_apri_in,
                                 xbea_y_apri_in,
                                 xbea_z_apri_in,
                                 xrec_in):

    xbea,ybea,zbea = sympy.symbols('xbea ybea zbea')
    xrec,yrec,zrec = sympy.symbols('xrec yrec zrec')
    
    dist = sympy.sqrt((xbea-xrec)**2 + (ybea-yrec)**2 + (zbea-zrec)**2)
    
    vx = (xbea - xrec) / dist
    vy = (ybea - yrec) / dist
    vz = (zbea - zrec) / dist
    
    vx_diff_x = sympy.diff(vx,xbea)
    vx_diff_y = sympy.diff(vx,ybea)
    vx_diff_z = sympy.diff(vx,zbea)

    vy_diff_x = sympy.diff(vy,xbea)
    vy_diff_y = sympy.diff(vy,ybea)
    vy_diff_z = sympy.diff(vy,zbea)

    vz_diff_x = sympy.diff(vz,xbea)
    vz_diff_y = sympy.diff(vz,ybea)
    vz_diff_z = sympy.diff(vz,zbea)
    
    dict_for_eval = dict()    

    dict_for_eval['xbea']=xbea_x_apri_in
    dict_for_eval['ybea']=xbea_y_apri_in
    dict_for_eval['zbea']=xbea_z_apri_in
    dict_for_eval['xrec']=xrec_in[0]
    dict_for_eval['yrec']=xrec_in[1]
    dict_for_eval['zrec']=xrec_in[2]

    vx_diff_x_eva = vx_diff_x.evalf(subs=dict_for_eval)
    vx_diff_y_eva = vx_diff_y.evalf(subs=dict_for_eval)
    vx_diff_z_eva = vx_diff_z.evalf(subs=dict_for_eval)
    
    vy_diff_x_eva = vy_diff_x.evalf(subs=dict_for_eval)
    vy_diff_y_eva = vy_diff_y.evalf(subs=dict_for_eval)
    vy_diff_z_eva = vy_diff_z.evalf(subs=dict_for_eval)
    
    vz_diff_x_eva = vz_diff_x.evalf(subs=dict_for_eval)
    vz_diff_y_eva = vz_diff_y.evalf(subs=dict_for_eval)
    vz_diff_z_eva = vz_diff_z.evalf(subs=dict_for_eval)
    
    out_tup = (vx_diff_x_eva,vx_diff_y_eva,vx_diff_z_eva,
               vy_diff_x_eva,vy_diff_y_eva,vy_diff_z_eva,
               vz_diff_x_eva,vz_diff_y_eva,vz_diff_z_eva)
    
    return out_tup
   

def fct_obs_vlbi(B_in,S_in,c_in):
    dt = -1 * (np.dot(B_in,S_in)) / c_in
    return dt

def fct_obs_vlbi_deriv_symbo(B_in,S_in,c_in):

    bx,by,bz = sympy.symbols('bx by bz')    
    sx,sy,sz = sympy.symbols('sx sy sz')
    c = sympy.symbols('c')
    
    dt = -1 * (bx*sx + by*sy + bz*sz) / c
    
    dt_diff_sx = sympy.diff(dt,sx)
    dt_diff_sy = sympy.diff(dt,sy)
    dt_diff_sz = sympy.diff(dt,sz)
    
    dict_for_eval = dict()    

    dict_for_eval['bx']=B_in[0]
    dict_for_eval['by']=B_in[1]
    dict_for_eval['bz']=B_in[2]
    dict_for_eval['sx']=S_in[0]
    dict_for_eval['sy']=S_in[1]
    dict_for_eval['sz']=S_in[2]
    dict_for_eval['c']=c_in
    
    dt_diff_sx_eva = np.float64(dt_diff_sx.evalf(subs=dict_for_eval))
    dt_diff_sy_eva = np.float64(dt_diff_sy.evalf(subs=dict_for_eval))
    dt_diff_sz_eva = np.float64(dt_diff_sz.evalf(subs=dict_for_eval))
            
    return dt_diff_sx_eva,dt_diff_sy_eva,dt_diff_sz_eva


def fct_obs_vlbi_deriv_ana(B_in,S_in,c_in):
    
    dt_diff_sx_eva = - B_in[0] / c_in
    dt_diff_sy_eva = - B_in[1] / c_in
    dt_diff_sz_eva = - B_in[2] / c_in
    
    dt_diff_c_eva  = np.dot(B_in,S_in) / (c_in**2)
        
    return dt_diff_sx_eva,dt_diff_sy_eva,dt_diff_sz_eva,dt_diff_c_eva


def direction_vector_finder_cost_fct(Sin,
                                     DFepocin,
                                     lever_arm_RPY_dic,
                                     cin,
                                     id_tdc_ref = 1,
                                     twtt_key='TWTT',
                                     debug=False,
                                     short_output=True,
                                     with_estim_c=False):
    if len(DFepocin) != 4:
        return np.ones(3) * np.nan
    
    if np.sum(Sin) == 0:
        print("WARN: direction_vector_finder: Sin == 0, it is critical !!!")
        
    B = []
    O = []
    M = []
    Aline_stk = []        

    Bnorm = []
    Aline_nrm_stk = []        

    DFepoc_ref = DFepocin[DFepocin['ID_TDC'] == id_tdc_ref]
    t_ref  = DFepoc_ref[twtt_key].values[0]    
    
    for ilin , lin in DFepocin.iterrows():
        
        id_tdc = lin['ID_TDC']
        
        if id_tdc == id_tdc_ref:
            continue

        t_tdc  = lin[twtt_key]

        dt_obs = t_tdc - t_ref
        Baseline_in = (lever_arm_RPY_dic["TDC" + str(id_tdc)] - lever_arm_RPY_dic["TDC" + str(id_tdc_ref)] )
        ##### NO NORMALISATION OF Baseline_in OF COURSE !!!
        #Baseline_in = Baseline_in / np.linalg.norm(Baseline_in)
    
        dt_mod = fct_obs_vlbi(Baseline_in,Sin,cin) * 10**6 ### Bring it to microsec
        
        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        #####dirty mod
        #Seed = np.random.RandomState(42)
        #dt_mod = dt_obs + Seed.randn() * 10
        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        b = dt_obs - dt_mod

        if debug:
            print(id_tdc_ref,id_tdc,"dt_obs,dt_mod", dt_obs,dt_mod,b)
        
        B.append(b)
        O.append(dt_obs)
        
        #Aline = fct_obs_vlbi_deriv_symbo(Baseline_in,Sin,cin)
        if with_estim_c:
            Aline = fct_obs_vlbi_deriv_ana(Baseline_in,Sin,cin*10**-6)
        else:
            Aline = fct_obs_vlbi_deriv_ana(Baseline_in,Sin,cin*10**-6)[:-1]
            

        #Aline_stk.append(Aline)
        Aline_stk.append(Aline)
    
    ###### norm
    n     = np.linalg.norm(Sin)
    
    bnorm = n - 1
    Bnorm.append(bnorm)
    
    xdiff = Sin[0] / n
    ydiff = Sin[1] / n
    zdiff = Sin[2] / n

    if with_estim_c:
        Aline_nrm_stk.append(np.array([xdiff,ydiff,zdiff,0.]))
    else:
        Aline_nrm_stk.append(np.array([xdiff,ydiff,zdiff]))
            
    O  = np.array(O)     
    B1 = np.array(B)
    B2 = np.array(Bnorm)
    B  = np.hstack((B1,B2))
    M  = np.array(M)
    A1 = np.vstack(Aline_stk)
    A2 = np.vstack(Aline_nrm_stk)
    A  = np.vstack((A1,A2))
    
    if debug:
        print("sum(B**2)",np.sum(B1**2))
    
    if short_output:
        return np.sum(B1**2)
    else:
        return B,O,M,A


def direction_vector_finder_minimize(Sin,
                                     DFepocin,
                                     lever_arm_RPY_dic,
                                     cin,
                                     id_tdc_ref = 1,
                                     twtt_key='TWTT',
                                     debug=True):
    # meth in ('Powell',"CG",'Nelder-Mead')
    meth = 'Nelder-Mead'
    
    RES = scipy.optimize.minimize(direction_vector_finder_cost_fct,
                                  Sin,
                                  args=(DFepocin,lever_arm_RPY_dic,cin),
                                  method=meth)
    if debug:
        print(RES)
    
    return RES.x

def direction_vector_finder_lsq(Sin,
                            DFepocin,
                            lever_arm_RPY_dic,
                            cin,
                            id_tdc_ref = 1,
                            twtt_key='TWTT',
                            debug=False,
                            with_estim_c=True):
    if len(DFepocin) != 4:
        return np.ones(3) * np.nan
    
    if np.sum(Sin) == 0:
        print("WARN: direction_vector_finder: Sin == 0, it is critical !!!")

    
    ############## COST FCT #####################################
    B,O,M,A = direction_vector_finder_cost_fct(Sin,
                                               DFepocin,
                                               lever_arm_RPY_dic,
                                               cin,
                                               id_tdc_ref,
                                               twtt_key,
                                               debug,
                                               short_output=False,
                                               with_estim_c=with_estim_c)
        
    ##############################################################
    
    At = np.transpose(A)

    if debug:
        print(A)
    
    
    _,_,P = geok.weight_mat([1,10**-9],[3,1])

    N = At.dot(P).dot(A)
    Ninv = scipy.linalg.inv(N)
    
    dX = Ninv.dot(At).dot(P).dot(B)
    
    if debug:
        print("dX",dX)
    
    
    Snew = Sin + dX[:3]
    if with_estim_c:
        cnew = cin + dX[3]*10**6
    
    #Snew = Snew / np.linalg.norm(Snew)
    
    if debug:
        Rtheo = B - A.dot(dX)
        print("Resid Theo:",Rtheo)
        
        M = []
        for ilin , lin in DFepocin.iterrows():
            id_tdc = lin['ID_TDC']
            if id_tdc == id_tdc_ref:
                continue
            t_tdc  = lin[twtt_key]
            Baseline_in = (lever_arm_RPY_dic["TDC" + str(id_tdc)] - lever_arm_RPY_dic["TDC" + str(id_tdc_ref)] )
            dt_mod_res = fct_obs_vlbi(Baseline_in,Snew,cin) * 10**6 ### Bring it to microsec
            M.append(dt_mod_res)
            
        Rnatu = np.array(O) - np.array(M)
        print("Resid Natu:",Rnatu)
        
        if debug:
            print("sum(B**2)",np.sum(Rnatu**2))
    
    
    print("DV inp",Sin)
    print("DV out",Snew)
    
    return Snew


######################## site azimut finder 

def ned2site_azim(ned):
    """
    Grewall p 460
    """
    n,e,d = np.array(ned)
    r    = np.linalg.norm(ned)
    site = np.arcsin(d/r)  #aka phi
    azim = np.arctan2(e,n)   #aka theta
    return np.squeeze(np.column_stack((site,azim,r)))


def site_azim2ned(sar):
    """
    Grewall p 460
    """
    site,azim,r = np.array(sar)
    n = r*np.cos(azim)*np.cos(site)
    e = r*np.sin(azim)*np.cos(site)
    d = r*np.sin(site)
    return np.squeeze(np.column_stack((n,e,d)))


def ned2dir_cos(ned,degrees_out=False):
    ned_norm = np.linalg.norm(ned)
    Dir_cos = ned / ned_norm
    Angle = np.arccos(Dir_cos)
    if degrees_out:
        Angle = np.rad2deg(Angle)
    return Dir_cos , Angle
    

def fct_obs_vlbi_siteazim(B_in,SitAzi_in,c_in):
    bx,by,bz = B_in
    if len(SitAzi_in) == 3: 
        si,az,r = SitAzi_in
    else:
        si,az = SitAzi_in
        r = 1.
        
    dt = - (bx * r * np.cos(az) * np.cos(si) + by * r * np.sin(az) * np.cos(si) + bz * r * np.sin(si)) / c_in
    return dt

def fct_obs_vlbi_siteazim_deriv_ana(B_in,SitAzi_in,c_in):
    bx,by,bz = B_in
    if len(SitAzi_in) == 3: 
        si,az,r = SitAzi_in
    else:
        si,az = SitAzi_in
        r = 1.
        
    dt_diff_azim_eva = -(-bx * r * np.sin(az) * np.cos(si) + by * r * np.cos(az) * np.cos(si)) / c_in
    dt_diff_site_eva = -(-bx * r * np.cos(az) * np.sin(si) - by * r * np.sin(az) * np.sin(si) + bz * r * np.cos(si)) / c_in
    dt_diff_c_eva    =   (bx * r * np.cos(az) * np.cos(si) + by * r * np.sin(az) * np.cos(si) + bz * r * np.sin(si)) / c_in ** 2
    
    return dt_diff_azim_eva, dt_diff_site_eva, dt_diff_c_eva


def siteazim_finder_cost_fct(SitAziin,
                             DFepocin,
                             lever_arm_RPY_dic,
                             cin,
                             id_tdc_ref = 1,
                             twtt_key='TWTT',
                             debug=False,
                             short_output=True,
                             with_estim_c=True):
    
    if len(DFepocin) != 4:
        return np.ones(3) * np.nan
        
    B = []
    O = []
    M = []
    Aline_stk = []        

    Bnorm = []
    Aline_nrm_stk = []        

    DFepoc_ref = DFepocin[DFepocin['ID_TDC'] == id_tdc_ref]
    t_ref  = DFepoc_ref[twtt_key].values[0]    
    
    for ilin , lin in DFepocin.iterrows():
        
        id_tdc = lin['ID_TDC']
        
        if id_tdc == id_tdc_ref:
            continue

        t_tdc  = lin[twtt_key]

        dt_obs = t_tdc - t_ref
        Baseline_in = (lever_arm_RPY_dic["TDC" + str(id_tdc)] - lever_arm_RPY_dic["TDC" + str(id_tdc_ref)] )
        ##### NO NORMALISATION OF Baseline_in OF COURSE !!!
        #Baseline_in = Baseline_in / np.linalg.norm(Baseline_in)
    
        dt_mod = fct_obs_vlbi_siteazim(Baseline_in,SitAziin,cin) * 10**6 ### Bring it to microsec
        
        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        #####dirty mod
        #Seed = np.random.RandomState(42)
        #dt_mod = dt_obs + Seed.randn() * 10
        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        b = dt_obs - dt_mod

        if debug:
            print(id_tdc_ref,id_tdc,"dt_obs,dt_mod", dt_obs,dt_mod,b)
        
        B.append(b)
        O.append(dt_obs)
        
        #Aline = fct_obs_vlbi_deriv_symbo(Baseline_in,Sin,cin)
        if with_estim_c:
            Aline = fct_obs_vlbi_siteazim_deriv_ana(Baseline_in,SitAziin,cin*10**-6)
        else:
            Aline = fct_obs_vlbi_siteazim_deriv_ana(Baseline_in,SitAziin,cin*10**-6)[:-1]
            

        #Aline_stk.append(Aline)
        Aline_stk.append(Aline)
            
    O  = np.array(O)     
    B  = np.array(B)
    M  = np.array(M)
    A  = np.vstack(Aline_stk)
    
    if debug:
        print("sum(B**2)",np.sum(B**2))
    
    if short_output:
        return np.sum(B**2)
    else:
        return B,O,M,A


def siteazim_finder_lsq(SitAziin,
                        DFepocin,
                        lever_arm_RPY_dic,
                        cin,
                        id_tdc_ref = 1,
                        twtt_key='TWTT',
                        debug=True,
                        with_estim_c=False):
    
    if len(DFepocin) != 4:
        return np.ones(3) * np.nan
    
    ############## COST FCT #####################################
    B,O,M,A = siteazim_finder_cost_fct(SitAziin,
                                       DFepocin,
                                       lever_arm_RPY_dic,
                                       cin,
                                       id_tdc_ref,
                                       twtt_key,
                                       debug,
                                       short_output=False,
                                       with_estim_c=with_estim_c)    
    ##############################################################
    
    At = np.transpose(A)

    if debug:
        print(A)
    
    
    #_,_,P = geok.weight_mat([1,10**-9],[3,1])
    
    P = np.eye(len(A))

    N = At.dot(P).dot(A)
    Ninv = scipy.linalg.inv(N)
    
    dX = Ninv.dot(At).dot(P).dot(B)
    
    if debug:
        print("dX",dX)
    
    
    SitAzinew = SitAziin + dX[:2]
    if with_estim_c:
        cnew = cin + dX[2]*10**6
    else:
        cnew = cin
        
    if debug:
        Rtheo = B - A.dot(dX)
        print("Resid Theo:",Rtheo)
        
        Bnatu,_,_,_ = siteazim_finder_cost_fct(SitAzinew,
                                               DFepocin,
                                               lever_arm_RPY_dic,
                                               cnew,
                                               id_tdc_ref,
                                               twtt_key,
                                               debug,
                                               short_output=False,
                                               with_estim_c=with_estim_c)    
        print("Resid Natu:",Bnatu)
        
    print("SiAz rad inp",SitAziin)
    print("SiAz rad out",SitAzinew)
    print("SiAz deg inp",np.rad2deg(SitAziin))
    print("SiAz deg out",np.rad2deg(SitAzinew))
    return SitAzinew


    
def c_empiric(DFepoc_in,lever_arm_RPY_dic1,twtt_key="TWTT"):
    print("C EMPIRIC IS A WORNG VALIDATION !!!!!")
    
    DFtwtt1 = DFepoc_in[["ID_TDC",twtt_key]].copy()
    DFtwtt1.set_index('ID_TDC',inplace=True)
    DFtwtt_min = DFtwtt1[DFtwtt1[twtt_key] == DFtwtt1[twtt_key].min()]
    idx_twtt_min = DFtwtt_min.index[0]
    TupStk = []
    for i,lin in DFtwtt1.iterrows():
        if i == idx_twtt_min:
            continue
        d = np.linalg.norm(lever_arm_RPY_dic1["TDC" + str(i)] - 
                           lever_arm_RPY_dic1["TDC" + str(idx_twtt_min)])
        t = (lin.values[0] - DFtwtt_min.values[0]) * 10**-6
        c = d/t
        tup = (idx_twtt_min,i,c[0])
        TupStk.append(tup)
    
    
    
    return pd.DataFrame(TupStk)


def generate_DFtwtt(DFepocin):
    """
    Simpled version of a DFepoc
    """
    DFtwtt = DFepocin[["TWTT","Ph_raw","ID_TDC"]].copy()
    DFtwtt.sort_values("TWTT",inplace=True)    
    DFtwtt.set_index("ID_TDC",inplace=True)
    return DFtwtt
##############################################################################
