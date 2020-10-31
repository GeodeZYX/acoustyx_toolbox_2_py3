#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 12:20:07 2020

@author: psakicki
"""

#### Import geodeZYX
from geodezyx import *                   # Import the GeodeZYX modules
from geodezyx.externlib import *         # Import the external modules
from geodezyx.megalib.megalib import *   # Import the legacy modules names

#### Import Seafloor positionning module
from pygoat.seafloorpos import ixblue_fcts

import numpy as np
import scipy

import matplotlib.pyplot as plt




export_dir_orig = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/01_/"
experience_name = "PAMELI_BREST_vJupyter_5_4transducers_new_lever_arms"
### create the dir if doesn't exists
export_dir = os.path.join(export_dir_orig,experience_name)
utils.create_dir(export_dir)



#### GAPS files
p_gaps_10 = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/01_DATA/07-2019_IXBLUE_BREST/0_RAW_DATA/data-gaps/data - jeudi 25 juillet 2019/10 - Repeater Gaps - 14h10 - au barycentre IIS.txt"
p_gaps_11 = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/01_DATA/07-2019_IXBLUE_BREST/0_RAW_DATA/data-gaps/data - jeudi 25 juillet 2019/11 - Repeater Gaps - 14h46 - au barycentre - CIS.txt"
p_gaps_01 = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/01_DATA/07-2019_IXBLUE_BREST/0_RAW_DATA/data-gaps/data - jeudi 25 juillet 2019/04 - Repeater Gaps - 11h58 - box in TP1.txt"
p_gaps_02 = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/01_DATA/07-2019_IXBLUE_BREST/0_RAW_DATA/data-gaps/data - jeudi 25 juillet 2019/02 - Repeater Gaps - 11h09 - box in TP2.txt"
p_gaps_03 = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/01_DATA/07-2019_IXBLUE_BREST/0_RAW_DATA/data-gaps/data - jeudi 25 juillet 2019/06 - Repeater Gaps - 12h45 - box in TP3.txt"

#### Position files
p_posi_02  = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/01_DATA/07-2019_IXBLUE_BREST/4_RESULTS/2019-07-25/Manip_5_1/2019_07_25_Manip5_1_gins_ppp.csv"
p_posi_02b = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/01_DATA/07-2019_IXBLUE_BREST/4_RESULTS/2019-07-25/Manip_5_1/2019_07_25_Manip5_1_nrcan.csv"
p_posi_02c = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/01_DATA/07-2019_IXBLUE_BREST/4_RESULTS/2019-07-25/2019_07_25_Manip5_123_nrcan.csv"
p_posi_03  = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/01_DATA/07-2019_IXBLUE_BREST/4_RESULTS/2019-07-25/Manip_4_1/2019_07_25_Manip4_1_gins_ppp.csv"
p_posi_03b = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/01_DATA/07-2019_IXBLUE_BREST/4_RESULTS/2019-07-25/Manip_4_1/2019_07_25_Manip4_1_nrcan.csv"
p_posi_03c = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/01_DATA/07-2019_IXBLUE_BREST/4_RESULTS/2019-07-25/2019_07_25_Manip4_123_nrcan.csv"
p_posi_01  = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/01_DATA/07-2019_IXBLUE_BREST/4_RESULTS/2019-07-25/Manip_3_1/2019_07_25_Manip3_1_gins_ppp.csv"
p_posi_01b = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/01_DATA/07-2019_IXBLUE_BREST/4_RESULTS/2019-07-25/Manip_3_1/2019_07_25_Manip3_1_nrcan.csv"
p_posi_01c = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/01_DATA/07-2019_IXBLUE_BREST/4_RESULTS/2019-07-25/2019_07_25_Manip3_123_nrcan.csv"

p_posi_11  = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/01_DATA/07-2019_IXBLUE_BREST/4_RESULTS/2019-07-25/Manip_2_2/2019_07_25_Manip2_2_gins_ppp.csv"
p_posi_11b = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/01_DATA/07-2019_IXBLUE_BREST/4_RESULTS/2019-07-25/Manip_2_2/2019_07_25_Manip2_2_rtklib_nrcancoord.csv"

p_gaps_tat = p_gaps_11

#### Barycentre
p_gaps = p_gaps_11
p_posi = p_posi_11b

#### Boxin 
p_gaps = p_gaps_01
p_posi = p_posi_01c

p_gaps = p_gaps_02
p_posi = p_posi_02c

p_gaps = p_gaps_03
p_posi = p_posi_03c

DFposi_post = pd.read_csv(p_posi,header=0,sep=",")

if 'Unnamed: 0' in DFposi_post.columns:
    DFposi_post = DFposi_post.rename(columns={'Unnamed: 0':'date'})

########### POSI EXTE
dsec = 18

DFposi_post["date"] = pd.to_datetime(DFposi_post['date']) -  dt.timedelta(seconds=dsec)
XYZpost             = conv.GEO2XYZ_vector(DFposi_post[["lat","lon","values"]])
DFposi_post["X"],DFposi_post["Y"],DFposi_post["Z"] = XYZpost[:,0],XYZpost[:,1],XYZpost[:,2]

if 1:
    print("WARN POST PROCESSING IN RGF93")
    HPparam = reffram.itrf_helmert_get_parameters("ITRF2014","ETRF2000",
                                              verbose=False,
                                              convert=True)
    XYZ_post_RGF93 = reffram.itrf_speed_calc(DFposi_post[["X"]],
                                             DFposi_post[["Y"]],
                                             DFposi_post[["Z"]],2019.6,
                                             -1.16615726390990e-02,
                                             1.76696025246000e-02,
                                             1.08543846164409e-02,
                                             2009.0)
    XYZ_post_RGF93 = np.column_stack(XYZ_post_RGF93)
    XYZ_post_RGF93 = reffram.itrf_helmert_trans(XYZ_post_RGF93,2009.,*HPparam)
    DFposi_post[["X","Y","Z"]] = XYZ_post_RGF93
    

xyz_ref = np.mean(XYZpost,axis=0)


########### POSI EXTE
string_PTSAG         = "^\$PTSAG"
DFgaps_posi_ixblue   = ixblue_fcts.read_ixblue_data(p_gaps, string_PTSAG,True)
DFposi_ixblue        = pd.DataFrame()

TTT = DFgaps_posi_ixblue.apply(ixblue_fcts.ixblue_time_conv,axis=1,args=(2,))
DFposi_ixblue["ID"]    = DFgaps_posi_ixblue[6]
DFposi_ixblue["date"]  = TTT
DFposi_ixblue["lat"]   = DFgaps_posi_ixblue.apply(ixblue_fcts.ixblue_coord_conv,axis=1,args=(7,))
DFposi_ixblue["lon"]   = DFgaps_posi_ixblue.apply(ixblue_fcts.ixblue_coord_conv,axis=1,args=(9,))
DFposi_ixblue["depth"] = DFgaps_posi_ixblue[12]


DFposi_ixblue = DFposi_ixblue[DFposi_ixblue.ID == 0].copy()


XYZixblue = conv.GEO2XYZ_vector(DFposi_ixblue[["lat","lon","depth"]])
DFposi_ixblue["X"],DFposi_ixblue["Y"],DFposi_ixblue["Z"] = XYZixblue[:,0],XYZixblue[:,1],XYZixblue[:,2]


ENUtmp = conv.XYZ2ENU_vector(DFposi_ixblue[["X","Y","Z"]], xyz_ref)
DFposi_ixblue["E"],DFposi_ixblue["N"],DFposi_ixblue["U"] = ENUtmp[:,0],ENUtmp[:,1],ENUtmp[:,2]

ENUtmp = conv.XYZ2ENU_vector(DFposi_post[["X","Y","Z"]], xyz_ref)
DFposi_post["E"],DFposi_post["N"],DFposi_post["U"] = ENUtmp[:,0],ENUtmp[:,1],ENUtmp[:,2]
  


Dpst = DFposi_post
Dixb = DFposi_ixblue

Tref = Dpst["date"]

I_XYZ_post   = ixblue_fcts.interpolator_position(Dpst[["X","Y","Z"]],Dpst["date"])
I_XYZ_ixblue = ixblue_fcts.interpolator_position(Dixb[["X","Y","Z"]],Dixb["date"])

I_ENU_post   = ixblue_fcts.interpolator_position(Dpst[["E","N","U"]],Dpst["date"])
I_ENU_ixblue = ixblue_fcts.interpolator_position(Dixb[["E","N","U"]],Dixb["date"])

##############################################################################
fig10,Ax10 = plt.subplots(3,1)

Ax10[0].plot(Dpst.date,Dpst.X)
Ax10[1].plot(Dpst.date,Dpst.Y)
Ax10[2].plot(Dpst.date,Dpst.Z)
Ax10[0].plot(Dixb.date,Dixb.X)
Ax10[1].plot(Dixb.date,Dixb.Y)
Ax10[2].plot(Dixb.date,Dixb.Z)


##############################################################################
fig11,Ax11 = plt.subplots(3,1)

Ax11[0].plot(Dpst.date,Dpst.E)
Ax11[1].plot(Dpst.date,Dpst.N)
Ax11[2].plot(Dpst.date,Dpst.U)
Ax11[0].plot(Dixb.date,Dixb.E)
Ax11[1].plot(Dixb.date,Dixb.N)
Ax11[2].plot(Dixb.date,Dixb.U)

##############################################################################

fig30,ax30 = plt.subplots(1,1)

XYZ_Tref_pst = I_XYZ_post(Tref).T
XYZ_Tref_ixb = I_XYZ_ixblue(Tref).T

ENU_Tref_pst = I_ENU_post(Tref).T
ENU_Tref_ixb = I_ENU_ixblue(Tref).T

dDxyz = XYZ_Tref_pst - XYZ_Tref_ixb
Dxyz  = np.linalg.norm(dDxyz,axis=1)

dDenu = ENU_Tref_pst - ENU_Tref_ixb
Denu  = np.linalg.norm(dDenu,axis=1)

ax30.plot(Tref,Dxyz,"x")
ax30.plot(Tref,Denu,"+")

dDen = ENU_Tref_pst[:,:2] - ENU_Tref_ixb[:,:2]
Den  = np.linalg.norm(dDen,axis=1)

fig31,ax31 = plt.subplots(1,1)
ax31.plot(Tref,Den,"+")
ax31.set_title("Den(Tref)")

##############################################################################

Res , A , l   = reffram.helmert_trans_estim(ENU_Tref_pst[:1000],
                                            ENU_Tref_ixb[:1000])
ENU_Tref_pst2 = reffram.helmert_trans_apply(ENU_Tref_pst,Res)
ENU_Tref_pst2 - ENU_Tref_pst


##############################################################################

fig40,Ax40 = plt.subplots(3,1)

Ax40[0].plot(Tref , dDenu[:,0]) #- np.nanmean(dDenu[:,0]))
Ax40[1].plot(Tref , dDenu[:,1]) #- np.nanmean(dDenu[:,1]))
Ax40[2].plot(Tref , dDenu[:,2]) #- np.nanmean(dDenu[:,2]))










