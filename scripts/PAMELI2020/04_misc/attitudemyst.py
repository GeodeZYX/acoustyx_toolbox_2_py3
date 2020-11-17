#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 12:13:38 2020

@author: psakicki
"""

#### Import star style
from geodezyx import *                   # Import the GeodeZYX modules
from geodezyx.externlib import *         # Import the external modules
from geodezyx.megalib.megalib import *   # Import the legacy modules names

############ Attitude management ##############
  
    
    

Head  = DFgaps[11].values
Roll  = DFgaps[12].values
Pitch = - DFgaps[13].values



### We switch to Head Pitch Roll (zyx) ########
Euler_rec = np.column_stack((Head,Pitch,Roll))
Euler_rec_for_slerp = list(np.vstack((Euler_rec[0],Euler_rec)))

# We generate the Rotation objects at the emission epoch using a Slerp interpo
# the _for_slerp objects are a manual (and dumb) extrapolation for the 1st value
Trec_posix_for_slerp = np.insert(Trec_posix, 0, Trec_posix[0] - 1000.)
Trec_posix_for_slerp,Euler_rec_for_slerp = utils.sort_binom_list(Trec_posix_for_slerp,Euler_rec_for_slerp)

euler_order = "zyx"

Rot_rec_for_slerp = Rotation.from_euler(euler_order,Euler_rec_for_slerp,degrees=False)
Islerp = Slerp(Trec_posix_for_slerp,Rot_rec_for_slerp)

Rot_emi = Islerp(Temi_posix[:])
Rot_rec = Rotation.from_euler(euler_order,Euler_rec,degrees=False)

Rot_PTSAY = Islerp(DF_PTSAY["Tposix"])


##### Here are defined the Lever Arms ##############
## AHD: Acoustic Head
## REF: reference point of the Attitude Unit
# print("WRONG LEVER ARM")
lever_arm_RPY_REF2AHD = np.array([0. ,0., 0.])

lever_arm_RPY_REF2GPS = np.array([0.2,0.,0.])
lever_arm_RPY_REF2GPS = np.array([0.02,0.,-1.930])
##################################################
iNED_AHD_rec   = Rot_rec.apply(lever_arm_RPY_REF2GPS)
iNED_AHD_emi   = Rot_emi.apply(lever_arm_RPY_REF2GPS)
iNED_AHD_PTSAY = Rot_PTSAY.apply(lever_arm_RPY_REF2GPS)

DFout["Head"]  = Head
DFout["Pitch"] = Roll  
DFout["Roll"]  = Pitch 

DFout["iN_AHD_rec"] = iNED_AHD_rec[:,0]
DFout["iN_AHD_emi"] = iNED_AHD_emi[:,0]
DFout["iE_AHD_rec"] = iNED_AHD_rec[:,1]
DFout["iE_AHD_emi"] = iNED_AHD_emi[:,1]
DFout["iD_AHD_rec"] = iNED_AHD_rec[:,2]
DFout["iD_AHD_emi"] = iNED_AHD_emi[:,2]

DFout["N_AHD_emi"] = DFout["N_GPS_emi"] + iNED_AHD_emi[:,0]
DFout["E_AHD_emi"] = DFout["E_GPS_emi"] + iNED_AHD_emi[:,1]
DFout["D_AHD_emi"] = DFout["D_GPS_emi"] + iNED_AHD_emi[:,2]
DFout["N_AHD_rec"] = DFout["N_GPS_rec"] + iNED_AHD_rec[:,0]
DFout["E_AHD_rec"] = DFout["E_GPS_rec"] + iNED_AHD_rec[:,1]
DFout["D_AHD_rec"] = DFout["D_GPS_rec"] + iNED_AHD_rec[:,2]

############## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DF_PTSAY["iN_AHD"] = iNED_AHD_PTSAY[:,0]
DF_PTSAY["iE_AHD"] = iNED_AHD_PTSAY[:,1]
DF_PTSAY["iD_AHD"] = iNED_AHD_PTSAY[:,2]

if False:
    DF_PTSAY["N_AHD"] = DF_PTSAY["iN_AHD"] +  DF_PTSAY["N_GPS"]
    DF_PTSAY["E_AHD"] = DF_PTSAY["iE_AHD"] +  DF_PTSAY["E_GPS"]
    DF_PTSAY["D_AHD"] = DF_PTSAY["iD_AHD"] +  DF_PTSAY["D_GPS"]
else:
    DF_PTSAY["N_AHD"] = DF_PTSAY["N_GPS"]
    DF_PTSAY["E_AHD"] = DF_PTSAY["E_GPS"]
    DF_PTSAY["D_AHD"] = DF_PTSAY["D_GPS"]


DF_PTSAY["N_BEA"] = DF_PTSAY["iN_AHD2BEA"] +  DF_PTSAY["N_AHD"]
DF_PTSAY["E_BEA"] = DF_PTSAY["iE_AHD2BEA"] +  DF_PTSAY["E_AHD"]
DF_PTSAY["D_BEA"] = DF_PTSAY["iD_AHD2BEA"] +  DF_PTSAY["D_AHD"]


DFout = DFout.dropna()
