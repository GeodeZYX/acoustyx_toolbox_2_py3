#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 17:49:48 2020

@author: psakicki
"""

# In[2]:


#### Import geodeZYX
from geodezyx import *                   # Import the GeodeZYX modules
from geodezyx.externlib import *         # Import the external modules
from geodezyx.megalib.megalib import *   # Import the legacy modules names

#### Import Seafloor positionning module
from pygoat.seafloorpos import ixblue_fcts

import numpy as np
import scipy

import matplotlib.pyplot as plt

# ## Define the experiment names and output directories

# In[3]:


export_dir_orig = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/01_/"
experience_name = "PAMELI_BREST_vJupyter_5_4transducers_new_lever_arms"
### create the dir if doesn't exists
export_dir = os.path.join(export_dir_orig,experience_name)
utils.create_dir(export_dir)


# ## Get the input file paths

# In[4]:


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

# ## Read the position

# ### Import

# In[5]:


DFposi_post = pd.read_csv(p_posi,header=0,sep=",")

if 'Unnamed: 0' in DFposi_post.columns:
    DFposi_post = DFposi_post.rename(columns={'Unnamed: 0':'date'})
    
# ### Conversion to XYZ

# In[6]:
    
    
########### POSI EXTE
DFposi_post["date"] = pd.to_datetime(DFposi_post['date']) - dt.timedelta(seconds=21)
XYZpost = conv.GEO2XYZ_vector(DFposi_post[["lat","lon","values"]])
DFposi_post["X"],DFposi_post["Y"],DFposi_post["Z"] = XYZpost[:,0],XYZpost[:,1],XYZpost[:,2]

xyz_ref = np.mean(XYZpost,axis=0)

########### POSI iXBLUE
string_PTSAG         = "^\$PTSAG"
DFgaps_posi_ixblue   = ixblue_fcts.read_ixblue_data(p_gaps, string_PTSAG,True)
DFgaps_posi_ixblue_3 = pd.DataFrame()

TTT = DFgaps_posi_ixblue.apply(ixblue_fcts.ixblue_time_conv,axis=1,args=(2,))
DFgaps_posi_ixblue_3["ID"] = DFgaps_posi_ixblue[6]
DFgaps_posi_ixblue_3["date"]  = TTT
DFgaps_posi_ixblue_3["lat"]   = DFgaps_posi_ixblue.apply(ixblue_fcts.ixblue_coord_conv,axis=1,args=(7,))
DFgaps_posi_ixblue_3["lon"]   = DFgaps_posi_ixblue.apply(ixblue_fcts.ixblue_coord_conv,axis=1,args=(9,))
DFgaps_posi_ixblue_3["depth"] = DFgaps_posi_ixblue[12]

DFgaps_posi_ixblue_2 = DFgaps_posi_ixblue_3[DFgaps_posi_ixblue_3.ID == 0].copy()

XYZixblue = conv.GEO2XYZ_vector(DFgaps_posi_ixblue_2[["lat","lon","depth"]])
DFgaps_posi_ixblue_2["X"],DFgaps_posi_ixblue_2["Y"],DFgaps_posi_ixblue_2["Z"] = XYZixblue[:,0],XYZixblue[:,1],XYZixblue[:,2]

#############
ENU = conv.XYZ2ENU_vector(DFgaps_posi_ixblue_2[["X","Y","Z"]], xyz_ref)
DFgaps_posi_ixblue_2["E"],DFgaps_posi_ixblue_2["N"],DFgaps_posi_ixblue_2["U"] = ENU[:,0],ENU[:,1],ENU[:,2]

ENU = conv.XYZ2ENU_vector(DFposi_post[["X","Y","Z"]], xyz_ref)
DFposi_post["E"],DFposi_post["N"],DFposi_post["U"] = ENU[:,0],ENU[:,1],ENU[:,2]


# plt.figure()
# plt.plot(DFposi_post["date"],DFposi_post["X"])
# plt.plot(DFgaps_posi_ixblue_2["date"],DFgaps_posi_ixblue_2["X"])

# plt.figure()
# plt.plot(DFposi_post["date"],DFposi_post["Y"])
# plt.plot(DFgaps_posi_ixblue_2["date"],DFgaps_posi_ixblue_2["Y"])

# plt.figure()
# plt.plot(DFposi_post["date"],DFposi_post["Z"])
# plt.plot(DFgaps_posi_ixblue_2["date"],DFgaps_posi_ixblue_2["Z"])

# plt.figure()
# plt.plot(DFposi_post["date"],DFposi_post["Z"] - DFposi_post["Z"].mean())
# plt.plot(DFgaps_posi_ixblue_2["date"],DFgaps_posi_ixblue_2["Z"] - DFgaps_posi_ixblue_2["Z"].mean())

# plt.figure()
# plt.plot(DFposi_post["date"],DFposi_post["E"])
# plt.plot(DFgaps_posi_ixblue_2["date"],DFgaps_posi_ixblue_2["E"])
# plt.figure()
# plt.plot(DFposi_post["date"],DFposi_post["E"])
# plt.plot(DFgaps_posi_ixblue_2["date"],DFgaps_posi_ixblue_2["E"])

# plt.figure()
# plt.plot(DFposi_post["date"],DFposi_post["N"])
# plt.plot(DFgaps_posi_ixblue_2["date"],DFgaps_posi_ixblue_2["N"])
# plt.figure()
# plt.plot(DFposi_post["date"],DFposi_post["N"])
# plt.plot(DFgaps_posi_ixblue_2["date"],DFgaps_posi_ixblue_2["N"])

plt.figure()
plt.plot(DFposi_post["date"],DFposi_post["E"])
plt.plot(DFgaps_posi_ixblue_2["date"],DFgaps_posi_ixblue_2["E"])
plt.figure()
plt.plot(DFposi_post["date"],DFposi_post["N"])
plt.plot(DFgaps_posi_ixblue_2["date"],DFgaps_posi_ixblue_2["N"])

############# CHOIX POSI IXBLUE OR EXTERN

DFposi = DFgaps_posi_ixblue_2
bool_lever_arm_correction = True

DFposi = DFposi_post
bool_lever_arm_correction = False

DFposi = DFgaps_posi_ixblue_2
bool_lever_arm_correction = False

DFposi = DFposi_post
bool_lever_arm_correction = True

estimate_ties = False 

# ### Here we define the arbitrary center of the network

# In[7]:




# ### Conversion to ENU


# In[9]:

I_XYZ = ixblue_fcts.interpolator_position(DFposi[["X","Y","Z"]],
                                          DFposi["date"])

I_ENU = ixblue_fcts.interpolator_position(DFposi[["E","N","U"]],
                                          DFposi["date"])

# ### A priori coorinates of the Beacons

# In[11]:


string_PIXOG_twtt = "^\$PIXOG,PPC,DETEC"
DFgaps_raw = ixblue_fcts.read_ixblue_data(p_gaps, string_PIXOG_twtt,True)

string_PTSAY = "^\$PTSAY"
DFgaps_ptsa  = ixblue_fcts.read_ixblue_data(p_gaps, string_PTSAY,True)


DFgaps_ptsa_2 = DFgaps_ptsa.copy() 
DFgaps_ptsa_3 = pd.DataFrame()

TT = DFgaps_ptsa_2.apply(ixblue_fcts.ixblue_time_conv,axis=1,args=(2,))
DFgaps_ptsa_3["ID_BEA"] = DFgaps_ptsa_2[6]
DFgaps_ptsa_3["date_rec"] = TT
DFgaps_ptsa_3["date_emi"] = TT
DFgaps_ptsa_3["dN"]     = DFgaps_ptsa_2[7]
DFgaps_ptsa_3["dE"]     = DFgaps_ptsa_2[8]
DFgaps_ptsa_3["dD"]     = DFgaps_ptsa_2[10]

TATdict = ixblue_fcts.read_ixblue_TAT(p_gaps_tat)


# ## Do preliminary processing on the DataFrame

# ### Generate a userfriendly DataFrame

# #### Average transducer = if False, keep independently each head (4 transducers) if True, average the 4 tranducers

# In[12]:

DFixblue1 = ixblue_fcts.iXblueDF_prepro_step1_basic(DFgaps_raw,TATdict,
                                                    average_transducers=False)

# ### Interpolate the position at emission and reception epochs

# In[13]:


### On XYZ component
DFixblue2 = ixblue_fcts.iXblueDF_prepro_step2_position(DFixblue1,I_XYZ,"XYZ","GPS")
### On ENU component
DFixblue2 = ixblue_fcts.iXblueDF_prepro_step2_position(DFixblue2,I_ENU,"ENU","GPS")

### On XYZ component
DFgaps_ptsa_4 = ixblue_fcts.iXblueDF_prepro_step2_position(DFgaps_ptsa_3,I_XYZ,"XYZ","GPS")
### On ENU component
DFgaps_ptsa_4 = ixblue_fcts.iXblueDF_prepro_step2_position(DFgaps_ptsa_4,I_ENU,"ENU","GPS")



# ### Generate the Down component

# In[14]:


DFgaps_ptsa_4 = ixblue_fcts.ixblueDF_prepro_step3_down_component(DFgaps_ptsa_4)


# ## Interpolate the attitude and apply the lever arms 

# ### Generate attitude interpolators

# In[15]:
I_Att = ixblue_fcts.interpolator_attitude(DFixblue2['head_rec'],
                                          DFixblue2['pitc_rec'],
                                          DFixblue2['roll_rec'],
                                          DFixblue2['date_rec'])


############# CHOIX POSI IXBLUE OR EXTERN

# DFposi = DFgaps_posi_ixblue_2
# bool_lever_arm_correction = True

# DFposi = DFposi_post
# bool_lever_arm_correction = False

# DFposi = DFgaps_posi_ixblue_2
# bool_lever_arm_correction = False

# DFposi = DFposi_post
# bool_lever_arm_correction = True


estimate_ties = False 


# ### Here we define the arbitrary center of the network

# In[7]:




# ### Conversion to ENU


# In[9]:

# I_XYZ = ixblue_fcts.interpolator_position(DFposi[["X","Y","Z"]],
#                                           DFposi["date"])

# I_ENU = ixblue_fcts.interpolator_position(DFposi[["E","N","U"]],
#                                           DFposi["date"])

I_XYZ_post   = ixblue_fcts.interpolator_position(DFposi_post[["X","Y","Z"]],
                                          DFposi_post["date"])

I_XYZ_ixblue = ixblue_fcts.interpolator_position(DFgaps_posi_ixblue_2[["X","Y","Z"]],
                                          DFgaps_posi_ixblue_2["date"])

I_ENU_post   = ixblue_fcts.interpolator_position(DFposi_post[["E","N","U"]],
                                          DFposi_post["date"])

I_ENU_ixblue = ixblue_fcts.interpolator_position(DFgaps_posi_ixblue_2[["E","N","U"]],
                                          DFgaps_posi_ixblue_2["date"])

plt.figure()
plt.plot(DFposi_post["date"],DFposi_post["X"])
plt.plot(DFgaps_posi_ixblue_2["date"],DFgaps_posi_ixblue_2["X"])

plt.figure()
plt.plot(DFposi_post["date"],DFposi_post["Y"])
plt.plot(DFgaps_posi_ixblue_2["date"],DFgaps_posi_ixblue_2["Y"])

plt.figure()
plt.plot(DFposi_post["date"],DFposi_post["Z"])
plt.plot(DFgaps_posi_ixblue_2["date"],DFgaps_posi_ixblue_2["Z"])






# ### A priori coorinates of the Beacons

# In[11]:
string_PIXOG_twtt = "^\$PIXOG,PPC,DETEC"
DFgaps_raw = ixblue_fcts.read_ixblue_data(p_gaps, string_PIXOG_twtt,True)

string_PTSAY = "^\$PTSAY"
DFgaps_ptsa  = ixblue_fcts.read_ixblue_data(p_gaps, string_PTSAY,True)


DFgaps_ptsa_2 = DFgaps_ptsa.copy() 
DFgaps_ptsa_3 = pd.DataFrame()

TT = DFgaps_ptsa_2.apply(ixblue_fcts.ixblue_time_conv,axis=1,args=(2,))
DFgaps_ptsa_3["ID_BEA"] = DFgaps_ptsa_2[6]
DFgaps_ptsa_3["date_rec"] = TT
DFgaps_ptsa_3["date_emi"] = TT
DFgaps_ptsa_3["dN"]     = DFgaps_ptsa_2[7]
DFgaps_ptsa_3["dE"]     = DFgaps_ptsa_2[8]
DFgaps_ptsa_3["dD"]     = DFgaps_ptsa_2[10]

TATdict = ixblue_fcts.read_ixblue_TAT(p_gaps_tat)


# ## Do preliminary processing on the DataFrame

# ### Generate a userfriendly DataFrame

# #### Average transducer = if False, keep independently each head (4 transducers) if True, average the 4 tranducers

# In[12]:

DFixblue1 = ixblue_fcts.iXblueDF_prepro_step1_basic(DFgaps_raw,TATdict,
                                                    average_transducers=False)

# ### Interpolate the position at emission and reception epochs

# In[13]:

### On XYZ component
DFixblue2 = ixblue_fcts.iXblueDF_prepro_step2_position(DFixblue1,I_XYZ,"XYZ","GPS")
### On ENU component
DFixblue2 = ixblue_fcts.iXblueDF_prepro_step2_position(DFixblue2,I_ENU,"ENU","GPS")

### On XYZ component
DFgaps_ptsa_4 = ixblue_fcts.iXblueDF_prepro_step2_position(DFgaps_ptsa_3,I_XYZ,"XYZ","GPS")
### On ENU component
DFgaps_ptsa_4 = ixblue_fcts.iXblueDF_prepro_step2_position(DFgaps_ptsa_4,I_ENU,"ENU","GPS")

# ### Generate the Down component

# In[14]:


DFgaps_ptsa_4 = ixblue_fcts.ixblueDF_prepro_step3_down_component(DFgaps_ptsa_4)


# ## Interpolate the attitude and apply the lever arms 

# ### Generate attitude interpolators

# In[15]:
I_Att = ixblue_fcts.interpolator_attitude(DFixblue2['head_rec'],
                                          DFixblue2['pitc_rec'] * 0,
                                          DFixblue2['roll_rec'] * 0,
                                          DFixblue2['date_rec'])

#### We check here the attitude extrapo veracity
plt.figure()

plt.plot(DFixblue2['date_rec'],DFixblue2['head_rec'],"-x")
Rot = I_Att(DFixblue2['date_rec'] + dt.timedelta(seconds=0.5))
plt.plot(DFixblue2['date_rec'] + dt.timedelta(seconds=0.5),Rot.as_euler("zyx")[:,0],"+-")

TTTref = DFixblue2['date_rec'] 

ENU_Work_post   = I_ENU_post(TTTref)
ENU_Work_ixblue = I_ENU_ixblue(TTTref)

XYZ_Work_post   = I_XYZ_post(TTTref)
XYZ_Work_ixblue = I_XYZ_ixblue(TTTref)

iENU_Work_post_in_ixblue_Ref = ENU_Work_post - ENU_Work_ixblue
dENU = iENU_Work_post_in_ixblue_Ref.T

dXYZ =  XYZ_Work_post - XYZ_Work_ixblue
dXYZ = dXYZ.T



D1ENU = np.linalg.norm(dENU[:,:2],axis=1)
D1XYZ = np.linalg.norm(dXYZ[:,:2],axis=1)


iENU_Work_post_in_ixblue_Ref = iENU_Work_post_in_ixblue_Ref.T

iRPY = Rot.apply(iENU_Work_post_in_ixblue_Ref)

plt.figure()
plt.plot(TTTref,iENU_Work_post_in_ixblue_Ref[:,0])
plt.plot(TTTref,iENU_Work_post_in_ixblue_Ref[:,1])
plt.figure()

Rotplt = I_Att(TTTref)
Head1 = Rotplt.as_euler("zyx")[:,0] 
plt.plot(TTTref,Head1,"+-",label="Head")
plt.plot(TTTref,iRPY[:,0],label="N")
plt.plot(TTTref,iRPY[:,1],label="E")

Head2 = np.arctan2(iENU_Work_post_in_ixblue_Ref[:,0],
                   iENU_Work_post_in_ixblue_Ref[:,1])
plt.plot(TTTref,Head2,label="Head2")

plt.legend()

plt.figure()
plt.plot(iENU_Work_post_in_ixblue_Ref[:,0],iENU_Work_post_in_ixblue_Ref[:,1],"x")

plt.figure()
D2 = np.sqrt(iENU_Work_post_in_ixblue_Ref[:,0]**2 + iENU_Work_post_in_ixblue_Ref[:,1]**2)
plt.plot(D2)
plt.suptitle("Dist")


plt.figure()
plt.plot(TTTref,np.sin(Head2) * iENU_Work_post_in_ixblue_Ref[:,1],"x",label="E")
plt.plot(TTTref,np.cos(Head2) * iENU_Work_post_in_ixblue_Ref[:,0],"+",label="E")


PtsStk = []
for head,xxx,yyy in zip(Head1,
                        iENU_Work_post_in_ixblue_Ref[:,0],
                        iENU_Work_post_in_ixblue_Ref[:,1]):
    
    M   = geok.C_2D(head,"rad")
    P   = np.array([xxx,yyy])
    P2 = M.dot(P)
    PtsStk.append(P2)
                      
PtsStk = np.vstack(PtsStk)

plt.plot(TTTref,PtsStk[:,0])
plt.plot(TTTref,PtsStk[:,1])


0.92654171*14.39733338 + 0.37619205*5.84556776


0.37619205*14.39733338 + -0.92654171*5.84556776





