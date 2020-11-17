#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# ## Import of the libraires

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
## Barycenter
p_gaps_10 = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/01_DATA/07-2019_IXBLUE_BREST/0_RAW_DATA/data-gaps/data - jeudi 25 juillet 2019/10 - Repeater Gaps - 14h10 - au barycentre IIS.txt"
p_gaps_11 = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/01_DATA/07-2019_IXBLUE_BREST/0_RAW_DATA/data-gaps/data - jeudi 25 juillet 2019/11 - Repeater Gaps - 14h46 - au barycentre - CIS.txt"

### Boxin Beacon 1
p_gaps_01 = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/01_DATA/07-2019_IXBLUE_BREST/0_RAW_DATA/data-gaps/data - jeudi 25 juillet 2019/04 - Repeater Gaps - 11h58 - box in TP1.txt"
### Boxin Beacon 2
p_gaps_02 = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/01_DATA/07-2019_IXBLUE_BREST/0_RAW_DATA/data-gaps/data - jeudi 25 juillet 2019/02 - Repeater Gaps - 11h09 - box in TP2.txt"
### Boxin Beacon 3
p_gaps_03 = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/01_DATA/07-2019_IXBLUE_BREST/0_RAW_DATA/data-gaps/data - jeudi 25 juillet 2019/06 - Repeater Gaps - 12h45 - box in TP3.txt"

#### Position files
p_posi_02  = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/01_DATA/07-2019_IXBLUE_BREST/4_RESULTS/2019-07-25/Manip_5_1/2019_07_25_Manip5_1_gins_ppp.csv"
p_posi_02b = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/01_DATA/07-2019_IXBLUE_BREST/4_RESULTS/2019-07-25/Manip_5_1/2019_07_25_Manip5_1_nrcan.csv"
p_posi_02c = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/01_DATA/07-2019_IXBLUE_BREST/4_RESULTS/2019-07-25/2019_07_25_Manip5_123_nrcan.csv"
p_posi_02d = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/01_DATA/07-2019_IXBLUE_BREST/4_RESULTS/2019-07-25/2019_07_25_Manip5_123_rtklib_ginscoords.csv"
## RTKlib diff
p_posi_02e = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/01_DATA/21-08-2020_ResultsRTKLib/pameli_5_GNSS.csv"

p_posi_03  = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/01_DATA/07-2019_IXBLUE_BREST/4_RESULTS/2019-07-25/Manip_4_1/2019_07_25_Manip4_1_gins_ppp.csv"
p_posi_03b = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/01_DATA/07-2019_IXBLUE_BREST/4_RESULTS/2019-07-25/Manip_4_1/2019_07_25_Manip4_1_nrcan.csv"
p_posi_03c = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/01_DATA/07-2019_IXBLUE_BREST/4_RESULTS/2019-07-25/2019_07_25_Manip4_123_nrcan.csv"
p_posi_03d = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/01_DATA/07-2019_IXBLUE_BREST/4_RESULTS/2019-07-25/2019_07_25_Manip4_123_rtklib_ginscoords.csv"
## RTKlib diff
p_posi_03e = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/01_DATA/21-08-2020_ResultsRTKLib/pameli_4_GNSS.csv"

p_posi_01  = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/01_DATA/07-2019_IXBLUE_BREST/4_RESULTS/2019-07-25/Manip_3_1/2019_07_25_Manip3_1_gins_ppp.csv"
p_posi_01b = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/01_DATA/07-2019_IXBLUE_BREST/4_RESULTS/2019-07-25/Manip_3_1/2019_07_25_Manip3_1_nrcan.csv"
p_posi_01c = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/01_DATA/07-2019_IXBLUE_BREST/4_RESULTS/2019-07-25/2019_07_25_Manip3_123_nrcan.csv"
p_posi_01d = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/01_DATA/07-2019_IXBLUE_BREST/4_RESULTS/2019-07-25/2019_07_25_Manip3_123_rtklib_ginscoords.csv"
## RTKlib diff
p_posi_01e = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/01_DATA/21-08-2020_ResultsRTKLib/pameli_3_GNSS.csv"

p_posi_11  = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/01_DATA/07-2019_IXBLUE_BREST/4_RESULTS/2019-07-25/Manip_2_2/2019_07_25_Manip2_2_gins_ppp.csv"
p_posi_11b = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/01_DATA/07-2019_IXBLUE_BREST/4_RESULTS/2019-07-25/Manip_2_2/2019_07_25_Manip2_2_rtklib_nrcancoord.csv"

p_gaps_tat = p_gaps_11

#### Boxin 
### Boxin Beacon 2
p_gaps = p_gaps_02
p_posi = p_posi_02c

### Boxin Beacon 3
idBea = 3
p_gaps = p_gaps_03
p_posi = p_posi_03d
p_posi = p_posi_03c
p_posi = p_posi_03e

### Boxin Beacon 1
idBea = 1
p_gaps = p_gaps_01
p_posi = p_posi_01c
p_posi = p_posi_01d
p_posi = p_posi_01e

### Boxin Beacon 2
idBea = 2 
p_gaps = p_gaps_02
p_posi = p_posi_02c
p_posi = p_posi_02e

#### Barycentre
p_gaps = p_gaps_11
p_posi = p_posi_11b


# ## Read the position

# ### Import

# In[5]:


DFposi_post = pd.read_csv(p_posi,header=0,sep=",")

if 'Unnamed: 0' in DFposi_post.columns:
    DFposi_post = DFposi_post.rename(columns={'Unnamed: 0':'date'})
    
# ### Conversion to XYZ

# In[6]:
    
########### POSI EXTE
if not "X" in DFposi_post.columns:
    XYZpost = conv.GEO2XYZ_vector(DFposi_post[["lat","lon","values"]])
    DFposi_post["X"],DFposi_post["Y"],DFposi_post["Z"] = XYZpost[:,0],XYZpost[:,1],XYZpost[:,2]


DFposi_post["date"] = pd.to_datetime(DFposi_post['date']) - dt.timedelta(seconds=18)


if False:
    print("WARN POST PROCESSING IN RGF93")
    HPparam = reffram.itrf_helmert_get_parameters("ITRF2014","ETRF2000",
                                              verbose=False,
                                              convert=True)
    XYZ_post_RGF93 = reffram.itrf_speed_calc(DFposi_post[["X"]],DFposi_post[["Y"]],DFposi_post[["Z"]],
                                             2019.6,
                                             -1.16615726390990e-02,
                                             1.76696025246000e-02,
                                             1.08543846164409e-02,
                                             2009.0)
    XYZ_post_RGF93 = np.column_stack(XYZ_post_RGF93)
    XYZ_post_RGF93 = reffram.itrf_helmert_trans(XYZ_post_RGF93,2009.,*HPparam)
    DFposi_post[["X","Y","Z"]] = XYZ_post_RGF93
    

xyz_ref = np.mean(DFposi_post[["X","Y","Z"]].values,axis=0)
xyz_ref = np.array([4236447.54165442, -329431.84160855, 4740657.55330567])
##### Posi central taken from the barycenter
xyz_ref = np.array([4236463.24737245, -329434.98845412, 4740641.89446747])



########### POSI iXBLUE
string_PTSAG         = "^\$PTSAG"
DFgaps_posi_ixblue   = ixblue_fcts.read_ixblue_data(p_gaps, string_PTSAG,True)
DFgaps_posi_ixblue_3 = pd.DataFrame()

TTT = DFgaps_posi_ixblue.apply(ixblue_fcts.ixblue_time_conv,axis=1,args=(2,))
DFgaps_posi_ixblue_3["ID"]    = DFgaps_posi_ixblue[6]
DFgaps_posi_ixblue_3["date"]  = TTT
DFgaps_posi_ixblue_3["lat"]   = DFgaps_posi_ixblue.apply(ixblue_fcts.ixblue_coord_conv,axis=1,args=(7,))
DFgaps_posi_ixblue_3["lon"]   = DFgaps_posi_ixblue.apply(ixblue_fcts.ixblue_coord_conv,axis=1,args=(9,))
DFgaps_posi_ixblue_3["depth"] = 0 #DFgaps_ptsa_2[12]

DFgaps_posi_ixblue_2 = DFgaps_posi_ixblue_3[DFgaps_posi_ixblue_3.ID == 0].copy()

XYZixblue = conv.GEO2XYZ_vector(DFgaps_posi_ixblue_2[["lat","lon","depth"]])
DFgaps_posi_ixblue_2["X"],DFgaps_posi_ixblue_2["Y"],DFgaps_posi_ixblue_2["Z"] = XYZixblue[:,0],XYZixblue[:,1],XYZixblue[:,2]

#############
ENU = conv.XYZ2ENU_vector(DFgaps_posi_ixblue_2[["X","Y","Z"]], xyz_ref)
DFgaps_posi_ixblue_2["E"],DFgaps_posi_ixblue_2["N"],DFgaps_posi_ixblue_2["U"] = ENU[:,0],ENU[:,1],ENU[:,2]

ENU = conv.XYZ2ENU_vector(DFposi_post[["X","Y","Z"]], xyz_ref)
DFposi_post["E"],DFposi_post["N"],DFposi_post["U"] = ENU[:,0],ENU[:,1],ENU[:,2]


############# CHOIX POSI IXBLUE OR EXTERN



DFposi = DFposi_post
bool_lever_arm_correction = True

DFposi = DFgaps_posi_ixblue_2
bool_lever_arm_correction = False

DFposi = DFposi_post
bool_lever_arm_correction = False

estimate_ties = 0

exp_title = "DFposi_post"
DFposi = DFposi_post
bool_lever_arm_correction = False

exp_title = "DFgaps_posi_ixblue_2"
DFposi = DFgaps_posi_ixblue_2
bool_lever_arm_correction = False

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
plt.plot(DFposi_post["date"],DFposi_post["E"],label="postpro")
plt.plot(DFgaps_posi_ixblue_2["date"],DFgaps_posi_ixblue_2["E"],label="RTK")
plt.title("East")
plt.legend()


plt.figure()
plt.plot(DFposi_post["date"],DFposi_post["N"],label="postpro")
plt.plot(DFgaps_posi_ixblue_2["date"],DFgaps_posi_ixblue_2["N"],label="RTK")
plt.title("Nord")
plt.legend()

plt.figure()
plt.plot(DFposi_post["E"],DFposi_post["N"],"x",label="Posi Post")
plt.plot(DFgaps_posi_ixblue_2["E"],DFgaps_posi_ixblue_2["N"],"x",label="Posi iXblue")
plt.axis("equal")


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
    


DFgaps_ptsa_4 = ixblue_fcts.iXblueDF_prepro_step1B_basic_AHD_BEA_vector(DFgaps_ptsa)

from pygoat.seafloorpos import ixblue_fcts

DFixblue1 = ixblue_fcts.iXblueDF_prepro_step1_basic(DFgaps_raw,TATdict,
                                                    average_transducers=False,
                                                    add_AHDBEA_vector=1,
                                                    DF_AHDBEA_vector=DFgaps_ptsa_4,
                                                    clean_nan = 1)

#%%

################## Merge vectors in the big dataframe

# from pandas.tseries.offsets import DateOffset  

# roundid = '100L'

# DFgaps_raw_bis = DFgaps_raw.copy()

# Tround1 = DFixblue1.date_rec.dt.round(roundid)
# DFgaps_raw_bis['ID'] = DFgaps_raw[7]

# DFixblue1_bis = DFixblue1.copy()




# DFgaps_ptsa_3_bis = DFgaps_ptsa_3.copy()

# Tround2 = DFgaps_ptsa_3_bis.date_rec.dt.round(roundid)

# DFgaps_raw_bis.set_index(['TT','ID'],inplace=True)
# DFgaps_ptsa_3_bis.set_index([Tround2,'ID_BEA'],inplace=True)

# DFgaps_raw_bis.sort_index(inplace=True)
# DFgaps_ptsa_3_bis.sort_index(inplace=True)

# I1 = DFgaps_raw_bis.index
# I2 = DFgaps_ptsa_3_bis.index

# Iinter = I1.intersection(I2)
# Iinter = Iinter.sort_values()

# DFgaps_raw_bis.loc[Iinter]
# DFgaps_ptsa_3_bis.loc[Iinter]



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

# plt.plot(DFixblue2['date_rec'],DFixblue2['head_rec'])




# ### Define the lever arms

# ##### Naming convention
# - GPS: GPS Antenna
# - AHD: Acoustic Head 
# - TDCx: Transducer x

# In[16]:

lever_arm_RPY_dic_orig   = dict()
lever_arm_RPY_dic_custom = dict()


RESstk = []

#for TUP in itertools.product(np.arange(-1.4,-1.2,0.01),np.arange(-0.2,0.2,0.01)):
# for TUP in itertools.product(np.arange(-1.5,-1.2,0.005),np.arange(-0.5,-0.2,0.005)):
# for TUP in itertools.product(np.arange(-1.4,-1.2,0.001),np.arange(-0.2,0.2,0.001)):

for TUP in [(1,1)]:    

    rahd , pahd = TUP
    
    rahd = np.round(rahd,4)
    pahd = np.round(pahd,4)
       
    # !!!!! CHECK LEVER ARMS and DIRECTIONs !!!!!
    # -0.039 m / +0.003 m / +1.481 m
    lever_arm_RPY_dic_orig["GPS"]  = np.array([-0.039,+0.003,-1.481])
    lever_arm_RPY_dic_orig["AHD"]  = np.array([0.,0.,0.3592])
    
    lever_arm_RPY_dic_orig["TDC1"] = np.array([ 0.10735,0.,0.5095])
    lever_arm_RPY_dic_orig["TDC2"] = np.array([-0.10735,0.,0.5095])
    lever_arm_RPY_dic_orig["TDC3"] = np.array([0.,-0.10735,0.5733])
    lever_arm_RPY_dic_orig["TDC4"] = np.array([0., 0.10735,0.5733])
    
    ############ CUSTOM VERSION 
    lever_arm_RPY_dic_custom["AHD"]  = np.array([-1.5,-0.2,0.3592])
    lever_arm_RPY_dic_custom["GPS"]  = np.array([-0.039,+0.003,-1.481])
    lever_arm_RPY_dic_custom["AHD"]  = np.array([-0.66997978,0.30292373,0.])
    lever_arm_RPY_dic_custom["AHD"]  = np.array([rahd,pahd,0.3592])
    lever_arm_RPY_dic_custom["AHD"]  = np.array([-1.5,-0.2,0.3592])
    
    lever_arm_RPY_dic_custom["AHD"]  = np.array([-1.28477607, -0.05694931,0])
    lever_arm_RPY_dic_custom["AHD"]  = np.array([-1.29801996, -0.0539056,0])
    lever_arm_RPY_dic_custom["AHD"]  = np.array([-0.00807833, -0.13654185,0])
    ### Based on STD - RTK
    lever_arm_RPY_dic_custom["AHD"]  = np.array([ 0.03843428, -0.09740861])
    ### Based on D   - RTK
    lever_arm_RPY_dic_custom["AHD"]  = np.array([-0.00807833, -0.13654185,0])
    ### Based on D   - POST
    lever_arm_RPY_dic_custom["AHD"]  = np.array([0.15670585, -0.05443877,0])
    ### Based on STD - POST
    lever_arm_RPY_dic_custom["AHD"]  = np.array([0.11316333, -0.09949707,0])
    
###############################

    lever_arm_RPY_dic_orig["AHD"]  = np.array([0.,0.,0.3592])
    lever_arm_RPY_dic_orig["GPS"]  = np.array([-0.18641225,  0.06124067,0])

    
    lever_arm_RPY_dic_custom["TDC1"] = np.array([ 0.10735,0.,0.5095])
    lever_arm_RPY_dic_custom["TDC2"] = np.array([-0.10735,0.,0.5095])
    lever_arm_RPY_dic_custom["TDC3"] = np.array([0.,-0.10735,0.5733])
    lever_arm_RPY_dic_custom["TDC4"] = np.array([0., 0.10735,0.5733])
    
    
    lever_arm_RPY_dic_null = dict()
    # !!!!! CHECK LEVER ARMS and DIRECTIONs !!!!!
    # -0.039 m / +0.003 m / +1.481 m
    lever_arm_RPY_dic_null["GPS"]  = np.array([0.,0.,0.])
    lever_arm_RPY_dic_null["AHD"]  = np.array([0.,0.,0.])
    
    lever_arm_RPY_dic_null["TDC1"] = np.array([0.,0.,0.])
    lever_arm_RPY_dic_null["TDC2"] = np.array([0.,0.,0.])
    lever_arm_RPY_dic_null["TDC3"] = np.array([0.,0.,0.])
    lever_arm_RPY_dic_null["TDC4"] = np.array([0.,0.,0.])
    
    
    if not bool_lever_arm_correction:
        lever_arm_RPY_dic_opera = lever_arm_RPY_dic_null
        print("WARN NO LEVER ARM CORRECTION")
    else:
        print("CUSTOM LEVER ARMS !!!!!!")
        lever_arm_RPY_dic_opera = lever_arm_RPY_dic_custom
        lever_arm_RPY_dic_opera = lever_arm_RPY_dic_orig
    
    # ### Interpolate the attitude and apply the lever arms
    
    # In[17]:
    
    
    DFgaps_ptsa_5 = ixblue_fcts.ixblueDF_prepro_step4_attitude(DFgaps_ptsa_4,I_Att,
                                                               lever_arm_RPY_dic_opera,
                                                               ref_dev = "GPS",
                                                               dev_tup = ("AHD",),
                                                               emirec_tup = ('emi','rec'),
                                                               coord_sys = "NED")
    DFgaps_ptsa_9 = DFgaps_ptsa_5.dropna().copy()
    
    
    DFgaps_ptsa_5_null = ixblue_fcts.ixblueDF_prepro_step4_attitude(DFgaps_ptsa_4,I_Att,
                                                                    lever_arm_RPY_dic_null,
                                                                    ref_dev = "GPS",
                                                                    dev_tup = ("AHD",),
                                                                    emirec_tup = ('emi','rec'),
                                                                    coord_sys = "NED")
    DFgaps_ptsa_9_null = DFgaps_ptsa_5.dropna().copy()
    
    dNcorr = 0
    dEcorr = 0
    dDcorr = 0
    
    dNcorr = DFgaps_ptsa_9["dN"]
    dEcorr = DFgaps_ptsa_9["dE"]
    dDcorr = DFgaps_ptsa_9["dD"]
    
    dNcorr = 1.*DFgaps_ptsa_9["dN"]
    dEcorr = 1.*DFgaps_ptsa_9["dE"]
    dDcorr = 1.*DFgaps_ptsa_9["dD"]
    
    dNEDcorr = DFgaps_ptsa_9[["dN","dE","dD"]].values
    
    N_AHD_rec = DFgaps_ptsa_9.N_AHD_rec
    E_AHD_rec = DFgaps_ptsa_9.E_AHD_rec
    D_AHD_rec = DFgaps_ptsa_9.D_AHD_rec
    
    Dist1 = np.linalg.norm(np.column_stack((dNcorr,dEcorr)),axis=1)
    Dist2 = np.linalg.norm(np.column_stack((N_AHD_rec,E_AHD_rec)),axis=1)
    
    dDist = Dist1 - Dist2
    
    DFgaps_ptsa_9["N_BEA"] = N_AHD_rec + dNcorr
    DFgaps_ptsa_9["E_BEA"] = E_AHD_rec + dEcorr
    DFgaps_ptsa_9["D_BEA"] = D_AHD_rec + dDcorr
    
    
    def lever_arm_estim_fct(lever_arm,DFgaps_in,out="STD"):
        
        dNEDcorr = DFgaps_in[["dN","dE","dD"]].values
            
        lever_arm_RPY_dic_in = dict()
        #lever_arm_RPY_dic_in["GPS"] = np.array([-0.039,+0.003,-1.481])
        #lever_arm_RPY_dic_in["AHD"] = lever_arm
        #lever_arm_RPY_dic_in["AHD"] = np.array((lever_arm[0],lever_arm[1],0))
        
        lever_arm_RPY_dic_in["GPS"] = np.array((lever_arm[0],lever_arm[1],0))        
        lever_arm_RPY_dic_in["AHD"] = np.array([0.,0.,0.3592])
        
        attitude_fct = ixblue_fcts.ixblueDF_prepro_step4_attitude
        
        DFgaps_ptsa_tmp = attitude_fct(DFgaps_in,I_Att,
                                       lever_arm_RPY_dic_in,
                                       ref_dev    = "GPS",
                                       dev_tup    = ("AHD",),
                                       emirec_tup = ('emi','rec'),
                                       coord_sys  = "NED")
        
        NED_AHD = DFgaps_ptsa_tmp[["N_AHD_rec",
                                   "E_AHD_rec",
                                   "D_AHD_rec"]].values
        
        NED_BEA = NED_AHD + dNEDcorr
        
        NED_BEA_mean = np.array([-2.214199812665709 ,
                                 -1.652074585238479 ,
                                 31.05983569])
        
        NED_BEA_mean = np.mean(NED_BEA,axis=0)

        if out == "D":
            D = np.linalg.norm(NED_BEA - NED_BEA_mean ,axis=1) 
            DD = np.mean(D)            
            OUT = DD
        
        elif out == "STD":            
            Nclean, NcleanBool = stats.outlier_mad(NED_BEA[:,0])
            Eclean, EcleanBool = stats.outlier_mad(NED_BEA[:,1])
            
            N_BEA_std = np.std(Nclean,axis=0)
            E_BEA_std = np.std(Eclean,axis=0)

            OUT = E_BEA_std
            
            OUT = np.sqrt(N_BEA_std**2 + E_BEA_std**2)
            
        print(OUT)
        return OUT
        
    if 0 and estimate_ties:
        SOL = scipy.optimize.least_squares(lever_arm_estim_fct,
                                    np.array([-1.5,-0.2,0.3592]),
                                    args=(DFgaps_ptsa_9_null,))  
    if 1 and estimate_ties:
        ## (-1.4,-1.2,0.001),np.arange(-0.2,0.2,0.001)
        rranges = (slice(-1.4,-1.2,0.1),
                   slice(-0.2,0.2,0.1),
                   slice(-1, 1, 0.5))
        rranges = ((-1.4,-1.2),
                   (-0.2,0.2)) #,(-1, 1))
        rranges =  ((-0.2,0.2),
                   (-0.2,0.2)) #,(-1, 1))
        SOL1     = scipy.optimize.brute(lever_arm_estim_fct,
                                   rranges,args=(DFgaps_ptsa_9_null,))
        
    if 1 and estimate_ties:
        rranges = ((-1.4,-1.2),(-0.2,0.2)) #,(-1, 1))
        rranges =  ((-0.2,0.2),
                   (-0.2,0.2)) #,(-1, 1))
        SOL2 = scipy.optimize.shgo(lever_arm_estim_fct,
                                  rranges,args=(DFgaps_ptsa_9_null,))        
    if 1 and estimate_ties:
        rranges = ((-1.4,-1.2),(-0.2,0.2)) #,(-1, 1))
        rranges =  ((-0.2,0.2),
                   (-0.2,0.2)) #,(-1, 1))
        SOL3 = scipy.optimize.differential_evolution(lever_arm_estim_fct,
                                  rranges,args=(DFgaps_ptsa_9_null,))        
                
    if 0 and estimate_ties:
        rranges = ((-1.4,-1.2),(-0.2,0.2)) #,(-1, 1))
        rranges =  ((-0.2,0.2),
                   (-0.2,0.2)) #,(-1, 1))
        SOL4 = scipy.optimize.dual_annealing(lever_arm_estim_fct,
                                  rranges,args=(DFgaps_ptsa_9_null,))        
    if False:
        plt.figure()
    
    for idbea in DFgaps_ptsa_9["ID_BEA"].unique():    
        DFgaps_ptsa_BEA = DFgaps_ptsa_9[DFgaps_ptsa_9["ID_BEA"] == idbea]
        
        E = DFgaps_ptsa_BEA["E_BEA"].values
        N = DFgaps_ptsa_BEA["N_BEA"].values
        
        print("len",len(DFgaps_ptsa_BEA))
        
        from geodezyx import stats
        
        Eclean, EcleanBool = stats.outlier_mad(E)
        Nclean, NcleanBool = stats.outlier_mad(N)
    
        Emean = np.nanmean(E)
        Estd  = np.nanstd(E)
        Emed  = np.nanmedian(E)
    
        Nmean = np.nanmean(N)
        Nstd  = np.nanstd(N)
        Nmed  = np.nanmedian(N)
        
        DistOutlier = np.linalg.norm(np.column_stack((N - Nmed,E - Emed)),axis=1)
        _, Booloutlier = stats.outlier_mad(DistOutlier)
        
        Eclean = E[Booloutlier]
        Nclean = N[Booloutlier]
        
        Eclean_mean = np.nanmean(Eclean)
        Eclean_std  = np.nanstd(Eclean)
    
        Nclean_mean = np.nanmean(Nclean)
        Nclean_std  = np.nanstd(Nclean)
    
        print("E",idbea,Emean,Estd)
        print("N",idbea,Nmean,Nstd)
        print("Eclean",idbea,Eclean_mean,Eclean_std)
        print("Nclean",idbea,Nclean_mean,Nclean_std)
            
        #plt.plot(E,N,"+")
        
        if True:
            plt.plot(E_AHD_rec,N_AHD_rec,"b+",label="Posi Interpoled")
            plt.quiver(E_AHD_rec,N_AHD_rec,dEcorr,dNcorr,scale_units="xy",scale=1)
            plt.plot(Eclean,Nclean,"y+",label="Posi Beacon")
            plt.xlim((-2,3))
            plt.ylim((-4,-1))
            # utils.join_improved(" ",rahd,pahd)
            plt.suptitle(exp_title)
            plt.title("BEA" + str(idBea) + "LeverArmCorr" + str(bool_lever_arm_correction))
            plt.axis("equal")
            plt.legend()
            print("stdev", np.nanstd(Eclean),np.nanstd(Nclean))
        
        print(rahd,pahd,np.std(Eclean),np.std(Nclean))
        
        RESstk.append((rahd,pahd,np.std(Eclean),np.std(Nclean)))
        
        # plt.figure()
        
        # T = DFgaps_ptsa_9.date_emi
        # V1 = DFgaps_ptsa_9['N_GPS_emi'] - DFgaps_ptsa_9['N_AHD_emi']
        # plt.plot(T,V1,".")
        # V2 = DFgaps_ptsa_9['E_GPS_emi'] - DFgaps_ptsa_9['E_AHD_emi']
        # plt.plot(T,V2,".")
        # V3 = DFgaps_ptsa_9['D_GPS_emi'] - DFgaps_ptsa_9['D_AHD_emi']
        # plt.plot(T,V3,".")
        
        # plt.plot(T,np.sqrt(V1**2 + V2**2 + V3**2),"x")
    
    if False:    
        plt.figure()
        plt.plot(DFposi.date,DFposi.N,"-")
        plt.plot(DFgaps_ptsa_9.date_emi,DFgaps_ptsa_9.N_GPS_emi,".")
    

A = np.vstack(RESstk)
A = pd.DataFrame(A)
print("***** FINAL MIN *****")
print(A[A[2] == A[2].min()])
print(A[A[3] == A[3].min()])

print("EXPERI:", exp_title)
print("BEACON:", idBea)
print("CORRLA:", bool_lever_arm_correction)
    
    
    
    # ## Export the DataFrame as a O-File
    
