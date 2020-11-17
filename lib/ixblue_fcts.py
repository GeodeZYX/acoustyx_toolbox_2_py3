# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 11:20:11 2020

@author: psakicki
"""

#### Import star style
from geodezyx import *                   # Import the GeodeZYX modules
from geodezyx.externlib import *         # Import the external modules
from geodezyx.megalib.megalib import *   # Import the legacy modules names

#### geodeZYX modules
from geodezyx import conv
from geodezyx import utils

#### Import extern modules
from io import StringIO
from scipy.spatial.transform import Rotation
import itertools
import numpy as np
import os
from matplotlib import pyplot as plt
import pandas as pd

def read_ixblue_data(file_in_path,message_prefix,regex=False):
    """
    Read a iXBlue file and get the data in a Pandas DataFrame

    Parameters
    ----------
    file_in_path : str
        path of the input file.
    message_prefix : str
        descriptor of the iXBlue message, e.g. 'PIXOG,PPC,DETEC', '^\$PTSAY.*' or '^\$PTSAG.*'
        c.f. iXBlue doc
    regex : bool, optional
        If yes, message_prefix is a regular expression like ^\$PTSAY.* or ^\$PTSAG.*
        The default is False.

    Returns
    -------
    DFout : pandas DataFrame
        a RAW DataFrame containing the data. Raw means that it looks like the iXBlue file
    """

    Greped_lines = utils.grep(file_in_path, message_prefix,regex=regex)
    GrepBloc = "".join(Greped_lines)
    DFout = pd.read_csv(StringIO(GrepBloc),header=None,sep=",")
    return DFout


def read_ixblue_TAT(file_in_path):
    """
    Read a iXBlue file and get the Turn Around Time in a dictionnary

    Parameters
    ----------
    file_in_path : str
        path of the input file.

    Returns
    -------
    TATdict : dict
        TAT dico, key is the beacon ID, value is the TAT in ms.

    """
    
    string_PIXOG_tat = "PIXOG,CONFIG,BAL"
    try:
        Greped_lines = utils.grep(file_in_path, string_PIXOG_tat)
        GrepBloc = "".join(Greped_lines)
        DFgaps_tat = pd.read_csv(StringIO(GrepBloc),header=None,sep=",")
        DFgaps_tat = DFgaps_tat[[3,10]]
        print(DFgaps_tat)
        
        # creation of a TAT Dict from a Serie from a DataFrame
        TATdict = DFgaps_tat[[3,10]].set_index(3)[10].to_dict()
    except:
        print("WARN: TAT at 0 !!!!!")
        for ii in DFgaps[7].unique():
            TATdict = {ii : 0}
            
    return TATdict
    

########### POSI iXBLUE
def read_ixblue_realtime_posi(file_in_path,
                              use_extern_up=False,
                              DF_external_U=None):
    """
    Parameters
    ----------
    file_in_path : TYPE
        DESCRIPTION.
    use_extern_up : bool
        If activated, use external Up values. 
    DF_external_U : 4 columns dataframe (columns name date,F,L,&U), optional
        Data for use_extern_up option The default is None.
        Only Up will be used

    Returns
    -------
    DFgaps_posi_ixblue_out : TYPE
        DESCRIPTION.

    """
    
    if use_extern_up:
        FLU  = DF_external_U[["F","L","U"]].values
        TimeInterp = DF_external_U["date"].values
        
        _,_,IU = interpolator_position(FLU,TimeInterp,
                              one_interpolator_per_component=True,
                              bounds_error=False)
    
    string_PTSAG         = "^\$PTSAG"
    DFgaps_posi_ixblue   = read_ixblue_data(file_in_path, string_PTSAG,True)
    DFgaps_posi_ixblue_tmp = pd.DataFrame()
    
    Time = DFgaps_posi_ixblue.apply(ixblue_time_conv,axis=1,args=(2,))
    DFgaps_posi_ixblue_tmp["ID"]    = DFgaps_posi_ixblue[6]
    DFgaps_posi_ixblue_tmp["date"]  = Time
    DFgaps_posi_ixblue_tmp["lat"]   = DFgaps_posi_ixblue.apply(ixblue_coord_conv,axis=1,args=(7,))
    DFgaps_posi_ixblue_tmp["lon"]   = DFgaps_posi_ixblue.apply(ixblue_coord_conv,axis=1,args=(9,))
    
    if use_extern_up:
        DFgaps_posi_ixblue_tmp["depth"] = IU(Time)
    else:
        DFgaps_posi_ixblue_tmp["depth"] = -1. * DFgaps_posi_ixblue[12]
    
    #### ID 0 is the 0 of the GAPS
    DFgaps_posi_ixblue_out = DFgaps_posi_ixblue_tmp[DFgaps_posi_ixblue_tmp.ID == 0].copy()
    
    XYZixblue = conv.GEO2XYZ_vector(DFgaps_posi_ixblue_out[["lat","lon","depth"]])
    DFgaps_posi_ixblue_out["X"],DFgaps_posi_ixblue_out["Y"],DFgaps_posi_ixblue_out["Z"] = XYZixblue[:,0],XYZixblue[:,1],XYZixblue[:,2]
    
    return DFgaps_posi_ixblue_out





### Internal fct for time conversion
def ixblue_time_conv(x,hms_index=3):
    """
    Convert the time from a iXBlue Data Frame row
    ----------
    x : TYPE
        A DataFrame row.
    hms_index : int, optional
        the index of the Hour Min Sec field. The default is 3.

    Returns
    -------
    datetime
        the converted date.
    """
    
    d  = x[hms_index+1]
    mo = x[hms_index+2]
    y  = x[hms_index+3]
    
    hmsdec = str(x[hms_index])
    sdec   = hmsdec.split(".")[-1]
    hms    = hmsdec.split(".")[0]   
    
    s  = float(hms[-2:] + "." + sdec)
    m  = int(hms[-4:-2])
    h  = int(hms[-6:-4])
    
    return conv.datetime_improved(y,mo,d,h,m,s)


    
### Internal fct for coordinate conversion
def ixblue_coord_conv(x,index=9,debug=False):
    """
    Convert a coordinate geographic component from a iXBlue DataFrame row
    ----------
    x : TYPE
        A DataFrame row.
    index : int, optional
        the index of the coordinate component. The default is 9.
    debug : bool
        return also the degree and minute values

    Returns
    -------
    coords
        the converted coordinate in a decimal form.
    """
    orig = str(x[index])
    posneg = x[index+1]
    
    if posneg in ("N","E"):
        koef = 1
    else:
        koef = -1

    decpart = orig.split(".")[-1]
    entpart = orig.split(".")[0]  
      
    try:
        m = float(entpart[-2:] + "." + decpart)
        d = int(entpart[-4:-2])
    except:
        m = np.nan
        d = np.nan
    
    if not debug:
        return koef * conv.dms2dec_num(d,m)
    else:
        return koef * conv.dms2dec_num(d,m) , d , m 
        

def interpolator_position(XYZ,Time,
                          one_interpolator_per_component=False,
                          bounds_error=False):
    """
    A Frontend function to get the coordinate interpolators.
    
    Manage DateTime !!! :)

    Parameters
    ----------
    XYZ : float 3-iterable (3xN)
        the known coordinates.
    Time : datetime iterable
        the known time.
    one_interpolator_per_component : boolean, optional
        If False will generate one interpolator for the 3 components.
        If True will generate 3 interpolators
        The default is False.
    bounds_error : bool, optional
        if False, gives NaN if the interpolated value is outside the range.
        The default is False.

    Returns
    -------
    Iposi : Scipy Interpolator
        Coordinates interpolator.

    """
    
    XYZ = utils.transpose_vector_array(XYZ)
    
    if one_interpolator_per_component:
        IX = conv.interp1d_time(Time,XYZ[0],bounds_error=bounds_error)        
        IY = conv.interp1d_time(Time,XYZ[1],bounds_error=bounds_error)        
        IZ = conv.interp1d_time(Time,XYZ[2],bounds_error=bounds_error)     
        Iposi = (IX,IY,IZ)
    else:
        Iposi = conv.interp1d_time(Time,XYZ,bounds_error=bounds_error)
    return Iposi


def interpolator_attitude(Head,Pitch,Roll,Time,
                          euler_order="zyx",
                          degrees=False):
    """
    A Frontend function to get the attitude interpolators
    
    Manage DateTime !!! :)

    Parameters
    ----------
    Head : float iterable
        Head.
    Pitch : float iterable
        Pitch.
    Roll : float iterable
        Roll.
    Time : datetime iterable
        Time.
    euler_order : str, optional
        c.f. scipy rotation doc. The default is "zyx".
    degrees : bool, optional
        input angles in rad/deg. The default is False.

    Returns
    -------
    Iattitude : Scipy Interpolator
        Attitude interpolator.

    """

    ### We switch to Head Pitch Roll (zyx) ######
    Euler = np.column_stack((Head,Pitch,Roll))

    # Drop duplicates
    DF = pd.DataFrame((np.array(Time),
                       Euler[:,0],
                       Euler[:,1],
                       Euler[:,2])).T

    DF.drop_duplicates(inplace = True)
    Time  = DF[0].values
    Euler = DF[[1,2,3]].values

    # We sort to avoid any conflict
    Time,Euler = utils.sort_binom_list(Time,Euler)
    
    # We generate the interpolator
    Rot_for_slerp = Rotation.from_euler(euler_order,Euler,degrees=degrees)

    Iattitude = conv.Slerp_time(Time,Rot_for_slerp)
    
    return Iattitude


def iXblueDF_prepro_step1B_basic_AHD_BEA_vector(DFgaps_ptsa_in,
                                                normalize=True,
                                                inp_type = "PTSAY"):
    """
    ADD PTSAY/X infos in a iXBlue DataFrame

    Parameters
    ----------
    DFgaps_ptsa_in : TYPE
        DESCRIPTION.
    normalize : TYPE, optional
        DESCRIPTION. The default is True.
    inp_type : TYPE, optional
        if PTSAY => vN,vE,vD.
        if PTSAX => vR,vP,vY.
        The default is "PTSAY".

    Raises
    ------
    Exception
        If inp_type is wrong.

    Returns
    -------
    DFgaps_ptsa_out : TYPE
        DataFrame enriched with PTSAY/X.
        
    Notes
    -----
    PTSAX:
        xxxxx.x X coordinate (+ forward) in meters
        yyyyy.y Y coordinate (+ starboard) in meters
        oooo.oo Calculated depth in meters
    PTSAY:
        xxxxx.x X coordinates (positive northwards) in meters
        yyyyy.y Y coordinates (positive eastwards) in meters
        oooo.o Calculated depth in meters

    """
    
    
    DFgaps_ptsa_out = pd.DataFrame()
    T = DFgaps_ptsa_in.apply(ixblue_time_conv,axis=1,args=(2,))
    DFgaps_ptsa_out["date_rec"] = T
    
    
    if inp_type == "PTSAY":
        cmpstr = ["vN","vE","vD"]
    elif inp_type == "PTSAX":
        cmpstr = ["vR","vP","vY"]
    else:
        print("ERR:iXblueDF_prepro_step1B_basic_AHD_BEA_vector: check inp_type")
        raise Exception
        
    
    ##### non-normalized
    DFgaps_ptsa_out["ID_BEA"]   = DFgaps_ptsa_in[6]
    DFgaps_ptsa_out[cmpstr[0]]  = DFgaps_ptsa_in[7]
    DFgaps_ptsa_out[cmpstr[1]]  = DFgaps_ptsa_in[8]
    DFgaps_ptsa_out[cmpstr[2]]  = DFgaps_ptsa_in[10]
    
    if normalize:
        N = np.linalg.norm(DFgaps_ptsa_out[[cmpstr[0],
                                            cmpstr[1],
                                            cmpstr[2]]],axis=1)
        DFgaps_ptsa_out[cmpstr[0] + "_nrm"] = DFgaps_ptsa_out[cmpstr[0]] / N
        DFgaps_ptsa_out[cmpstr[1] + "_nrm"] = DFgaps_ptsa_out[cmpstr[1]] / N
        DFgaps_ptsa_out[cmpstr[2] + "_nrm"] = DFgaps_ptsa_out[cmpstr[2]] / N        
    
    return DFgaps_ptsa_out



def iXblueDF_prepro_step1_basic(DFin,TATdict,add_attitude=True,
                                average_transducers=False,
                                add_AHDBEA_vector=True,
                                DF_AHDBEA_vectors_tuple=None,
                                clean_nan=True):
    """
    Generate a userfriendly DataFrame from a Raw DataFrame
    it means:
        * Rename the columns
        * compute the emission epoch

    Parameters
    ----------
    DFin : pandas DataFrame
        DataFrame generated by read_ixblue_data.
    TATdict : dictionnary
        TAT dictionnary generated by read_ixblue_TAT.
    add_attitude : bool, optional
        add attitude in the output DataFrame. The default is True.
    average_transducers : bool, optional
        average the 4 transducers.
        If not, return all values
        
    OPTION Description ARE MISSING

    Returns
    -------
    DFout
        a userfriendly DataFrame.
    """

    DFout = pd.DataFrame()
    
    DFin_orig = DFin.copy()
    
    DFout["date_rec"] = DFin.apply(ixblue_time_conv,axis=1)
    
    # new column for ID, friendly readable
    DFout["ID_BEA"] = DFin[7]
    
    ### Extraction of the 4 TWTT (because 4 acoustic heads on GAPS)
    # null TWTT are replaced with NaN ...
    col_transduct_list       = [17,21,25,29]
    col_transduct_phase_list = [18,22,26,30]

    if add_AHDBEA_vector:
        for DF_AHDBEA_vector in DF_AHDBEA_vectors_tuple:
            
            if len(DF_AHDBEA_vector) == 0:
                continue
        
            roundd = "100L" # 100 milisec
            
            
            if "vN" in DF_AHDBEA_vector.columns:
                print("INFO: iXblueDF_prepro_step1_basic: consider NED dir. vectors")
                col_vect_keys_lst = ['vN',
                                     'vE',
                                     'vD',
                                     'vN_nrm',
                                     'vE_nrm',
                                     'vD_nrm']
                
            elif "vR" in DF_AHDBEA_vector.columns:
                print("INFO: iXblueDF_prepro_step1_basic: consider RPY dir. vectors")
                col_vect_keys_lst = ['vR',
                                     'vP',
                                     'vY',
                                     'vR_nrm',
                                     'vP_nrm',
                                     'vY_nrm']
            else:
                print("ERR: iXblueDF_prepro_step1_basic: check DF_AHDBEA_vector column names")
                raise Exception
                
            
            DFout_tmp_vect = DFout.copy()
            DFout_tmp_vect[col_vect_keys_lst] = np.nan
            
            DF_AHDBEA_vector_tmp = DF_AHDBEA_vector.copy()
    
            Tround1 = DFout_tmp_vect.date_rec.dt.round(roundd)        
            Tround2 = DF_AHDBEA_vector_tmp.date_rec.dt.round(roundd)
            
            ID1 = DFout_tmp_vect['ID_BEA']
            ID2 = DF_AHDBEA_vector_tmp['ID_BEA']
            
            DFout_tmp_vect.set_index([Tround1,ID1],inplace=True)
            DF_AHDBEA_vector_tmp.set_index([Tround2,ID2],inplace=True)
            
            DFout_tmp_vect.sort_index(inplace=True)
            DF_AHDBEA_vector_tmp.sort_index(inplace=True)
            
            I1 = DFout_tmp_vect.index
            I2 = DF_AHDBEA_vector_tmp.index
            
            Iinter = I1.intersection(I2)
            Iinter = Iinter.sort_values()
            
            print("INFO: Indexes len for vect assignation", len(I1),len(I2),len(Iinter))
            
            DFout_tmp_vect.loc[Iinter,col_vect_keys_lst] = \
                DF_AHDBEA_vector_tmp.loc[Iinter,col_vect_keys_lst].values.astype(np.float64)
                    
            DFout = DFout_tmp_vect.copy()
        

    DF_TWTT_nan  = DFin[col_transduct_list].replace(0,np.nan)
    DF_Phase_nan = DFin[col_transduct_phase_list].replace(0,np.nan)
    

    
    if average_transducers:     # ... and TWTT averaged with nanmean
        DFout["TWTT_raw"]   = DF_TWTT_nan.apply(np.nanmean,axis=1)
        DFout["Phase_raw"] = DF_Phase_nan.apply(np.nanmean,axis=1)

    else: # Give all the transducer values separatly
        n = len(col_transduct_list)
        
        # heare we duplicate n times the DataFrame
        DFout_tmp_tdc = pd.DataFrame(np.repeat(DFout.values,n,axis=0))
        DFout_tmp_tdc.columns = DFout.columns
        
        DFout_tmp_tdc['TWTT_raw'] = DF_TWTT_nan.values.flatten()
        DFout_tmp_tdc['Phase_raw'] = DF_Phase_nan.values.flatten()
        
        Id_tdc = np.tile(np.arange(1,len(col_transduct_list)+1),len(DFin))

        DFout_tmp_tdc['ID_TDC'] = Id_tdc
        DFout = DFout_tmp_tdc
        
        #if independant transducers, DFin must be extended too
        DFin = pd.DataFrame(np.repeat(DFin.values,n,axis=0))
        
    ##### new column for ID, friendly readable DEFINED ABOVE
    ####DFout["ID_BEA"] = DFin[7]
        
    # # let's keep the date rec first...
    DFout["date_rec_raw"] = DFout["date_rec"].values 
    
    # because we need to correct them based on the 1st ping received
    DFGRP = DFout.groupby("date_rec")
    
    TWTT_raw_min = DFGRP.min()["TWTT_raw"]
    TWTT_raw_min = TWTT_raw_min[DFout["date_rec"]]
    TWTT_raw_min.index = DFout.index
    TWTT_raw_min = pd.to_timedelta(TWTT_raw_min,'us')
    TWTT_raw     = pd.to_timedelta(DFout["TWTT_raw"],'us')
    
    # ### rec is the original rec corrected with the 1st ping reception
    DFout["date_rec"] = DFout["date_rec"] + TWTT_raw - TWTT_raw_min
    # get the emission epoch
    DFout["date_emi"] = DFout["date_rec"] - TWTT_raw
    
    if add_attitude:
        # get attitude
        DFout["head_rec"] =   DFin[11].values
        DFout["roll_rec"] =   DFin[12].values
        DFout["pitc_rec"] = - DFin[13].values
        
            
    
    def internal_twtt_corr(x):
        return x['TWTT_raw'] - 1000. * TATdict[x['ID_BEA']]
    DFout['TWTT'] = DFout.apply(internal_twtt_corr,axis=1)
    
    #remove NaN vals
    if clean_nan:
        DFout.dropna(inplace=True)
        
    return DFout

def iXblueDF_prepro_step2_position(DFin,
                                   InterpolatorPosition,
                                   component_name="XYZ",
                                   component_suffix="GPS",
                                   recemi_tup = ("emi","rec")):
    """
    Interpolate coordinates for the emission and reception epochs
    
    Parameters
    ----------
    DFin : pandas DataFrame
        DataFrame generated at step 1.
    InterpolatorPosition : Scipy Interpolator
        interpolator generated with interpolator_position.
    component_name : str, optional
        coordinate components considered ("XYZ", "ENU", "NED"). 
        The default is "XYZ".
    component_suffix : str, optional
        DESCRIPTION. The default is "GPS".
    recemi_tup : tuple, optional
        The default is ("emi","rec").

    Returns
    -------
    DFout : pandas DataFrame
            DataFrame with interpolated coordinates for
            the emission and reception epochs.
            columns with have the name <component_name>_<component_suffix>_recemi
            e.g. X_GPS_emi
    """
    
    

    DFout = DFin.copy(deep=True)
    
    for recemi in recemi_tup:
        suffix = "_".join(("",component_suffix,recemi))
        PosiInterpolated = InterpolatorPosition(DFout["_".join(("date",recemi))])
        DFout[component_name[0]+suffix] = PosiInterpolated[0]
        DFout[component_name[1]+suffix] = PosiInterpolated[1]
        DFout[component_name[2]+suffix] = PosiInterpolated[2]
    
    return DFout


def ixblueDF_prepro_step3_down_component(DFin):
    """
    Generate Down component columns 

    Parameters
    ----------
    DFin : pandas DataFrame
        DataFrame generated at step 2.

    Returns
    -------
    DFout : pandas DataFrame
    """
    
    DFout = DFin.copy(deep=True)
    Up_col = [e for e in DFout.columns if e.startswith("U_")]
    for col in Up_col:
        DFout["D" + col[1:]] = - DFout[col]
    return DFout

def ixblueDF_prepro_step4_attitude(DFin,InterpolatorAttitude,
                                   lever_arm_RPY_dic,
                                   ref_dev = "GPS",
                                   dev_tup = ("AHD",),
                                   emirec_tup = ('emi','rec'),
                                   coord_sys = "NED"):
    """
    Interpolate attitude for the emission and reception epoch
    And apply the lever arms

    Parameters
    ----------
    DFin : pandas DataFrame
        DataFrame generated at step 3.
    InterpolatorAttitude : Scipy Interpolator
        Attitude interpolator generated with interpolator_attitude.
    lever_arm_RPY_dic : dict
        a dictionnary containg the lever arm ties.
        key is the device (GPS, AHD for Acoustic Head ...)
        value is a 3D array with the ties in the given coord_sys
    ref_dev : str, optional
        the device witch gives the reference position. The default is "GPS".
    dev_tup : str, optional
        the devices on which the position will be transfered.
        The default is ("AHD",).
    emirec_tup : tuple, optional
        The default is ('emi','rec').
    coord_sys : str, optional
        the coordinate system ("NED", "ENU" ...). The default is "NED".
        
    Note
    ----
    
    We do: NED_device = NED - iNED_reference (GPS) + iNED_device

    Returns
    -------
    DFout : pandas DataFrame

    """
    DFout = DFin.copy(deep=True)
        
    for emirec,dev in itertools.product(emirec_tup,dev_tup):        
        Rot  = InterpolatorAttitude(DFout['date_' + emirec])
        ######## !!!!!!!!!!!!!!!!!!!!!!!
        ## Dirty replacement
        ##iABC_REF = lever_arm_RPY_dic[ref_dev]
        ##iABC_DEV = lever_arm_RPY_dic[dev]    
        ######## !!!!!!!!!!!!!!!!!!!!!!!
        iABC_REF = Rot.apply(lever_arm_RPY_dic[ref_dev])
        iABC_DEV = Rot.apply(lever_arm_RPY_dic[dev])
    
        col_name_REF = ["_".join((coord_sys[0],ref_dev,emirec)),
                        "_".join((coord_sys[1],ref_dev,emirec)),
                        "_".join((coord_sys[2],ref_dev,emirec))]
        
        col_name_DEV = ["_".join((coord_sys[0],dev,emirec)),
                        "_".join((coord_sys[1],dev,emirec)),
                        "_".join((coord_sys[2],dev,emirec))]
        
        ABC = DFout[col_name_REF].values
        ABC_DEV = ABC - iABC_REF + iABC_DEV
        
        #print("*******************************")
        #print(emirec,dev)
        #print(ABC,ABC_DEV,- iABC_REF + iABC_DEV)
        #print("*******************************")
            
        DFout[col_name_DEV[0]],DFout[col_name_DEV[1]],DFout[col_name_DEV[2]] = ABC_DEV.T
        
    return DFout


def ixblueDF_prepro_step4b_ixblue_rt_posi_manager(DFin):
    DFout = DFin.copy(deep=True)

    DFout["N_AHD_emi"] = DFin["N_GPS_emi"]
    DFout["E_AHD_emi"] = DFin["E_GPS_emi"]
    DFout["U_AHD_emi"] = DFin["U_GPS_emi"]
    DFout["D_AHD_emi"] = DFin["D_GPS_emi"]
    
    DFout["N_AHD_rec"] = DFin["N_GPS_rec"]
    DFout["E_AHD_rec"] = DFin["E_GPS_rec"]
    DFout["U_AHD_rec"] = DFin["U_GPS_rec"]
    DFout["D_AHD_rec"] = DFin["D_GPS_rec"]  
    
    DFout["N_TDC_rec"] = DFin["N_GPS_rec"]
    DFout["E_TDC_rec"] = DFin["E_GPS_rec"]
    DFout["U_TDC_rec"] = DFin["U_GPS_rec"]
    DFout["D_TDC_rec"] = DFin["D_GPS_rec"]      
    return DFout

def ixblueDF_prepro_step4c_TDC_posi(DFin):
    
    DFout = DFin.copy(deep=True)

    for emirec in ('_emi','_rec'):
        Valstk = []
        for irow , row in DFout.iterrows():            
            tdcstr = "TDC" + str(row["ID_TDC"])
            cols_tdc_in   = ['N_' + tdcstr + emirec,
                             'E_' + tdcstr + emirec,
                             'D_' + tdcstr + emirec]
            Valstk.append(row[cols_tdc_in])
            
        cols_tdc_out =  ['N_' + "TDC" + emirec,
                         'E_' + "TDC" + emirec,
                         'D_' + "TDC" + emirec]
        
        DFout[cols_tdc_out] = np.vstack(Valstk)

    return DFout


def ixblueDF_prepro_step5_Ofile_export(DFin,TATdict,
                                   export_dir,
                                   experience_name,
                                   apriori_coordinates):
    """
    Export the Dataframe to a O-file
    
    Parameters
    ----------
    DFin : pandas DataFrame
        DataFrame generated at step 5.
    TATdict : dict
        TAT dictionnary generated by read_ixblue_TAT.
    export_dir : str
        export directory.
    experience_name : str
        name of the experiment.
    apriori_coordinates : dict
        key = Beacon ID, value a 3D array

    Returns
    -------
    ofile_path : str
        path of the exported O-file.

    """

    DFout = DFin.copy(deep=True)    
    DFout = DFout.dropna()


    for id_beacon in DFout['ID_BEA'].unique():
        
        DFout_beacon = DFout[DFout['ID_BEA'] == id_beacon]
        DFout_beacon.reset_index(inplace=True)
        
        ############## HARDCODED CLEANING
        #DFout_beacon = DFout_beacon[DFout_beacon["TWTT"] < 70000.]
        ##############
        
        DFout_Ofile = pd.concat((
        DFout_beacon["date_emi"].apply(conv.dt2posix),
        DFout_beacon["E_AHD_emi"],
        DFout_beacon["N_AHD_emi"],
        DFout_beacon["D_AHD_emi"],
        DFout_beacon["date_emi"].dt.year,
        DFout_beacon["date_emi"].dt.month,
        DFout_beacon["date_emi"].dt.day,
        DFout_beacon["date_emi"].dt.hour,
        DFout_beacon["date_emi"].dt.minute,
        DFout_beacon["date_emi"].dt.second,
        DFout_beacon["date_emi"].dt.microsecond,
        
        DFout_beacon["date_rec"].apply(conv.dt2posix),
        DFout_beacon["E_AHD_rec"],
        DFout_beacon["N_AHD_rec"],
        DFout_beacon["D_AHD_rec"],
        DFout_beacon["date_rec"].dt.year,
        DFout_beacon["date_rec"].dt.month,
        DFout_beacon["date_rec"].dt.day,
        DFout_beacon["date_rec"].dt.hour,
        DFout_beacon["date_rec"].dt.minute,
        DFout_beacon["date_rec"].dt.second,
        DFout_beacon["date_rec"].dt.microsecond,
        
        DFout_beacon["TWTT_raw"] * 10**-6,
        pd.Series([TATdict[id_beacon]* 10**-3] * len(DFout_beacon["date_emi"]),name='TAT'),
        DFout_beacon["TWTT"] * 10**-6),
        axis=1)
            
        ofile_path = os.path.join(export_dir,
                                  experience_name+".PXP"+str(id_beacon)+".O.dat")
    
        header = 'pxp_coords : ' + str(apriori_coordinates[id_beacon])
    
        np.savetxt(ofile_path,DFout_Ofile,header=header)

    return ofile_path
    
