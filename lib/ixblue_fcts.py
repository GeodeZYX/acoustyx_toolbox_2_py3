#!/usr/bin/env python3
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

from io import StringIO
from scipy.spatial.transform import Rotation


import itertools
import numpy as np
import os
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation
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
        # creation of a TAT Dict from a Serie from a DataFrame
        TATdict = DFgaps_tat[[3,10]].set_index(3)[10].to_dict()
    except:
        print("WARN: TAT at 0 !!!!!")
        for ii in DFgaps[7].unique():
            TATdict = {ii : 0}
            
    return TATdict
    

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


    
### Internal fct for time conversion
def ixblue_cood_conv(x,index=9):
    """
    Convert a coordinate geographic component from a iXBlue DataFrame row
    ----------
    x : TYPE
        A DataFrame row.
    index : int, optional
        the index of the coordinate component. The default is 9.

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
    
    m = float(entpart[-2:] + "." + decpart)
    d  = int(entpart[-4:-2])
        
    return koef * conv.dms2dec_num(d,m)



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
        Coordinates interpolator.

    """

    ### We switch to Head Pitch Roll (zyx) ########
    Euler = np.column_stack((Head,Pitch,Roll))
    
    # We sort to avoid any conflict
    Time,Euler = utils.sort_binom_list(Time,Euler)
    
    # We generate the interpolator
    Rot_for_slerp = Rotation.from_euler(euler_order,Euler,degrees=degrees)
    Iattitude = conv.Slerp_time(Time,Rot_for_slerp)
    return Iattitude


def iXblueDF_prepro_step1_basic(DFin,TATdict,add_attitude=True):
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
    add_attitude : TYPE, optional
        add attitude in the output DataFrame. The default is True.

    Returns
    -------
    DFout
        a userfriendly DataFrame.
    """

    DFout = pd.DataFrame()
    
    DFout["date_rec"] = DFin.apply(ixblue_time_conv,axis=1)
    ### Extraction of the 4 TWTT (because 4 acoustic heads on GAPS)
    # null TWTT are replaced with NaN ...
    DF_TWTT_nan = DFin[[17,21,25,29]].replace(0,np.nan)
    # ... and TWTT averaged with nanmean
    DFout["TWTTraw"] = DF_TWTT_nan.apply(np.nanmean,axis=1)
    # new column for ID, friendly readable
    DFout["ID_BEA"] = DFin[7]
    # get the emission epoch
    DFout["date_emi"] = DFout["date_rec"] - 1. * pd.to_timedelta(DFout['TWTTraw'],'us')
    
    if add_attitude:
        # get attitude
        DFout["head_rec"] =   DFin[11].values
        DFout["roll_rec"] =   DFin[12].values
        DFout["pitc_rec"] = - DFin[13].values
    
    def internal_twtt_corr(x):
        return x['TWTTraw'] - 1000. * TATdict[x['ID_BEA']]
    
    DFout['TWTT'] = DFout.apply(internal_twtt_corr,axis=1)
    
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
                                   emirec_tup=('emi','rec'),
                                   coord_sys = "NED"):
    """
    Interpolate attitude for the emission and reception epoch
    And apply the lever arms

    Parameters
    ----------
    DFin : pandas DataFrame
        DataFrame generated at step 3.
    InterpolatorAttitude : Scipy Interpolator
        Attitude interpolatorgenerated with interpolator_attitude.
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

        DFout[col_name_DEV[0]],DFout[col_name_DEV[1]],DFout[col_name_DEV[2]] = ABC_DEV.T
        
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
        
        DFout_beacon["TWTTraw"] * 10**-6,
        pd.Series([TATdict[id_beacon]* 10**-3] * len(DFout_beacon["date_emi"]),name='TAT'),
        DFout_beacon["TWTT"] * 10**-6),
        axis=1)
            
        ofile_path = os.path.join(export_dir,
                                  experience_name+".PXP"+str(id_beacon)+".O.dat")
    
        header = 'pxp_coords : ' + str(apriori_coordinates[id_beacon])
    
        np.savetxt(ofile_path,DFout_Ofile,header=header)

    return ofile_path





# def plot_3d(df,df_col=['time','x','y','z']):
#     """
#     a = np.random.rand(2000, 3)*10
#     t = np.array([np.ones(100)*i for i in range(20)]).flatten()
#     df = pd.DataFrame({"time": t ,"x" : a[:,0], "y" : a[:,1], "z" : a[:,2]})
#     """
#     def update_graph(num):
#         print(num)
#         data=df.iloc[num]
#         graph._offsets3d = (data[df_col[1]], 
#                             data[df_col[2]], 
#                             data[df_col[3]])
#         title.set_text('3D Test, time={}'.format(data[df_col[0]]))
    
    
#     fig = plt.figure()
#     ax = fig.add_subplot(111, projection='3d')
#     title = ax.set_title('3D Test')
    
#     data=df[df[df_col[0]]==0]
#     graph = ax.scatter(data[df_col[1]], 
#                        data[df_col[2]], 
#                        data[df_col[3]])
    
#     ani = matplotlib.animation.FuncAnimation(fig, update_graph, 1000, 
#                                    interval=40, blit=False)
    
#     plt.show()
    
    
