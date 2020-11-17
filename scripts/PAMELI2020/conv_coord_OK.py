#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  4 11:21:14 2020

@author: psakicki
"""

#### Import star style
from geodezyx import *                   # Import the GeodeZYX modules
from geodezyx.externlib import *         # Import the external modules
from geodezyx.megalib.megalib import *   # Import the legacy modules names

p="/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/04_Doc_Reports/ixBlue_coords.xlsx"
DF = pd.read_excel(p)

refflhdic = dict()

refflhdic[1]=(48.319824833333335, -4.446147         , -40.07)
refflhdic[2]=(48.319821833333336, -4.446828166666666, -40.66)
refflhdic[3]=(48.32009133333333 , -4.446444666666666, -39.68)


DF['latdec'] = DF.lat.apply(conv.dms2dec,onlyDM=True)
DF['londec'] = DF.lon.apply(conv.dms2dec,onlyDM=True)
DF['h'] = -1 * DF['h']

DF[['X','Y','Z']] = np.array(conv.GEO2XYZ_vector(DF[['latdec','londec','h']]))

xyz_ref = np.array([4236463.24737245, -329434.98845412, 4740641.89446747])

DF[['E','N','U']]  = conv.XYZ2ENU_vector(DF[['X','Y','Z']],xyz_ref)

DF['U'] = DF['h']

utils.pickle_saver(DF,"/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING","ixblue_coords")

