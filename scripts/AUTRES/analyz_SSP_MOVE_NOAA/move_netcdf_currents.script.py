# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 15:57:32 2016

@author: psakicki
"""


from netCDF4 import Dataset
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import acouclass 
import glob
import os
import geodetik as geok
import genefun
import SSP as ssp
import pandas as pd
import MOVEclass as mcls
import time
import scipy


path = '/home/psakicki/THESE/bd_CTD/NEW_DATA_2015/MOVE/MOVE3/OTHERS/*CURRENT*'

Lfiles = glob.glob(path)

D = [Dataset(e) for e in Lfiles]

Tgoodstk   = []
Cgoodstk   = []
Dirgoodstk = []

refdatelis = [ dt.datetime(2010, 2, 25, 0) ,
dt.datetime(2009, 10, 2, 0, 0)  ,
dt.datetime(2007, 4, 24) ,
dt.datetime(2012, 12, 4) ]


for refdate in refdatelis:
    print(refdate)

    for d in D:
        T = d['TIME'][:]
        C = np.array(d['CSPD'][:])
        Dir = np.array(d['CDIR'][:])
        
        T = np.array(geok.jjulCNES2dt(T))
        
        booltru = (geok.dt2posix(refdate - dt.timedelta(days=1)) <= np.array(geok.dt2posix(T))) * (np.array(geok.dt2posix(T)) <= geok.dt2posix(refdate + dt.timedelta(days=1.1)))
        if np.sum(booltru) > 0:
            print(refdate)
            Tgood   = T[booltru]
            Cgood   = C[booltru]
            Dirgood = Dir[booltru]
            
            Tgoodstk.append(Tgood)
            Cgoodstk.append(Cgood)
            
            Dirgoodstk.append(Dirgood)


Cgoodstk[-2] = np.delete(Cgoodstk[-2],1,1)  / (100)

Cmean = [np.mean(e,1) for e in Cgoodstk]
Cstd  = [ np.std(e,1) for e in Cgoodstk]

INTERPOS = [scipy.interpolate.interp1d(geok.dt2posix(TT),CC) for TT , CC in zip(Tgoodstk,Cmean)]

col4interpo = []
for TT , CC in zip(Tgoodstk,Cmean):
    col4interpo.append(np.column_stack((np.array(geok.dt2posix(TT)) , CC)))

genefun.pickle_saver(col4interpo,'/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/scripts/analyz_SSP_MOVE_NOAA','CURRENTS_INTERPOLATORS')



#print Cmean
