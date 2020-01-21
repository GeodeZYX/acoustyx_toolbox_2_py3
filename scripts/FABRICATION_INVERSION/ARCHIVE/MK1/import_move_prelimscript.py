# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 12:55:55 2015

@author: psakicki

VERSION OBSOLETE
Utiliser analyz_SSP_MOVE_NOAA/move_netcdf2simple.script.py


"""

from netCDF4 import Dataset
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import acouclass 
import glob
import os
import geodetik as geok
import SSP as ssp
import pandas as pd
import MOVEclass as mcls

#filname = 'OS_MOVE2_01_D_MICROCAT-PART01.nc'
#filname = 'OS_MOVE2_01_D_MICROCAT-PART11.nc'
#filname = 'OS_MOVE2_02_D_MINILOGGER-PART05.nc'
#filname = 'OS_MOVE2_05_D_CURRENTMETER.nc'

filname = 'OS_MOVE3_06_D_CURRENTMETER.nc'
filname = 'OS_MOVE3_06_D_MINILOGGER-PART01.nc'

# MICRO CAT CONTIENT DES DONNÉES CTD => FABRIQUER UN SSP A PARTIR DE çA
filname = 'OS_MOVE3_08_D_MICROCAT.nc'

fillist = sorted(glob.glob('/home/psakicki/THESE/bd_CTD/NEW_DATA_2015/MOVE/MOVE3/*MICROCAT*'))
BIGTIME = []

M_lis = []
sens_grp = mcls.SensorGroup()
for fil in fillist:
    
    bool_part = False
    
    fil_name   = os.path.basename(fil)
    id_campagn = int(fil_name.split('_')[2])
    
    if 'PART' in fil:
        bool_part = True
        continue
          
    d = Dataset(fil)
    
    TIME = [dt.datetime(1950,1,1) + dt.timedelta(days=t) for t in d.variables['TIME'][:]]
    
    SALINITY = np.squeeze(np.array(d.variables['PSAL'][:]))
    PRESURE  = np.squeeze(np.array(d.variables['PRES'][:]))
    TEMP     = np.squeeze(np.array(d.variables['TEMP'][:]))
    if not bool_part:
        DEPTH = np.squeeze(np.array(d.variables['DEPTH'][:]))
    else:
        DEPTH = np.array(d.variables['DEPTH'][:])
    
    if not bool_part:
        for di,d in enumerate(DEPTH):
            for ti,t in enumerate(TIME):
                M = mcls.Mesure()
                M.depth = d
                M.pres  = PRESURE[ti,di] 
                M.temp  = TEMP[ti,di]
                M.sali  = SALINITY[ti,di] 
                M.time = t
                
                M.campagn = id_campagn
                M.sensor = di+1
                
                sens_grp(M.sensor).append_data(M)
#                sens_grp.append_mes_in_rigth_sens(M)
                  

path = '/home/psakicki/THESE/bd_CTD/NEW_DATA_2015/MOVE/MOVE3_clean'

for sens in sens_grp.grp:
    print(sens.id)
    sens.check()
    sens.calc_sndspd()
    
    filobj = open(os.path.join(path,'data_'+str(sens.id)),'w')
    filobj2 = open(os.path.join(path,'ssp_'+str(sens.id)),'w')
    
    for m in sens.Data:
        final_str = ''
        for a in ['time' , 'sensor' ,  'depth' , 'pres' , 'sali' , 'sndspd' , 'bool_valid']:
            final_str = final_str + str(getattr(m,a)) + ' '
        filobj.write(final_str + '\n')
        if m.bool_valid:
            filobj2.write(str(geok.dt2posix(m.time)) + ' ' + str(m.sndspd) + ' ' + str(m.depth) +   '\n')


    filobj.close()
    filobj2.close()
    

print("aaaa")
