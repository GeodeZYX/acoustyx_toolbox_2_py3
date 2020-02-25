# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 12:55:55 2015

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

#filname = 'OS_MOVE2_01_D_MICROCAT-PART01.nc'
#filname = 'OS_MOVE2_01_D_MICROCAT-PART11.nc'
#filname = 'OS_MOVE2_02_D_MINILOGGER-PART05.nc'
#filname = 'OS_MOVE2_05_D_CURRENTMETER.nc'

#filname = 'OS_MOVE3_06_D_CURRENTMETER.nc'
#filname = 'OS_MOVE3_06_D_MINILOGGER-PART01.nc'

# MICRO CAT CONTIENT DES DONNÉES CTD => FABRIQUER UN SSP A PARTIR DE çA
#filname = 'OS_MOVE3_08_D_MICROCAT.nc'

fillist = sorted(glob.glob('/home/psakicki/THESE/bd_CTD/NEW_DATA_2015/MOVE/MOVE3/MICROCAT/*OS_MOVE3_0*_D_MICROCAT-PART*'))
fillist = sorted(glob.glob('/home/psakicki/THESE/bd_CTD/NEW_DATA_2015/MOVE/MOVE3/MICROCAT/*MICROCAT*'))

# TIME WINDOW
start = dt.datetime(2007,1,1)
end   = dt.datetime(2099,1,1)

generik_depth_list = np.array([   38.,    87.,   137.,   234.,   384.,   585.,
                                 835.,  1092.,  1343.,  1594.,  1831.,  2132.,
                                2433.,  2734.,  3043.,  3364.,  3685.,  4004., 
                                4325.,  4645.,  4904.])

M_lis = []
sens_grp = mcls.SensorGroup()

idsens4manu = -1

LONGITUDE = []
LATITUDE  = []
TIME_STK = []
# PART 1 : NETCDF => OBJECTS

for fil in fillist:
        
    print('reading : ' , os.path.basename(fil))
    
    bool_part = False
    
    fil_name   = os.path.basename(fil)
    id_campagn = int(fil_name.split('_')[2])
    
    if 'PART' in fil:
        is_part = True

    D = Dataset(fil)
    
    TIME = [dt.datetime(1950,1,1) + dt.timedelta(days=t) for t in D.variables['TIME'][:]]
    TIME_STK = TIME_STK + TIME
    LONGITUDE.append(D.variables['LONGITUDE'][:][0])
    LATITUDE.append(D.variables['LATITUDE'][:][0])
    
    
    print(len(TIME) , 'mes. inside')
    
    Nsensor = np.array(D.variables['DEPTH'][:]).shape[0]
    Ndim    = len(np.array(D.variables['TEMP'][:]).shape)
    
    print(Nsensor, 'sensors')
    
    
    if Nsensor == 1:
        is_only_1sens = True
        is_manu_search_required = True
    elif Nsensor == 21:
        is_only_1sens = False
        is_manu_search_required = False
    else:
        is_only_1sens = False
        is_manu_search_required = True

    if Nsensor == 1 and Ndim == 2:
        is_need_sqeeze = True
    else:
        is_need_sqeeze = False
   

    # LOADING IN THE ARRAYS
    DEPTH    = np.array(D.variables['DEPTH'])
    if is_need_sqeeze:
        SALINITY = np.squeeze(np.array(D.variables['PSAL'][:]))
        PRESURE  = np.squeeze(np.array(D.variables['PRES'][:]))
        TEMP     = np.squeeze(np.array(D.variables['TEMP'][:]))
    else:
        SALINITY = np.array(D.variables['PSAL'][:])
        PRESURE  = np.array(D.variables['PRES'][:])
        TEMP     = np.array(D.variables['TEMP'][:])        


    # ARRAYS 2 FILES
    # case severals sensors in the array
    if is_only_1sens:
        for ti,t in enumerate(TIME):
            if ti % 20000 == 0:
                print(ti)
            M = mcls.Mesure()
            M.depth = DEPTH[0]
            M.pres  = PRESURE[ti] 
            M.temp  = TEMP[ti]
            M.sali  = SALINITY[ti] 
            M.time = t
            
            M.campagn = id_campagn
            
            _, idsensproto = genefun.find_nearest(generik_depth_list,DEPTH[0])
            M.sensor =  idsensproto+1
            

            sens_grp(int(M.sensor)).append_data(M)


    # case severals sensors in the array
    else:
        for di,d in enumerate(DEPTH):
            for ti,t in enumerate(TIME):
                if ti % 20000 == 0:
                    print(ti)
                M = mcls.Mesure()
                M.depth = d
                M.pres  = PRESURE[ti,di] 
                M.temp  = TEMP[ti,di]
                M.sali  = SALINITY[ti,di] 
                M.time = t
                
                M.campagn = id_campagn
                if not is_manu_search_required:
                    M.sensor = di+1
                else:
                    _, idsens = genefun.find_nearest(generik_depth_list,d)
                    M.sensor =  idsens+1
                    
                    if idsens4manu != idsens:
                        idsens4manu = idsens
                        print('most probable sens.',M.sensor) 
                        

                sens_grp(int(M.sensor)).append_data(M)

    D.close()


# PART 2 : OBJECTS => FILES

path = '/home/psakicki/THESE/bd_CTD/NEW_DATA_2015/MOVE/MOVE3_clean'
path = "/home/psakicki/THESE/bd_CTD/NEW_DATA_2015/MOVE/MOVE_reform_1507/SIMPLE"

for sens in sens_grp.grp:
    print(sens.id)
    sens.check()
    sens.calc_sndspd()
    
    filobj  = open(os.path.join(path,'id'+str(sens.id)+'.simple.big.dat'),'w')
    filobj2 = open(os.path.join(path,'id'+str(sens.id)+'.simple.ssp.dat'),'w')
    
    for m in sens.Data:
        if not (start <= m.time <= end):
            continue
        final_str = ''
        for a in ['time' , 'sensor' ,  'depth' , 'pres' , 'sali' , 'sndspd' , 'bool_valid']:
            final_str = final_str + str(getattr(m,a)) + ' '
        filobj.write(final_str + '\n')
        if m.bool_valid:
            filobj2.write(str(geok.dt2posix(m.time)) + ' ' + str(m.sndspd) + ' ' + str(m.depth) +   '\n')

#    filobj.close()
    filobj2.close()


lon = np.unique(LONGITUDE)
lat = np.unique(LATITUDE)

lon_mean = np.mean(LONGITUDE)
lat_mean = np.mean(LATITUDE)

lon_std = np.std(LONGITUDE)
lat_std = np.std(LATITUDE)

TIME_STK = np.array(TIME_STK)
print(np.min(TIME_STK) , np.max(TIME_STK))

print(lon_mean , lat_mean , lon_std , lat_std)
    


