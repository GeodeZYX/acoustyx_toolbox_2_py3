# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 11:11:38 2015

@author: psakicki

OBSOLETE => analyz_SSP_MOVE_NOAA/
"""

import acouclass as acls
from netCDF4 import Dataset
import numpy as np
import datetime
import SSP
import datetime as dt

path='/home/psakicki/THESE/bd_CTD/NEW_DATA_2015/ocldb1425484660.15199.OSD.csv'
path='/home/psakicki/THESE/bd_CTD/NEW_DATA_2015/ocldb1425484660.15199.CTD.csv'

export_dir = '/home/psakicki/THESE/bd_CTD/NEW_DATA_2015/NOAA_exported_ssp'
export_data_2_txt = True 

f=open(path)

def str2float_smart(str_in):
    try:
        out = float(str_in)
    except ValueError:
        out = np.nan
    return out


def str2int_smart(str_in):
    try:
        out = int(str_in)
    except ValueError:
        out = np.nan
    return out
    
def get_season(now):
    from datetime import date, datetime
    
    seasons = [('winter', (date(1,  1,  1),  date(1,  3, 20))),
           ('spring', (date(1,  3, 21),  date(1,  6, 20))),
           ('summer', (date(1,  6, 21),  date(1,  9, 22))),
           ('autumn', (date(1,  9, 23),  date(1, 12, 20))),
           ('winter', (date(1, 12, 21),  date(1, 12, 31)))]
    
    if isinstance(now, datetime):
        now = now.date()
    now = now.replace(year=1)
    for season, (start, end) in seasons:
        if start <= now <= end:
            return season
    assert 0, 'never happens'
    
def color_of_season(datein):
    season = get_season(datein)
    if season == 'winter':
        outcolor = 'b'
    elif season == 'summer':
        outcolor = 'r'
    elif season == 'spring':
        outcolor = 'g'
    elif season == 'autumn':
        outcolor = 'k'
    return outcolor

out_ctd_list = []

for l in f:
    f = l.split(',')
    # New CTD
    if l[0] == '#':
        Zlis, Tlis , Slis , Plis = [] , [] , [] , []
        ctd_curr = acls.CTD(Zlis, Tlis , Slis ,Plis)
        out_ctd_list.append(ctd_curr)
        datazone = False
        y,m,d,td = np.nan,np.nan,np.nan,np.nan

    if 'END OF VARIABLES SECTION' in f[0]:
        datazone = False
        
    if 'Latitude' in f[0]:
        ctd_curr.lat = str2float_smart(f[2])
    if 'Longitude' in f[0]:
        ctd_curr.lon = str2float_smart(f[2])
    if 'Year' in f[0]:
        y = str2int_smart(f[2])
    if 'Month' in f[0]:
        m = str2int_smart(f[2])
    if 'Day' in f[0]:
        d = str2int_smart(f[2])
    if 'Time' in f[0]:
        td = datetime.timedelta(str2float_smart(f[2]))
    if ('VARIABLES' in f[0]) and not ('END OF VARIABLES SECTION' in f[0]):
        i_P = f.index('Pressure  ')   
        i_Z = f.index('Depth     ')
        i_T = f.index('Temperatur')
        i_S = f.index('Salinity  ')
    if datazone: 
        Zlis.append(str2float_smart(f[i_Z]))
        Tlis.append(str2float_smart(f[i_T]))
        Slis.append(str2float_smart(f[i_S]))
        Plis.append(str2float_smart(f[i_P]))

    if 'Prof-Flag' in f[0]:
        # the header is ended,
        # the date can be set
        if td is datetime.timedelta:
            ctd_curr.t = datetime.datetime(y,m,d) + td
        else:
            ctd_curr.t = datetime.datetime(y,m,d)      
        
        datazone = True

out_ctd_list[-10].t

out_ssp_lis = []
for ctd in out_ctd_list:
    out_ssp_lis.append(acls.CTD_2_SSP(ctd))
    
plt.clf()   

for ssp in out_ssp_lis:
#    if np.max(ssp.Z) > 2000 and ssp.t > datetime.datetime(2000,1,1,0,0,0):
#        ssp.plot()
    if ssp.t > datetime.datetime(2000,1,1,0,0,0): # and get_season(ssp.t) in 'autumn':
        ssp.plot(color_of_season(ssp.t)+'+',alpha=0.25)
        if export_data_2_txt:
            acls.export_ssp_file(ssp,export_dir,'SSP_NOAA')
            
f = plt.gcf()

#f.savefig('/home/psakicki/ssp.pdf')
f.show()
    
    

    

