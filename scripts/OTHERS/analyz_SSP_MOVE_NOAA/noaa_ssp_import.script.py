# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 11:11:38 2015

@author: psakicki

"""

import acouclass as acls
from netCDF4 import Dataset
import numpy as np
import datetime
import geodetik as geok
import genefun
import SSP
import datetime as dt

path='/home/psakicki/THESE/bd_CTD/NEW_DATA_2015/NOAA/ocldb1425484660.15199.OSD.csv'
path='/home/psakicki/THESE/bd_CTD/NEW_DATA_2015/NOAA/ocldb1425484660.15199.CTD.csv'


export_dir = '/home/psakicki/THESE/bd_CTD/NEW_DATA_2015/NOAA/NEW_ssp/'
export_data_2_txt = True 

out_ctd_list = acls.read_NOAA_file_2_CTD_lis(path)

for i,ctd in enumerate(out_ctd_list):
    print(i,ctd.t)
    
    
# 297 ID de ref
plt.plot(out_ctd_list[297].Temp ,  -np.array(out_ctd_list[297].Z))
plt.xlabel('temperature (deg)')
plt.ylabel('depth (m)')

#ZTSP         = acls.read_NOAA_file_2_CTD_lis(path,True)

#pour bien montrer que c'est un objet à l'interieur
out_ctd_list[-10].t

out_ssp_lis = []
for ctd in out_ctd_list:
    out_ssp_lis.append(acls.CTD_2_SSP(ctd))
    
plt.clf()   

for ssp in out_ssp_lis:
#    if np.max(ssp.Z) > 2000 and ssp.t > datetime.datetime(2000,1,1,0,0,0):
#        ssp.plot()
    if ssp.t > datetime.datetime(2000,1,1,0,0,0): # and get_season(ssp.t) in 'autumn':
        ssp.plot(geok.color_of_season(ssp.t)+'+',alpha=0.25)
        if export_data_2_txt:
            acls.export_ssp_file(ssp,export_dir,'NOAA_CTD')
    
    
            
f = plt.gcf()
#f.savefig('/home/psakicki/ssp.pdf')
f.show()
    
    

    

