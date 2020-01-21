# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 18:39:17 2015

@author: psakicki
"""

import glob
import numpy as np
import matplotlib.pyplot as plt
import geodetik as geok
from scipy.fftpack import fft
import scipy
import raytrace as rt

import datetime as dt
import MOVEclass as mcls

from scipy.optimize import curve_fit

path= '/home/psakicki/THESE/bd_CTD/NEW_DATA_2015/MOVE/MOVE3_clean/ssp*6'
path_fig = '/home/psakicki/THESE/bd_CTD/NEW_DATA_2015/MOVE/MOVE3_clean/'
fil_lis = glob.glob(path)

sens_grp = mcls.read_multi_ssp_file(path_fig,'decssp*')     
period_lis = mcls.search_multi_period_sensor(sens_grp(16),dt.timedelta(days=10))

a_sensor = sens_grp(16)

strt = dt.datetime(2010,1,3)
end  = dt.datetime(2010,1,13,0,0,0)

win_epoc_lis = a_sensor.get_windowed_epoch_list(strt,end)

f = plt.figure()
    
epoch_mes = sens_grp.get_mesure_epoch(end)

SSP = [e.sndspd for e in epoch_mes]
Z   = [e.depth  for e in epoch_mes]
    
modkriter = 25
for i,e in enumerate(win_epoc_lis):
    if not np.mod(i,modkriter) ==0:
        continue
    print(i,e)
    epoch_mes = sens_grp.get_mesure_epoch(e)

    SSP = np.array([m.sndspd for m in epoch_mes])
    Z   = np.array([m.depth  for m in epoch_mes])
    
    outZC = np.vstack((Z,SSP)).T 
    
    print(outZC.shape)
    
    dirout = '/home/psakicki/THESE/bd_CTD/NEW_DATA_2015/MOVE/MOVE3_epoch'    
    e_str = e.strftime('%Y%m%d%H%M%S')
    filout = 'SSP_epoch_'+e_str
    pathout = os.path.join(dirout,filout)
    np.savetxt(pathout,outZC)    
    
#    plt.plot(SSP,-np.array(Z),'+')
    mcls.plot_list_of_mesures(epoch_mes,f)

#coef_pseudo_munk, _ = curve_fit(mcls.exp_dble_fct,Z,SSP,mcls.apri_for_exp_dble())
#
#Zfit = np.arange(0,5000)
#Cfit = mcls.exp_dble_fct(Zfit,*coef_pseudo_munk)
#
#plt.figure()
#plt.plot(SSP,-np.array(Z))
#plt.plot(Cfit,-Zfit)
