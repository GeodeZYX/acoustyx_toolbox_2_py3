# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 18:56:55 2016

ce script a pour but de normaliser les plots entre l'atal et le meteor
version simplifiée pour les maps

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
import SSP as ssp
import vincenty
import mpl_toolkits
from mpl_toolkits.basemap import Basemap , shiftgrid, cm
import matplotlib
# matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
import  os
import genefun as gf
import itertools
import matplotlib.ticker as mtick

#########################################################   
#plt.ioff()

save_plot = 0

PLOT_OUTPUT = "/home/psakicki/THESE/RENDU/1608_Plots_Gradients/PLOT_NORMALIZED"

pickle_path_lis = [
    '/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/scripts/du_travail_sur_les_gradients/ctd_list_pickles/out_ctd_list_ATALwithonemore.pik',

    '/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/scripts/du_travail_sur_les_gradients/ctd_list_pickles/out_ctd_list_ATALjustgwada.pik',    
    '/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/scripts/du_travail_sur_les_gradients/ctd_list_pickles/out_ctd_list_METEORjustgwada.pik',

    '/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/scripts/du_travail_sur_les_gradients/ctd_list_pickles/out_ctd_list_ATAL.pik',
    '/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/scripts/du_travail_sur_les_gradients/ctd_list_pickles/out_ctd_list_METEOR.pik']


pickle_path_lis = ['/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/scripts/du_travail_sur_les_gradients/ctd_list_pickles/out_ctd_list_ATALjustgwada.pik']
pickle_path_lis = ['/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/scripts/du_travail_sur_les_gradients/ctd_list_pickles/out_ctd_list_METEORjustgwada.pik']
pickle_path_lis = ['/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/scripts/du_travail_sur_les_gradients/ctd_list_pickles/out_ctd_list_METEOR.pik']        
    
val_compar = []

for CTSboukle , pickle_path in itertools.product(['C','T','S'],pickle_path_lis):

    ship = pickle_path.split('.')[0].split('_')[-1]
    
    out_ctd_list_work_lis = gf.pickle_loader(pickle_path)
    
    #%%
    
    out_ctd_list_work = out_ctd_list_work_lis[0]
       
    if CTSboukle == 'C':
        upper  = (16.6,-60.5)
        lower  = (15,-62)
        centre = (16,-61.2)
        
        plt.figure()
        mp = Basemap(urcrnrlat=upper[0], urcrnrlon=upper[1] ,
                     llcrnrlat=lower[0], llcrnrlon=lower[1] ,
                     resolution='f',projection='tmerc',
                     lat_0=centre[0],lon_0=centre[1])
        
        parallels = np.arange(0.,90,1)
        mp.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
        # draw meridians
        meridians = np.arange(180.,360.,1.)
        mp.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
        mp.drawcoastlines(color = '0.15')
        mp.fillcontinents(color='tan',lake_color='aqua') #, zorder = 0)
        mp.drawmapboundary(fill_color='aqua')
        
        
        out_ctd_map = out_ctd_list_work
            
        # PLOT ET SELECTION PRIMAIRE
        for iid , ctd in enumerate(out_ctd_map):
            X, Y = mp(ctd.lon,ctd.lat)
            mp.scatter(X,Y,color='k')
            #plt.text(X + 10000 , Y , str(ctd.id),size=15)
        
        plt.savefig(os.path.join(PLOT_OUTPUT,'map_' + ship + '.svg' ))
        plt.savefig(os.path.join(PLOT_OUTPUT,'map_' + ship + '.pdf' ))
        plt.savefig(os.path.join(PLOT_OUTPUT,'map_' + ship + '.png' ))