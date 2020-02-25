# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 18:56:55 2016

ce script a pour but de normaliser les plots entre l'atal et le meteor

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
    
    
val_compar = []

for CTSboukle , pickle_path in itertools.product(['C','T','S'],pickle_path_lis):


    ship = pickle_path.split('.')[0].split('_')[-1]
    
    out_ctd_list_work_lis = gf.pickle_loader(pickle_path)
    
    figmean      , axmean    = plt.subplots()
    figall       , axall     = plt.subplots()
    figcontour   , axcontour = plt.subplots()

    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    axmean.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
    axall.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
    axcontour.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2f'))

    figall.set_size_inches((figall.get_size_inches()*1.35))
    figmean.set_size_inches((figmean.get_size_inches()*1.35))
    figcontour.set_size_inches((figcontour.get_size_inches()*1.35))
    
    #%%
    
    out_ctd_list_work = out_ctd_list_work_lis[0]
       
    if CTSboukle == 'C':
        upper  = (17,-59.6)
        lower  = (13.4,-62)
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
            plt.text(X + 10000 , Y , str(ctd.id),size=15)
        
        plt.savefig(os.path.join(PLOT_OUTPUT,'map_' + ship + '.svg' ))
        plt.savefig(os.path.join(PLOT_OUTPUT,'map_' + ship + '.pdf' ))
        plt.savefig(os.path.join(PLOT_OUTPUT,'map_' + ship + '.png' ))
        
        out_ctd_list_work_orig = list(out_ctd_list_work)
    
    
    giclage_faible_prof = False
    if giclage_faible_prof:
        out_ctd_list_work_new = []
        for ctd in out_ctd_list_work:
            if np.max(ctd.Z) > 1000:
                out_ctd_list_work_new.append(ctd)
        out_ctd_list_work = out_ctd_list_work_new
    
    minz = np.min([np.max(ctd.Z) for ctd in out_ctd_list_work])
    maxz = np.max([np.max(ctd.Z) for ctd in out_ctd_list_work])
    
    
    for ctd in out_ctd_list_work:
        ctd.ssp_calc()
    
    ctdmaxz = [ctd for ctd in out_ctd_list_work if maxz == np.max(ctd.Z)][0]
    ctdmaxZ = np.array(ctdmaxz.Z)
    
    maxz_lis = [-np.max(ctd.Z) for ctd in out_ctd_list_work]

    CTScontour = []    
    
    for ctd in out_ctd_list_work:
        ctd.C    = np.array(list(ctd.C) + list(np.ones(len(ctdmaxZ) - len(ctd.Z)) * np.nan))
        ctd.Temp = np.array(list(ctd.Temp) + list(np.ones(len(ctdmaxZ) - len(ctd.Z)) * np.nan))
        ctd.Sali = np.array(list(ctd.Sali) + list(np.ones(len(ctdmaxZ) - len(ctd.Z)) * np.nan))
        ctd.Zssp = np.array(ctdmaxZ)
        
        if   CTSboukle == 'C':
            CTScontour.append(ctd.C)    
        elif CTSboukle == 'T':
            CTScontour.append(ctd.Temp)   
        elif CTSboukle == 'S':
            CTScontour.append(ctd.Sali)
        
    
    latlon  = [(ctd.lat,ctd.lon) for ctd in out_ctd_list_work]
    lat_lis = [(ctd.lat) for ctd in out_ctd_list_work]
    
    try:
        import matplotlib.cm as cm
                
        if   CTSboukle == 'C':
            normbound = (1480,1550)
            levels = np.arange(1480,1551,5)
        elif CTSboukle == 'T':
            normbound = (3,30.5)
            levels = None
        elif CTSboukle == 'S':
            normbound = (34,38)
            levels = None

#        normbound = (np.nanmin(np.column_stack(CTScontour)),
#                     np.nanmax(np.column_stack(CTScontour)))

        tup_4_val_compar = (os.path.basename(pickle_path) , CTSboukle , normbound)
        print(tup_4_val_compar)
        val_compar.append(tup_4_val_compar)
        
        
        CM     =  matplotlib.colors.Colormap('viridis')
        NormCM =  matplotlib.colors.Normalize(normbound[0],normbound[1])

        Xcontour, Ycontour = np.meshgrid(lat_lis, -ctdmaxZ)
        CS = axcontour.contour(Xcontour, Ycontour, 
                               np.column_stack(CTScontour),cmap=cm.viridis,
                               norm=NormCM,levels = levels)
        plt.clabel(CS, inline=1, fontsize=10)
        #figcontour.set_title('Simplest default with labels')
        axcontour.plot(lat_lis,maxz_lis,'ro')
        axcontour.set_ylabel('Depth (m)')
        axcontour.set_xlabel('Latitude')
        axcontour.set_ylim((-1250,10))
    except:
        pass
    
    
    D = []
    for i in range(len(latlon)-1):
        D.append(vincenty.vincenty(latlon[i],latlon[i+1]) * 1000)
        
    GRADlis = []
    PROFlis = []
    
    GRAD_C_lis = []
    GRAD_T_lis = []
    GRAD_S_lis = []
        
    for i in range(len(D)):
        ctd1 = out_ctd_list_work[i]
        ctd2 = out_ctd_list_work[i+1]
            
    #    boolis1 = np.array(len(ctd1.Z) * [True])
    #    boolis2 = np.array(len(ctd2.Z) * [True])

        C1 = ctd1.C
        C2 = ctd2.C
    
        T1 = ctd1.Temp
        T2 = ctd2.Temp
        
        S1 = ctd1.Sali
        S2 = ctd2.Sali
        
        GRAD_C_lis.append((C2 - C1) / D[i])
        GRAD_T_lis.append((T2 - T1) / D[i])
        GRAD_S_lis.append((S2 - S1) / D[i])
     
    #    print np.all(ctd1.Zssp[boolis1] == ctd2.Zssp[boolis2])
        
        CTS = CTSboukle
        
        if  CTS == 'C':
            GRADlis = GRAD_C_lis
            axmean.set_xlabel('Speed Gradient  (m/s/m)')  
            axall.set_xlabel('Speed Gradient (m/s/m)')
    
        elif CTS == 'T':
            GRADlis = GRAD_T_lis
            axmean.set_xlabel('Temp Gradient  (Celsius/m)')  
            axall.set_xlabel('Temp Gradient  (Celsius/m)')  
    
        elif CTS == 'S':
            GRADlis = GRAD_S_lis
            axmean.set_xlabel('Salinity Gradient  (PSU/m)')  
            axall.set_xlabel('Salinity Gradient  (PSU/m)') 
    
        #GRADlis = GRAD_T_lis
    
        PROFlis.append(ctd1.Zssp)
        deltat = ctd2.t - ctd1.t
        deltat = str(deltat)[:-3]

        label = genefun.join_improved( '' ,'d=',np.round(D[i] *.001),'km, b/w ',
                                      ctd1.id , ' & ' , ctd2.id, ', dt = ' , (deltat))
        axall.plot(GRADlis[-1],-PROFlis[-1],label=label)
        axall.set_ylabel('Depth (m)')
    axall.legend()
    
    #ylimbound = np.array(genefun.ylim_easy(-PROFlis[-1],0.05))
    #ylimbound[0] = -4000
    #xlimbound = (-7 * 10**-5 , 1.5 * 10**-4 )
    #axall.set_ylim(ylimbound)
    #axall.set_xlim(xlimbound)
    
    GRADarr   = np.column_stack(GRADlis)
    Grad_mean = np.nanmean(GRADarr,1)
    Grad_std  = np.nanstd(GRADarr,1)
    
    axmean.plot(Grad_mean,-PROFlis[-1])
    #axmean.errorbar(Grad_mean,-PROFlis[-1],xerr=Grad_std)
    axmean.set_ylabel('Depth (m)')
    
    print("=========== FINS DES GRAPHS ===========")
    
    if save_plot:
        for ext in ('.svg','.pdf','.png'): 
            def namfct(nam):
                return '_'.join((ship,CTSboukle,nam,ext))
        
            figall.savefig( os.path.join(PLOT_OUTPUT,namfct('all')))
            figmean.savefig( os.path.join(PLOT_OUTPUT,namfct('mean')))
            figcontour.savefig(os.path.join(PLOT_OUTPUT,namfct('contour')))
            
    print("=========== FINS DES SAVES ===========")

print(val_compar)

#array([[ 1485.86029467,  1542.22034796],
#       [ 1487.5518    ,  1541.31339484],
#       [ 1487.32876513,  1545.47721572],
#       [ 1485.86029467,  1542.22034796],
#       [ 1486.76225119,  1546.21722782],
#       [    3.1592    ,    27.9468    ],
#       [    3.1592    ,    27.4992    ],
#       [    3.9311    ,    29.8219    ],
#       [    3.1592    ,    27.9468    ],
#       [    3.1756    ,    30.0922    ],
#       [   34.6247    ,    37.2601    ],
#       [   34.7278    ,    36.9604    ],
#       [   34.7362    ,    37.1573    ],
#       [   34.6247    ,    37.2601    ],
#       [   34.5324    ,    37.2086    ]])