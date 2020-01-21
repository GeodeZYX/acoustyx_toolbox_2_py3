# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 11:11:38 2015

@author: psakicki

Discontinué pour la fabrication des pickles CTD
METEOR661 fait le boulot pour les 2
Et en plus, les pickles qui sortent sont buggués, alors ...
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

path='/home/psakicki/THESE/bd_CTD/NEW_DATA_2015/NOAA/ocldb1425484660.15199.OSD.csv'
path='/home/psakicki/THESE/bd_CTD/NEW_DATA_2015/NOAA/ocldb1425484660.15199.CTD.csv'
path='/home/psakicki/THESE/bd_CTD/ocldb1467049521.874ATALANTE.CTD.csv'

export_dir        = '/home/psakicki/THESE/bd_CTD/NEW_DATA_2015/NOAA/NEW_ssp/'
export_pickle_ctd = '/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/scripts/du_travail_sur_les_gradients/ctd_list_pickles'

export_2_pickle   = 1
export_data_2_txt = True 
just_gwada        = True

out_ctd_list = acls.read_NOAA_file_2_CTD_lis(path)
[ e.Temp for e in out_ctd_list ]

for i,ctd in enumerate(out_ctd_list):
    print(i,ctd.t)

out_ctd_list_2 = []

upper  = (17,-58)
lower  = (13,-62)
centre = (16,-61)

plt.figure()
mp = Basemap(urcrnrlat=upper[0], urcrnrlon=upper[1] ,llcrnrlat=lower[0], llcrnrlon=lower[1],
             resolution='l',projection='tmerc',lat_0=centre[0],lon_0=centre[1])

mp.fillcontinents(color='#cc9955', lake_color='aqua', zorder = 0)
mp.fillcontinents(color='tan',lake_color='aqua', zorder = 0)
mp.drawcoastlines(color = '0.15')
mp.drawmapboundary(fill_color='aqua')

#%%

for iid , ctd in enumerate(out_ctd_list):
    if dt.datetime(2003,4,1) <= ctd.t <= dt.datetime(2003,4,30):
        ctd.id = iid
        out_ctd_list_2.append(ctd)
        #plt.plot(ctd.lon,ctd.lat,'ob')
        #plt.xlabel('longitude')
        #plt.ylabel('latitude')
        X, Y = mp(ctd.lon,ctd.lat)
        mp.scatter(X,Y,color='k')
        

# le profil horiz est lat > 16.15
# le profil vert est entre -61.72 & -61.62
    
################### selection des Zones ###################   

#### horiz
out_ctd_list_horiz = []
for ctd in out_ctd_list_2:
    if ctd.lat > 16.15 :
        out_ctd_list_horiz.append(ctd)

#### vert
out_ctd_list_vert = []
for ctd in out_ctd_list_2:
    print(ctd.lon) 
    if (-60.72 < ctd.lon < -60.62) and (ctd.lat < 16.50 ):
        out_ctd_list_vert.append(ctd)
        
# le 1er, au sud, proche d'un autre, pourri => del
out_ctd_list_vert.remove(out_ctd_list_vert[0])
# le second, mais premier du coup, au nord, proche d'un autre, pourri => del
out_ctd_list_vert.remove(out_ctd_list_vert[0])
# EVENTUELLEMENT, retrait du troisieme qui a été fait avec un dt de 2 jours :
out_ctd_list_vert.remove(out_ctd_list_vert[0])


if just_gwada:
    out_ctd_list_vert = out_ctd_list_vert[:3]
    just_gwada_suffix = 'justgwada'
else:
    just_gwada_suffix = ''


#### les petits horiz
out_ctd_list_horiz_small = []
for ctd in out_ctd_list_2:
    print(ctd.lon) 
    if (-60.64 < ctd.lon < -60.51) and ( ctd.lat > 16.15 ):
        out_ctd_list_horiz_small.append(ctd)
        
#### selesction par id
out_ctd_list_id_lis = []

for id_lis in ((23,29) , (21,28) , (20,27)):
    out_ctd_list_id = []
    for ctd in out_ctd_list_2:
        if ctd.id in id_lis:
            out_ctd_list_id.append(ctd)   
    out_ctd_list_id_lis.append(out_ctd_list_id)

#########################################################   

out_ctd_list_work_lis = out_ctd_list_id_lis



if 1:
    out_ctd_list_work_lis = [out_ctd_list_horiz]
    out_ctd_list_work_lis = [out_ctd_list_vert]

figmean , axmean = plt.subplots()
figall   , axall = plt.subplots()

if export_2_pickle:
    outpik = gf.pickle_saver(out_ctd_list_work_lis,export_pickle_ctd,'out_ctd_list_ATALwithonemore' + just_gwada_suffix)

#%%

for out_ctd_list_work in out_ctd_list_work_lis:
    
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
    ctdmaxZ = ctdmaxz.Z
    
    for ctd in out_ctd_list_work:
        ctd.C    = np.array(list(ctd.C) + list(np.ones(len(ctdmaxZ) - len(ctd.Z)) * np.nan))
        ctd.Zssp = np.array(ctdmaxZ)
        ctd.Temp = np.array(list(ctd.Temp) + list(np.ones(len(ctdmaxZ) - len(ctd.Z)) * np.nan))
        ctd.Sali = np.array(list(ctd.Sali) + list(np.ones(len(ctdmaxZ) - len(ctd.Z)) * np.nan))
        
    
    latlon = [(ctd.lat,ctd.lon) for ctd in out_ctd_list_work]
    
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
    #    
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
        
        CTS = 'C'
        if  CTS == 'C':
            GRADlis = GRAD_C_lis
            axall.set_xlabel('Gradient (m/s/m)')

        elif CTS == 'T':
            GRADlis = GRAD_T_lis
            axall.set_xlabel('Temp Gradient  (Celsius/m)')  

        elif CTS == 'S':
            GRADlis = GRAD_S_lis
            axall.set_xlabel('Salinity Gradient  (PSU/m)')  
        
        #GRADlis = GRAD_T_lis

        PROFlis.append(ctd1.Zssp)
        
        label = genefun.join_improved( '' ,'d=',D[i] *.001,'km, b/w ', ctd1.t , ' id' , ctd1.id , ' and ', ctd2.t , ' id' , ctd2.id )
        axall.plot(GRADlis[-1],-PROFlis[-1],label=label)
        axall.set_ylabel('Depth (m)')

    
    plt.legend()
    
    #ylimbound = np.array(genefun.ylim_easy(-PROFlis[-1],0.05))
    #ylimbound[0] = -4000
    #xlimbound = (-7 * 10**-5 , 1.5 * 10**-4 )
    #axall.set_ylim(ylimbound)
    #axall.set_xlim(xlimbound)
    
    GRADarr = np.column_stack(GRADlis)
    Grad_mean = np.nanmean(GRADarr,1)
    
    axmean.plot(Grad_mean,-PROFlis[-1])
    axmean.set_ylabel('Depth (m)')
    axmean.set_xlabel('Speed Gradient  (m/s/m)')  
    axmean.set_xlabel('Temp Gradient  (Celsius/m)')  
    axmean.set_xlabel('Salinity Gradient  (PSU/m)')  
    
    print("=========== FINS DES GRAPHS ===========")

    Zgrad = PROFlis[-1] 
    Grad = Grad_mean
    
    Grad[np.isnan(Grad)] = 0
    Grad[Zgrad >= 4000] = 0
    
    
    polykoefs = np.polyfit(Zgrad,Grad,1)
    Zgradpoly = np.arange(np.min(Zgrad) , np.max(Zgrad),0.01)
    GradPoly  = np.polyval(polykoefs,Zgradpoly)
    #GradPoly[np.isnan(Grad)] = 0
    #GradPoly[Zgradpoly >= 4000] = 0
    Ind_nullify_zone = np.squeeze(np.argwhere((3700 < Zgradpoly) * (Zgradpoly < 4000)))
    foundind = False
    
    GradPoly_of_nullifyzone = GradPoly[Ind_nullify_zone]
    indZe = np.where( np.min(np.abs(GradPoly_of_nullifyzone)) == np.abs(GradPoly) )[0]
    GradPoly_argsort = np.argsort(np.abs(GradPoly))
    
#    print GradPoly_argsort    
#    
#    Ind_nullify_zone_sorted = []
#    
#    for ind in GradPoly_argsort:
#        if ind in Ind_nullify_zone:
#            Ind_nullify_zone_sorted.append(ind)
#            
#    print np.array(Ind_nullify_zone_sorted)
    print(indZe) 
    
    GradPoly2 = np.array(GradPoly)
    GradPoly2[indZe:] = 0
    
    plt.figure()
    plt.plot(Zgrad , Grad          , '*b')
    plt.plot(Zgradpoly , GradPoly  , '+g')
    plt.plot(Zgradpoly , GradPoly2 , 'xk')
    
    fig = plt.gcf()
    fig.set_size_inches(18.5, 10.5)

    ZGrad_arr = np.column_stack((Zgrad, Grad))    
    ZGrad_arr = np.column_stack((Zgradpoly, GradPoly2))    
    
    if 0:
        print("strt save")
        outdirKfile  = '/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/exemple/Kfiles_kourents_gradients'
        outnameKfile = 'kfiletestpoly.K'
        outpathKfile = os.path.join(outdirKfile,outnameKfile)
        np.savetxt(outpathKfile,ZGrad_arr)
        
        plt.savefig(outdirKfile + '/plotgrad.png')
    
#ZTSP         = acls.read_NOAA_file_2_CTD_lis(path,True)

#pour bien montrer que c'est un objet à l'interieur
#out_ctd_list[-10].t

#out_ssp_lis = []
#for ctd in out_ctd_list:
#    out_ssp_lis.append(acls.CTD_2_SSP(ctd))
    
#plt.clf()   

#for ssp in out_ssp_lis:
##    if np.max(ssp.Z) > 2000 and ssp.t > datetime.datetime(2000,1,1,0,0,0):
##        ssp.plot()
#    if ssp.t > datetime.datetime(2000,1,1,0,0,0): # and get_season(ssp.t) in 'autumn':
#        ssp.plot(geok.color_of_season(ssp.t)+'+',alpha=0.25)
#        if export_data_2_txt:
#            acls.export_ssp_file(ssp,export_dir,'NOAA_CTD')
    
    
            
f = plt.gcf()
#f.savefig('/home/psakicki/ssp.pdf')
f.show()
    
    

    

