# -*- coding: utf-8 -*-
"""
Created on Fri Aug  5 11:29:15 2016

@author: psakicki
"""

#def read_ctd_CCHDO_csv(pathin):

from megalib import *
import vincenty
if 1:
    from mpl_toolkits.basemap import Basemap , shiftgrid, cm
import matplotlib
matplotlib.rcParams.update({'font.size': 12})


export     = 0
export_ctd = 0
export_plots = 1
plot_map = 1

ship = 'atal'
ship = 'meteor'
ship = 'japon'


just_gwada = 0
decimate_meteor = False

export_pickle = "/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/scripts/du_travail_sur_les_gradients/ctd_list_pickles"

if ship == 'meteor':
    pathin = '/home/psakicki/THESE/DATA/1608_CTD_METEOR/06MT200508_00398_00001_ct1.csv'
    pathin = '/home/psakicki/THESE/DATA/1608_CTD_METEOR/06MT200508_00398_00001_ct1.csv'
    
    filis = sorted(glob.glob('/home/psakicki/THESE/DATA/1608_CTD_METEOR/*csv'))
    filis = sorted(glob.glob('/home/adminuser/aaa_FOURBI/1608_CTD_METEOR/*csv'))
    
    ctd_lis_raw = [acls.read_CTD_CCHDO_csv(f) for f in filis]
    if decimate_meteor:
        [ctd.decimate() for ctd in ctd_lis_raw]
    out_ctd_list = []
    for iid , ctd in enumerate(ctd_lis_raw):
        ctd.id = iid
        if not (dt.datetime(2005,8,29,3) <= ctd.t <= dt.datetime(2005,9,1,0)):
            continue
        out_ctd_list.append(ctd)
    # RECUPERATION DU TOP
    out_ctd_list = out_ctd_list[:12]
    
    if just_gwada:
        out_ctd_list = out_ctd_list[:4]
        just_gwada_suffix = 'justgwada'
    else:
        just_gwada_suffix = ''
    
elif ship == 'atal':
    path='/home/psakicki/THESE/bd_CTD/ocldb1467049521.874ATALANTE.CTD.csv'
    out_ctd_list_1 = acls.read_NOAA_file_2_CTD_lis(path)

    out_ctd_list_2 = []
    for iid , ctd in enumerate(out_ctd_list_1):
        if dt.datetime(2003,4,1) <= ctd.t <= dt.datetime(2003,4,30):
            ctd.id = iid
            out_ctd_list_2.append(ctd)
        
    #### horiz
    out_ctd_list_horiz = []
    for ctd in out_ctd_list_2:
        if ctd.lat > 16.15 :
            out_ctd_list_horiz.append(ctd)
    
    #### vert
    out_ctd_list_vert = []
    for ctd in out_ctd_list_2:
        print(ctd.lon) 
        if (-60.72 < ctd.lon < -60.62) and (ctd.lat < 16.5 ): #and ( ctd.t > dt.datetime(2003,4,22)):
            out_ctd_list_vert.append(ctd)
    # RETRAIT MANUEL EXTRMEMENT DIRTY D'AUTANT QUE LE RETRAIT DE CE DERNIER
    # MODIFIER LE COMPORTEMENT LA MUNKISATION PREALABLEMENT EFFECTUEE
    # A REFAIRE INCH ALLAH
    # le 1er, au sud, proche d'un autre, pourri => del
    out_ctd_list_vert.remove(out_ctd_list_vert[0])
    # le second, mais premier du coup, au nord, proche d'un autre, pourri => del
    out_ctd_list_vert.remove(out_ctd_list_vert[0])
    # EVENTUELLEMENT, retrait du troisieme qui a été fait avec un dt de 2 jours :
    # out_ctd_list_vert.remove(out_ctd_list_vert[0])

    
    if just_gwada:
        out_ctd_list_vert = out_ctd_list_vert[:4]
        just_gwada_suffix = 'justgwada'
    else:
        just_gwada_suffix = ''
    
    
    #### les petits horiz
    out_ctd_list_horiz_small = []
    for ctd in out_ctd_list_2:
        print(ctd.lon) 
        if (-60.64 < ctd.lon < -60.51) and ( ctd.lat > 16.15 ):
            out_ctd_list_horiz_small.append(ctd)
            
    out_ctd_list = out_ctd_list_vert


elif ship == 'japon':
    path='/home/adminuser/Téléchargements/ocldb1513460594.11356.OSD.csv'
    out_ctd_list_1 = acls.read_NOAA_file_2_CTD_lis(path)
    out_ctd_list = out_ctd_list_1


#
#>>> out_ctd_map[1].lon
#136.85
#>>> out_ctd_map[1].lat
#33.75
#>>> out_ctd_map[0].lon
#136.7167
#>>> out_ctd_map[0].lat
#33.5

[ctd.ssp_calc() for ctd in out_ctd_list]




# Fin chargement

# CARTE
if plot_map:
    upper  = (17,-58)
    lower  = (13,-62)
    centre = (16,-61)
    
    #Japon
    upper  = (40,141)
    lower  = (28,130)
    centre = (34,135)
        
    figmap ,axmap = plt.subplots()
    mp = Basemap(urcrnrlat=upper[0], urcrnrlon=upper[1] ,llcrnrlat=lower[0], llcrnrlon=lower[1],
                 resolution='l',projection='tmerc',lat_0=centre[0],lon_0=centre[1])
    
    mp.fillcontinents(color='#cc9955', lake_color='aqua', zorder = 0)
    mp.fillcontinents(color='tan',lake_color='aqua', zorder = 0)
    mp.drawcoastlines(color = '0.15')
    mp.drawmapboundary(fill_color='aqua')
    print('!!! selection of some CTD in the whole list !!!')
    out_ctd_map = out_ctd_list[-19:-17]

    #out_ctd_map = out_ctd_list_vert
    
    # PLOT ET SELECTION PRIMAIRE
    for iid , ctd in enumerate(out_ctd_map):
        X, Y = mp(ctd.lon,ctd.lat)
        mp.scatter(X,Y,color='k')

if export_ctd:
    piksav = gf.pickle_saver([out_ctd_list],export_pickle,'out_ctd_list_' + ship.upper() + just_gwada_suffix)
    print(piksav)
    
        
#plt.figure()
#latlis = [ctd.lat for ctd in out_ctd_list]
#sm = plt.cm.ScalarMappable(cmap='hot_r', norm=plt.Normalize(vmin=np.min(latlis), vmax=np.max(latlis)))
#cbar = sm.set_array(latlis)
#cbarax = plt.gcf().colorbar(sm,cax=cbar)

#%%

print("munkization")
for iii , ctd in enumerate(out_ctd_list):
    print(iii)
    Ztmp , Ctmp = ssp.munk_pseudo_coef_calc(ctd.Zssp , ctd.C)
    #limdepth1 = 1500
    #limdepth2 = 4500
#
#    limdepth1 = 1500
#    limdepth2 = 4500

    limdepth1 = 1500
    limdepth2 = 1600

    Ctmp2 = Ctmp[Ztmp < limdepth1]
    Ztmp2 = Ztmp[Ztmp < limdepth1]
    #Zm , Cm = ssp.munk(6000)
    ZC = np.loadtxt("/home/adminuser/Documents/CODES/acoustyx_toolbox_2_py3/exemple/input_data_for_simulation/SSP/SSP_NOAA_dep5645_20030608000000")
    Zm , Cm = ZC[:,0] , ZC[:,1]
    Cm = Cm[ limdepth2 <= Zm ]
    Zm = Zm[ limdepth2 <= Zm ]
    
    Ctmp3 = np.hstack((Ctmp2,Cm))
    Ztmp3 = np.hstack((Ztmp2,Zm))
    
    ctd.anex['Ctmp3'] = Ctmp3
    ctd.anex['Ztmp3'] = Ztmp3

    Ztmp4 , Ctmp4 = ssp.munk_pseudo_coef_calc(Ztmp3 , Ctmp3)
    ctd.anex['Zssporig'] = np.array(ctd.Zssp)
    ctd.anex['Corig'] = np.array(ctd.C)
    
    
    #ctd.Zssp , ctd.C = Ztmp , Ctmp
    ctd.Zssp , ctd.C = Ztmp4 , Ctmp4
print("fin munkization")

plt.figure()
if ship == 'meteor':
    ctd_work_lis  = [ctd for ctd in out_ctd_list if (np.max(ctd.Z) > 1500) and (15 < ctd.lat < 16)]   
    ctd_work_lis  = ctd_work_lis[0:2]
else:
    ctd_work_lis  = [ctd for ctd in out_ctd_list if (np.max(ctd.Z) > 1500)]   
    ctd_work_lis  = ctd_work_lis[0:2]

if 1:
    print('reboot for review raytracing article')
    A = ctd_work_lis[0]
    B = ctd_work_lis[1]

    
    
    
    
    
    
Zmax_work_lis = [np.max(ctd.Z)    for ctd in ctd_work_lis] #np.max(ctd.Z) > 2100  and
Zmax_work_lis = [np.max(ctd.Zssp) for ctd in ctd_work_lis] #np.max(ctd.Z) > 2100  and
latlis_work      = np.array([ctd.lat for ctd in ctd_work_lis])

plt.xlim(1485,1545)
plt.ylim(-2500,100)


dist_from_latlis = latlis_work * 1000 * 111.11111
dist_from_latlis = dist_from_latlis
ctd_work_lis     = ctd_work_lis

if plot_map:
    for iid , ctd in enumerate(ctd_work_lis):
        X, Y = mp(ctd.lon,ctd.lat)
        axmap.scatter(X,Y,color='r')

#ctd_work_lis     = ctd_work_lis[1:]
#dist_from_latlis = dist_from_latlis[1:]

dist_from_latlis  = np.flipud(dist_from_latlis - dist_from_latlis[-1])
      
color_dark = ('darkblue' ,'darkred')
color_light = ('dodgerblue' ,'tomato')

[ssp.SSP_plot(ctd.anex['Zssporig'], ctd.anex['Corig'],color=c) for ctd , c in zip(ctd_work_lis,color_dark)]
[ssp.SSP_plot(ctd.Zssp, ctd.C,color=c) for ctd , c in zip(ctd_work_lis,color_light)]
#[ssp.SSP_plot(ctd.anex['Ztmp3'], ctd.anex['Ctmp3'],color='-') for ctd in ctd_work_lis]

#plt.gcf().set_size_inches(gf.Aformat(6))
plt.tight_layout()

if export_plots:
    plt.savefig('/home/psakicki/THESE/RENDU/1608_Plots_Gradients/PLOT_MUNKIZED/' + "SSPs_" + ship + '.pdf')


pseudo_grad = np.diff([ctd.C for ctd in ctd_work_lis],axis=0) / np.diff(dist_from_latlis)[0]
plt.figure()
plt.plot(np.squeeze(pseudo_grad),-ctd_work_lis[0].Zssp)
plt.ylabel('depth (m)')
plt.xlabel('gradient (m/s/m)')

#plt.gcf().set_size_inches(gf.Aformat(6))

import matplotlib.ticker as mtick

#plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0),useOffset=False)
ax = plt.gca()
ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.0e'))
plt.tight_layout()

if export_plots:
    plt.savefig('/home/psakicki/THESE/RENDU/1608_Plots_Gradients/PLOT_MUNKIZED/' + "GRAD_" + ship + '.pdf')


plt.figure()
Zmun,Cmun = ssp.munk(6000)
ssp.SSP_plot(Zmun, Cmun ,color=c)
plt.suptitle('Munk profile')
plt.tight_layout()
plt.subplots_adjust(top=0.90)

if export_plots:
    plt.savefig('/home/psakicki/THESE/RENDU/1608_Plots_Gradients/PLOT_MUNKIZED/' + "MUNK" + '.pdf')


#%%

Zref = ctd_work_lis[0].Zssp

C_grid_stk = [ctd.C for ctd in ctd_work_lis]
C_grid     = np.column_stack(C_grid_stk)


print("begin ssf")
ssf = acls.make_SSF3D_from_SSP_n_distance((dist_from_latlis,Zref,C_grid), 
                                           xmin=-4000, xmax=4000, 
                                           ymin=-4000, ymax=4000)
ssf.Cgrad
print("end   ssf")
 
# CE BLOC EST USELESS, PASSEZ VOTRE CHEMIN !                                   
Zmunk  , Cmunk  = ssp.munk(6000)
Cmunk2 = ssp.munk_pseudo(Zmunk,zref1=1550)
Xtup = (np.array([0.,25000.]),Zmunk,np.column_stack((Cmunk,Cmunk2)))
ssf_munk = acls.make_SSF3D_from_SSP_n_distance(Xtup)
(Cmunk[1] - Cmunk2[1] )/ 25000
#############################################

ssf_munk.Cgrad
AAA = ssf.Cgrad[0]

Udic = dict()
Udic['X']       = dist_from_latlis
Udic['Z_X']     = Zref
Udic['CgridXZ'] = C_grid

if export:
    gf.create_dir('/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/exemple/Udic/')
    picksav = gf.pickle_saver(Udic,'/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/exemple/Udic/','UdicFINAL_' + ship)
    print(picksav)

print(ssf.C)
