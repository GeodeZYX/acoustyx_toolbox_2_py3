# -*- coding: utf-8 -*-

from megalib import  *

path = "/home/psakicki/THESE/bd_CTD/NEW_DATA_2015/MOVE/MOVE_reform_1507/SIMPLE/"
pathplot  = '/home/psakicki/THESE/bd_CTD/NEW_DATA_2015/MOVE/MOVE_reform_1507/SIMPLE/plots'
pathspace = '/home/psakicki/THESE/bd_CTD/NEW_DATA_2015/MOVE/MOVE_reform_1507/SIMPLE/spatial'
pathspace = '/home/psakicki/THESE/bd_CTD/NEW_DATA_2015/MOVE/MOVE_reform_1512/SIMPLE/spatial2'

wildcard = "**ssp.dat"

# PART 0 : CHARGEMENT DES DONNEES MOVE
if 0:
    SensGrp = mcls.read_multi_ssp_file(path,wildcard)

# PART 1 : PLOTs TEMPORELS SUR UN SEUL SENSEUR
# 1a : fenêtrons individuellement un senseur independant

if 0:
    a_sensor = SensGrp(7)
    strt = dt.datetime(2010,1,3)
    end  = dt.datetime(2011,12,13,0,0,0)
    a_sensor_windowed , _ = mcls.search_period_sensor(a_sensor,strt,end)
    #T,Z,C = a_sensor_windowed.get_TZC()
    #plt.clf()
    #plt.plot(T,C)
    mcls.plot_TC_frontend(a_sensor_windowed,path=path)
    
    
    # 1b : faisons des petits bout de period_size j sur l'intégralité de la periode
    period_size = 10
    tdelta = dt.timedelta(days=period_size)
    windoweds_asensor_lis = mcls.search_multi_period_sensor(a_sensor,tdelta)
    
    for winsens in windoweds_asensor_lis[:]:
        mcls.plot_TC_frontend(winsens,1,path=path)

# PART 2 : PLOTs TEMPORELS SUR TOUS LES SENSEURS (OPERATIONNEL)

reload(mcls)

if 0:
    period_size_lis = [100,10,1]
    for period_size in period_size_lis:
        tdelta = dt.timedelta(days=period_size)
        for sensor in SensGrp.grp:
            print('sens' , sensor.id , 'period' , period_size)
            windoweds_sensor_lis = mcls.search_multi_period_sensor(sensor,tdelta)
            for winsens in windoweds_sensor_lis[:]:
                pathplot2 = os.path.join(pathplot,str(period_size),str(winsens.id))
                gf.create_dir(pathplot2)
                mcls.plot_TC_frontend(winsens,1,path=pathplot2,minmax=[sensor.min_sv,sensor.max_sv])

# PART 3 : SPATIALs SSP

if 0:
    reload(mcls)

    # On definit un senseur comme la reference pour les epochs
    ref_sensor = SensGrp(7)

    # On definit éventuellement une fenetre mais c'est facultatif
    strt = dt.datetime(2010,1,3)
    end  = dt.datetime(2010,1,13,0,0,0)
    
    strt = dt.datetime(1000,1,3)
    end  = dt.datetime(2099,12,13,0,0,0)

    # On trouve toutes les epochs
    ref_epoc_lis = ref_sensor.get_windowed_epoch_list(strt,end)

    # et pour tous les senseurs, on trouve les mesures associées à l'époque
    for i,e in enumerate(ref_epoc_lis):
        epoc_mes_lis = SensGrp.get_mesure_epoch(e)
        Z,C = mcls.ZC_from_mesure_lis(epoc_mes_lis)
        ZC  = np.vstack((Z,C)).T
        epoch_str=e.strftime("%Y%m%d_%H%M")
        gf.create_dir(pathspace)
        filepath = os.path.join(pathspace,epoch_str+'.space.ssp.dat')
        header = "epoch : " + str(e) 
        np.savetxt(filepath,ZC,header=header)
        
# Part 4 : on recharge les SSP en on en fait des plot
if 1:
    ssp_fil_lis = sorted(glob.glob('/home/psakicki/THESE/bd_CTD/NEW_DATA_2015/MOVE/MOVE_reform_1512/SIMPLE/spatial/*ssp.dat'))
    Zstk = []
    Cstk = []
    
    SSPmean = []
    datelis = []  
    dateposixlis = []    
    
    strt = dt.datetime(1000,1,3)
    end  = dt.datetime(2099,12,13,0,0,0)
    
    for sspfil in ssp_fil_lis:
             
        d = geok.date_string_2_dt(gf.read_comments(sspfil)[0][8:])

        if not (strt < d < end):
            continue

        datelis.append(d)
        print(d)
        
        ZC = np.loadtxt(sspfil)
        Z = ZC[:,0]
        C = ZC[:,1]
        ssp.SSP_plot(Z,C)
        SSPmean.append(ssp.SSP_mean(Z,C,0,zmax=np.max(Z)))
        dateposixlis.append(geok.dt2posix(d))
        
        Zstk.append(Z)
        Cstk.append(C)
  
A = np.column_stack((dateposixlis , SSPmean))

meanpath = '/home/psakicki/THESE/bd_CTD/NEW_DATA_2015/MOVE/MOVE_reform_1512/SIMPLE/spatial/MEAN.ssp'
#np.savetxt(meanpath,A)
A = np.loadtxt(meanpath)

plt.clf()
ax = plt.gca()
ax.ticklabel_format(useOffset=False)

import geodetik as geok
A2 , bbA = geok.outiler_mad(A[:,1])

#plt.plot(geok.posix2dt(A[:,0]),A[:,1],'+r')
plt.plot(geok.posix2dt(A[:,0][bbA]),A2,'.')
plt.xlabel('Time')
plt.ylabel('Sound Speed (m/s)')

            
    
    


        
        