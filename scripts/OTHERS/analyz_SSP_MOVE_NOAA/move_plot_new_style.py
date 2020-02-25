# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 15:01:15 2016

@author: psakicki
"""

from   megalib import *
import matplotlib.dates as mdates
import genefun as gf
import matplotlib.pyplot as plt


dicpath  = '/home/psakicki/THESE/bd_CTD/NEW_DATA_2015/MOVE/MOVE_reform_1507/SIMPLE/20160408_145816_sens_grp_minidico.pik'
dicpath  = "/media/psakicki/2FDEABBF02892416/ArchiveThese/Geo1_2To/CALIPSO/psakicki/THESE/bd_CTD/NEW_DATA_2015/MOVE/MOVE_reform_1507/SIMPLE/20160408_145816_sens_grp_minidico.pik"

### Old path Calipso
pathplot = '/home/psakicki/THESE/bd_CTD/NEW_DATA_2015/MOVE/MOVE_reform_1507/SIMPLE/plots_newstyle'
dirplot2 = "/home/psakicki/THESE/RENDU/1605_MANUSCRIT/FIG/MOVE_plots"

### new paths
pathplot = '/home/psakicki/GFZ_WORK/RENDU/1802_Article_LSQ_Simu_GNSSA/03_REVIEW_1/plot_ssp'
dirplot2 = "/home/psakicki/GFZ_WORK/RENDU/1802_Article_LSQ_Simu_GNSSA/03_REVIEW_1/plot_ssp"


if 1:
    dico = gf.pickle_loader(dicpath)
    print('LOADING FINISHED')

min_gal        = np.min([np.min(d[:,-1]) for d in list(dico.values())])
max_gal        = np.max([np.max(d[:,-1]) for d in list(dico.values())])
minmax_general = [min_gal , max_gal]

#plt.ioff()

with_plot = True

period_size_list = reversed([10,100,1000])

# Boucle de plot pour l'ensemble de la periode
#
if 0:
    for period_size in period_size_list:
        period_dt   = dt.timedelta(days=period_size)
        
        for iid , data in dico.items():
            print("sensor " , iid)
            strt = data[0,0]
            end  = strt + period_dt
            
            minmax = [np.min(data[:,-1]) , np.max(data[:,-1])]
        
            while end <= data[-1,0]:
                print('strt , end : ' , strt , end)
                bbool = (data[:,0] >= strt) * (data[:,0] < end)
                dataplot = data[bbool,:]
                T = dataplot[:,0]
                C = dataplot[:,-1]
                
                fig, ax1 = plt.subplots()
                ax2 = ax1.twinx()
                fig.set_size_inches(16.53,11.69)  
        
                y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
                ax1.yaxis.set_major_formatter(y_formatter)
                ax2.yaxis.set_major_formatter(y_formatter)
                
                ax1.plot(T,C,'.b')
                ax2.plot(T,C,'.g') 
                
                stats = ['Sound Speed (m/s) :'       , 
                         'mean : '    , np.mean(C)   ,
                         ', std. : '  , np.std(C)    ,
                         ', median : ', np.median(C) , 
                         ', min : '   , np.min(C)    ,
                         ', max : '   , np.max(C)    ]
                
                stattxt = gf.join_improved(' ' , *stats)
                plt.text( 0.00, -0.18 , stattxt , transform = ax1.transAxes)
        
                if minmax != []:
                    ax2.set_ylim(minmax)
                if minmax_general != []:
                    ax1.set_ylim(minmax_general)        
                
                ax1.xaxis.set_major_formatter(mdates.DateFormatter('%Y/%m/%d'))
                fig.autofmt_xdate()
                
                ax1.set_ylabel('general scale (blue)')
                ax2.set_ylabel('sensor maxi. amplitude scale  (green)')
                fig.suptitle( 'Sound speed (m/s) of sensor ' + str(iid) + ' at depth ' + str(np.mean(data[:,1])) + 'm., \n between ' + str(strt) + ' & ' + str(end))
                
                pathplot2 = os.path.join(pathplot,str(period_size),str(iid))
                gf.create_dir(pathplot2)
                if with_plot :
                    fig.savefig(os.path.join(pathplot2,gf.join_improved('_',iid,strt,end) + '.png' )) 
                    plt.close(fig.number)
                
                strt_ind = np.argwhere(bbool == True)[-1][-1] + 1
                strt     = data[:,0][strt_ind] 
                end      = data[:,0][strt_ind] + period_dt

# fabrication d'un statistic dico
# ie une valeur moyenne un std et une amplitude
# pour chaque jour & chaque senseur
if 0:
    statdayzbig = dict()
    period_dt   = dt.timedelta(days=1)

    linesallstk         = []
    
    for iid , data in dico.items():
        print("sensor " , iid)

        linesstk         = []
        
        strt = datetime.datetime.combine( data[0,0].date() , dt.time.min )
        end  = strt + period_dt
        
        minmax = [np.min(data[:,-1]) , np.max(data[:,-1])]
    
        while end <= data[-1,0]:
            print('strt , end : ' , iid , strt , end)
            bbool = (data[:,0] >= strt) * (data[:,0] < end)
            if not np.any(bbool):
                strt     = end
                end      = strt + period_dt
                print("SKIP !!!")
                continue
                
            dataplot = data[bbool,:]
            T = dataplot[:,0]
            C = dataplot[:,-1]
            
            date = list(geok.dt2ymdhms(geok.dates_middle(strt,end),0))      
            
            if len(C) < 10:
                strt     = end
                end      = strt + period_dt
                print("SKIP !!! bc len(C) < 10")
                continue
                
            line = [iid] + date + [np.mean(C), np.std(C), np.median(C), np.min(C), 
                           np.max(C), np.max(C) - np.min(C) , len(C)]

            linesstk.append( line )
            linesallstk.append( line )
            
            strt     = end
            end      = strt + period_dt
        
        finaltab         = np.vstack(linesstk)
        finaltabpd       = pd.DataFrame(finaltab)
        renamdic = gf.pandas_column_rename_dic(*['id','y','m','d','h',
        'mi','s','mean','std','median','min','max','amplitude','len'])
        finaltabpd = finaltabpd.rename(columns = renamdic)
        statdayzbig[iid] = finaltabpd

    finaltab         = np.vstack(linesallstk)
    finaltabpd       = pd.DataFrame(finaltab)
    renamdic = gf.pandas_column_rename_dic(*['id','y','m','d','h',
    'mi','s','mean','std','median','min','max','amplitude','len'])
    finaltabpd = finaltabpd.rename(columns = renamdic)
    statdayzbig['all'] = finaltabpd   
    gf.pickle_saver(statdayzbig,os.path.dirname(dicpath),'statistics_dico')

#%%
#### CETTE PARTIE EXPLOITE LES DONNÉES de l'exploit stat
if 1:
    #plt.ioff()
    
    TabStk = []
    std_residu_stk_stk = []
    
    datlis = [dt.datetime(2007, 5, 22, 0, 0),
              dt.datetime(2009, 10, 2, 0, 0),
              dt.datetime(2013, 3, 10, 0, 0)]
                     
    datlis = datlis + [dt.datetime(2010, 2, 25, 0, 12),
                       dt.datetime(2009, 6, 23, 0, 12),
                       dt.datetime(2007, 4, 24, 0, 12)]

    datlis =          [dt.datetime(2010, 2, 25, 0, 12) ,
                       dt.datetime(2009, 10, 2, 0, 0)  ,
                       dt.datetime(2007, 4, 24, 0, 12) ,
                       dt.datetime(2012, 12, 4,0,12)   ]
    
    nbdays = 1
        
        #    datlis = [dt.datetime(2007,1,1) + dt.timedelta(days=i) for i in np.arange(1,365)]
    
    goodbadstrlis = ["best 1","worst 1","median 1"] + ["best 2","worst 2","median 2"]
    goodbadstrlis = ["best","worst","median","median 2"]
    
    
    outdicolis = []
    
    
    for strt,goodbadstr in zip(datlis,goodbadstrlis):
        Dstk  = []
        Dstk2 = []
        std_residu_stk  = []
        std_global_stk  = []
    
        # 1st identification mode seems not the best
        #    best
        #    2007-05-23 00:12:00
        #    worst
        #    2009-10-02 00:12:00
        #    median
        #    2013-03-10 00:12:00
    
        # 2nd identification mode, maybe best
        #>>> Tbest
        #datetime.datetime(2010, 2, 25, 0, 12)
        #>>> Tworst
        #datetime.datetime(2009, 6, 23, 0, 12)
        #>>> Tmedian
        #datetime.datetime(2007, 4, 24, 0, 12)
        
        # median en utilisant la 3eme methode 
        # somme des 5 premiers senseurs
        # memes max et min
        # 2012-12-04 00:12:00 7.25803353652
        # datetime.datetime(2012, 12, 4,0,12)
        end = strt + dt.timedelta(days=nbdays)
    
        plt.ion()
        
        # =========== PLOT of C(time) / depth as colorbar ==============
        fig , ax = plt.subplots()    
          
        def fctcossin(x,kk,ks,kc,kp):
            return kk * x + (np.sin(ks*x) + np.cos(kc*x))^kp
    
        def fctsin2(x,A,w,p):
            return A * np.sin(w*x + p)
              
        for k , data in dico.items(): 
            bbool    = (data[:,0] >= strt) * (data[:,0] < end)
            datatrue = data[bbool,:]
            Dstk2 = Dstk2 + list(set(datatrue[:,1]))
            
        sm = plt.cm.ScalarMappable(cmap='winter_r',
                                   norm=plt.Normalize(vmin=np.min(Dstk2), 
                                                      vmax=np.max(Dstk2)))
        cbar   = sm.set_array(Dstk2)
        cbarax = fig.colorbar(sm) #,cax=cbar)
        cbarax.ax.invert_yaxis()
        
        cbarax.ax.set_title('depth',size='medium')
    
        ax.set_xlabel('time')
        ax.set_ylabel('Sound Speed (m/s)')
        fig.suptitle('Temporal Sound Speed Profile \n of MOVE3 for day ' + str(strt.date()) + ' (\"' + goodbadstr + '\")')
        
        D_4_graph2 = []
        C_4_graph2 = []
        T_4_graph2 = []
        
        outdico = dict()
        outdicolis.append(outdico)
        
        outdico['date'] = strt       
    
        for k in reversed(list(dico.keys())):
            
            data = dico[k]
            
            bbool = (data[:,0] >= strt) * (data[:,0] < end)
            datatrue = data[bbool,:]
            
            T = datatrue[:,0]
            D = datatrue[:,1]
            C = datatrue[:,2]        
            Dstk = Dstk + list(np.unique(D))
            
            D_4_graph2.append(D[0])
            C_4_graph2.append(C)
            T_4_graph2.append(T)
            
            Tposix = np.array(geok.dt2posix(T))
            Tposixreduc = Tposix - Tposix[0]
            CD = C + k
            
            C = np.array(C, dtype='float')
            
            Poly = np.polyfit(Tposixreduc,C,15)                
            Cpoly = np.polyval(Poly,Tposixreduc)
            std_residu_stk.append(np.std(C - Cpoly))
            std_global_stk.append(np.std(C))
    
            print('sensor' , k , std_residu_stk[-1] , 'm/s')
            
            coul = sm.to_rgba(np.unique(D)[0])
            
            ax.plot(T,C    ,'+' ,c=coul)
            if 0 : # ici on plot le polynome
                ax.plot(T,Cpoly,'x-',c=coul)

        outdico['D'] = D_4_graph2
        outdico['C'] = C_4_graph2
        outdico['T'] = T_4_graph2

        std_residu_stk_stk.append(std_residu_stk)

        
        my_dpi=96
        fig.set_dpi(my_dpi)
        fig.set_size_inches(800/my_dpi, 600/my_dpi)
        fig.tight_layout()
        plt.subplots_adjust(top=0.90)
                
        for ext in ('.svg','.pdf','.png'):
            pathplot2 = os.path.join(dirplot2,"MOVE3_varia_" + str(strt.date()) + ext)
            if with_plot:
                plt.savefig(pathplot2)
            
        stdstk = std_residu_stk        
        stdstk = std_global_stk
         
        bound_zone_lis = list(reversed(gf.middle(Dstk)))
        sigma_zones    = list(reversed(stdstk))
    
        # =========== Moustache box ==============
        fig , ax  = plt.subplots()
        Cmean_plt =   np.array([np.mean(CC) for CC in C_4_graph2])
        D_plt     = - np.array(D_4_graph2)
        ax.plot(Cmean_plt , D_plt)
        ax.boxplot(list(C_4_graph2) , vert=0 , 
                   positions = - np.array(D_4_graph2) , widths = 200)
        ax.set_xlim((1485,1550))        
        ax.set_ylim((-5100,100))        
        ax.set_yticks(np.arange( -5100,100, 500))
        ax.set_yticklabels([str(e) for e in np.arange( -5100,100, 500)])
        ax.set_ylabel('Depth (m)')
        ax.set_xlabel('Sound Speed (m/s)')
        if 0:
            fig.suptitle('Variation of the Sound Speed Profile \n of MOVE3 for day '  + str(strt.date()) + ' (\"' + goodbadstr + '\")')
        else:
            print("SIMPLIFIED TITLE")
            fig.suptitle('Variations of the Sound Speed Profile')

        my_dpi=96
        fig.set_dpi(my_dpi)
        fig.set_size_inches(800/my_dpi, 600/my_dpi)
        fig.tight_layout()
        plt.subplots_adjust(top=0.90)

        for ext in ('.svg','.pdf','.png'):
            dirplot2  = gf.create_dir(dirplot2)
            pathplot2 = os.path.join(dirplot2,"MOVE3_varia_moustachbox_" + str(strt.date()) + ext)
            if with_plot:
                plt.savefig(pathplot2)
#%%           
        #plt.close('all')

        # =========== PLOT of Depth(time) / SV - SV mean as colorbar ==============

        fig , ax = plt.subplots()    
          
        Cminmaxstk = []
        for k , data in dico.items(): 
            bbool    = (data[:,0] >= strt) * (data[:,0] < end)
            datatrue = data[bbool,:]
            Dstk2    = Dstk2 + list(set(datatrue[:,1]))
            dC       = datatrue[:,2] - np.mean(datatrue[:,2])
            Cminmaxstk.append((np.min(dC),np.max(dC)))
        
        Cminmaxstk = np.vstack(Cminmaxstk)
        
        cmintot = np.min(Cminmaxstk[:,0])
        cmaxtot = np.max(Cminmaxstk[:,1])
        
        cextrema     = np.max(np.abs(Cminmaxstk))
        cextremamarg = cextrema + cextrema *.1
        
        sm = plt.cm.ScalarMappable(cmap='coolwarm_r',
                                   norm=plt.Normalize(vmin=-cextremamarg, 
                                                      vmax=cextremamarg))
        cbar   = sm.set_array(Dstk2)
        cbarax = fig.colorbar(sm) #,cax=cbar)
        cbarax.ax.invert_yaxis()
        
        cbarax.ax.set_title('depth',size='medium')
    
        ax.set_xlabel('time')
        ax.set_ylabel('Sound Speed (m/s)')
        fig.suptitle('Temporal Sound Speed Profile \n of MOVE3 for day ' + str(strt.date()) + ' (\"' + goodbadstr + '\")')
        ax.set_yscale("log")
        
        D_4_graph2 = []
        C_4_graph2 = []
    
        for k in reversed(list(dico.keys())):
            
            data = dico[k]
            
            bbool = (data[:,0] >= strt) * (data[:,0] < end)
            datatrue = data[bbool,:]
            
            T = datatrue[:,0]
            D = datatrue[:,1]
            C = datatrue[:,2]        
            Dstk = Dstk + list(np.unique(D))
            
            D_4_graph2.append(D[0])
            C_4_graph2.append(C)
            
            Tposix = np.array(geok.dt2posix(T))
            Tposixreduc = Tposix - Tposix[0]
            CD = C + k
            
            print('sensor' , k , std_residu_stk[-1] , 'm/s')
            
            for tt,dd,cc in zip(T,D,C):
                coul = sm.to_rgba(cc - np.mean(C))
                ax.scatter(tt,dd,marker='+',c=coul)
        
        ax.set_xlim((np.min(T),np.max(T)))
        ax.invert_yaxis()

        my_dpi=96
        fig.set_dpi(my_dpi)
        fig.set_size_inches(800/my_dpi, 600/my_dpi)
        fig.tight_layout()
        plt.subplots_adjust(top=0.90)

                
        for ext in ('.svg','.pdf','.png'):
            pathplot2 = os.path.join(dirplot2,"MOVE3_varia_" + str(strt.date()) + ext)
            if with_plot:
                plt.savefig(pathplot2)
            
        stdstk = std_residu_stk        
        stdstk = std_global_stk
         
        bound_zone_lis = list(reversed(gf.middle(Dstk)))
        sigma_zones    = list(reversed(stdstk))



        # =========== TABLE ==============    
        Tab = []
    
        for k in list(dico.keys()):
            data     = dico[k]
            bbool    = (data[:,0] >= strt) * (data[:,0] < end)
            datatrue = data[bbool,:]
            
            T = datatrue[:,0]
            D = datatrue[:,1]
            C = datatrue[:,2]    
            
            Tab.append((k,D[0] , np.mean(C) , np.std(C) ,  np.max(C) - np.min(C) ))

        print(strt)
        head = [str(strt.date()),'Profondeur \n (m)','SV moyen \n (m/s)','Ecart type \n (m/s)','Amplitude (m/s)']
        TabStk.append(tabulate.tabulate(Tab,headers=head,tablefmt='latex',
                                        floatfmt=".3f"))
    
    for t in TabStk:
        print(t)
    
    
    gf.pickle_saver(outdicolis,'/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/scripts/analyz_SSP_MOVE_NOAA','SSPT_MOVE_DICLIS')

# ======= TRAVAIL SUR LES COURENTS =======

#%%
from megalib import *
if 1:
    OUTDICLIS = gf.pickle_loader('/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/scripts/analyz_SSP_MOVE_NOAA/SSPT_MOVE_DICLIS.pik')
    COUR      = [0.06190411, 0.15694521, 0.16195206, 0.24940379]
    COURnam   = ["best","worst","median","median 2"]
    CURRENTS_COLUMNS_4_INTERPO = gf.pickle_loader("/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/scripts/analyz_SSP_MOVE_NOAA/CURRENTS_INTERPOLATORS.pik")
    
    T_interpo_mini_max = []
    CURRENTS_INTERP = []
    for ccc in CURRENTS_COLUMNS_4_INTERPO:
        print('date interpolators' , geok.posix2dt(np.mean(ccc[:,0]),1))
        CURRENTS_INTERP.append(scipy.interpolate.interp1d(ccc[:,0] , ccc[:,1]))  
        T_interpo_mini_max.append( ( np.min(ccc[:,0]) , np.max(ccc[:,0]) ) )
    
    for iexp , ( dic , cour , cour_interp , cournam , cour_col_4_interpo , t_minmax) in enumerate(zip(OUTDICLIS,COUR , CURRENTS_INTERP , COURnam , CURRENTS_COLUMNS_4_INTERPO , T_interpo_mini_max)):
        D = np.array(dic['D'])
        T = np.array(dic['T'])
        C = np.array(dic['C'])

        #print 'date data' , geok.posix2dt(np.mean(geok.dt2posix(T,1)),1)
        
        print("new day")
         
        dCdLmean_stk = []
        dCdL_stk     = []
        
        for CC , TT , d  in zip(C,T,D):
            print([geok.posix2dt(t_minmax[0]) < tt < geok.posix2dt(t_minmax[1]) for tt in TT])

            cour_interpoled = cour_interp(geok.dt2posix(TT,1))
            print(np.mean(cour_interpoled))
            

            dLL      = np.array([e.seconds for e in np.diff(TT)]) * cour  #cour_interpoled[:-1]
            dCC      = np.diff(CC)
            dCdL     = dCC / dLL
            dCdLmean = np.mean(dCdL)  
            dCdLmean_stk.append(dCdLmean)
            dCdL_stk.append(dCdL)

            print(d , np.mean(dCdL) , np.median(dCdL)  , np.std(dCdL))
            
            #plt.plot(dCdT_stk , -np.array([d] * len(dCdT_stk)) )
        
        #plt.figure()
plt.plot(dCdLmean_stk ,-np.array(D),label=cournam)
plt.legend()
plt.xlim((-0.0002,0.0002))

# Reprise du code ATALANTE2003 pour détermination du gradient "polynomé"

Zgrad = D
Grad  = np.array(dCdLmean_stk)

Grad[np.isnan(Grad)] = 0
Grad[Zgrad >= 4000] = 0

polykoefs = np.polyfit(Zgrad,Grad,10)
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
plt.plot( Grad      ,   -Zgrad   , '*b')
plt.plot(  GradPoly , -Zgradpoly , '+g')
plt.plot(GradPoly2  , -Zgradpoly  , 'xk')
plt.xlim((-0.0002,0.0002))


fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)

ZGrad_arr = np.column_stack((Zgrad, Grad))    
ZGrad_arr = np.column_stack((Zgradpoly, GradPoly2))    

print("strt save")
outdirKfile  = '/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/exemple/Kfiles_kourents_gradients'
outnameKfile = 'kfile_poly_proxy_MOVE.K'
outpathKfile = os.path.join(outdirKfile,outnameKfile)
np.savetxt(outpathKfile,ZGrad_arr)

plt.savefig(outdirKfile + '/plotgrad.png')        
