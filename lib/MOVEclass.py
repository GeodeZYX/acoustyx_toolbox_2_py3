# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 17:30:59 2015

@author: psakicki
"""
import glob
import numpy as np
import matplotlib.pyplot as plt
import geodetik as geok
from scipy.fftpack import fft
import scipy
import glob
import os
import SSP as ssp
import datetime as dt
import genefun
import matplotlib
import sys

sys.dont_write_bytecode = True


from scipy.optimize import curve_fit

class Mesure(object):
    def __init__(self):
        self.time = np.nan
        self.pres = np.nan
        self.depth = np.nan
        self.temp = np.nan
        self.sali = np.nan
        
        self.bool_valid = None
        
        self.campagn = np.nan
        self.sensor  = np.nan
        
        self.sndspd = np.nan
        
    def calc_sndspd(self):
        if self.bool_valid:
            self.sndspd = ssp.soundspeed(self.sali, self.temp, self.pres, 'del_grosso')     
        
class Sensor(object):
    def __init__(self,inid = np.nan):
        self.id = inid
        self.depth = np.nan
        self.Data = []
        self.T    = []
        
        self.T_uptodate = False
        
        self.__start = dt.datetime(1980,1,1)
        self.__end   = dt.datetime(2099,1,1)
        self.bool_start_updt = False
        self.bool_end_updt   = False
        self.bool_sorted     = False
        self.bool_min_max    = False
        self.__min_sv = np.nan
        self.__max_sv = np.nan

    @property
    def start(self):
        if not self.bool_start_updt:
            self.__start = np.min([e.time for e in self.Data])
            self.bool_start_updt = True
        return self.__start
    
    @property
    def end(self):
        if not self.bool_end_updt:
            self.__end = np.max([e.time for e in self.Data])
            self.bool_end_updt = True
        return self.__end

    @property
    def max_sv(self):
        if not self.bool_min_max:
            self.__max_sv = np.max(self.get_list('sndspd'))
            self.__min_sv = np.min(self.get_list('sndspd'))
            self.bool_min_max = True
        return self.__max_sv

    @property
    def min_sv(self):
        if not self.bool_min_max:
            self.__max_sv = np.max(self.get_list('sndspd'))
            self.__min_sv = np.min(self.get_list('sndspd'))
            self.bool_min_max = True
        return self.__min_sv
        
    def updateT(self):
        self.T = []
        self.T_uptodate = False

        for d in self.Data:
            self.T.append(geok.dt2posix(d.time))
        self.T = np.array(self.T)
        self.T_uptodate = True
        return None
            
    def __getitem__(self,ind):
        return self.Data[ind]
        
    def title_plot_str(self):
        s_str = self.start.strftime('%Y/%m/%d-%H:%M:%S')
        e_str = self.end.strftime('%Y/%m/%d-%H:%M:%S')
        return 'sens. ' + str(self.id) + ' depth ' + str(self.depth) + 'm , ' + s_str + ' => ' + e_str
    
    def filename_str(self):
        s_str = self.start.strftime('%Y%m%d%H%M%S')
        e_str = self.end.strftime('%Y%m%d%H%M%S')   
        
        daydelta = (self.end - self.start).days
        return '_'.join((str(self.id) , s_str ,  e_str , str(daydelta) ))
        
    def append_data(self,indata):
        self.Data.append(indata)
        self.bool_start_updt = False
        self.bool_end_updt   = False
        self.bool_min_max    = False
        if np.isnan(self.id):
            self.id = indata.sensor
        else:
            if self.id != indata.sensor:
                print("WARN : ça chie sur l id ...")
                print(self.id , indata.sensor)
                
    def get_list(self,datatype):
        out = [getattr(m,datatype) for m in self.Data]
        return out
        
    def get_TZC(self):
        T = np.array(self.get_list('time'))
        Z = np.array(self.get_list('depth'))
        C = np.array(self.get_list('sndspd'))
        return T,Z,C
        
                
    def check(self):
        datatyp_lis = ['pres','temp','sali'] #depth
        boolis_final = np.ones(len(self.Data))
        for dtype in datatyp_lis:
            _ , boollis = geok.outiler_mad(self.get_list(dtype))

            boolis_final = boolis_final * boollis
         
#        for i in range(len(boollis)):
#            boolis_final[i] = boolis_final[i] * boollis[i]
        
        for b , d in zip(boolis_final , self.Data):
            d.bool_valid = b
        return None
        
    def clean(self):
        for m in self.Data:
            if m.bool_valid  == False:
                self.Data.remove(m)
                self.bool_start_updt = False
                self.bool_end_updt   = False
                self.bool_min_max    = False
        return None
        
    def stats(self):
        ssp = self.get_list('sndspd')
        sspmean = np.mean(ssp)
        sspstd  = np.std(ssp)
        sspmedian = np.median(ssp)
        
        return sspmean , sspstd , sspmedian , self.min_sv , self.max_sv
        
    def calc_sndspd(self):
        for m in self.Data:
            m.calc_sndspd()
        return None
        
    def plot(self,f=0):
        if f == 0:
            f = plt.figure()
        prof_lis  = list(set([e.depth for e in self.Data]))
        prof_mean = np.nanmean(prof_lis)
        
        T = [e.time for e in self.Data]
        SSP = [e.sndspd for e in self.Data]
        
        plt.plot(T,SSP,'-')
        plt.plot(T,SSP,'+r')
        
        return None
        
    def diffT(self):
        return list(set(np.diff(geok.dt2posix([e.time for e in self.Data]))))
        
    def sortT(self):
        self.T_uptodate = False
        self.Data.sort(key=lambda x: x.time)      
        self.bool_sorted = True
        
    def get_windowed_epoch_list(self,start=dt.datetime(1980,1,1),
                                end=dt.datetime(2099,1,1)):
        T = [e.time for e in self.Data]
        T.sort()
        _,i = genefun.find_nearest(T,start)
        _,j = genefun.find_nearest(T,end)
        return T[i:j]
    
class SensorGroup(object):
    def __init__(self):
        self.grp = []
        self.id_list = []
        
    def __call__(self,inid):
        out = None
        
        if inid in self.id_list:
            try:
                if self.grp[inid-1].id == inid:
                    out = self.grp[inid-1]
                else:
                    for e in self.grp:
                        if e.id == inid:
                            out = e                    
            except:
                for e in self.grp:
                    if e.id == inid:
                        out = e
                                        
        if out == None:
            e = Sensor(inid)
            self.append_sensor(e)
            out = e
        
        return out
        
    def append_sensor(self,sens_in):
        self.grp.append(sens_in)
        self.grp.sort(key=lambda x: x.id)
        self.id_list.append(sens_in.id)
        return None
        
    def append_mes_in_rigth_sens(self,mes_in):
        self(mes_in.sensor).append_data(mes_in)
        return None
        
    def get_mesure_epoch(self,epoch,marge=dt.timedelta(minutes=5)):
        out_mes_list = []
        guessmode = False
        for sens in self.grp:
#            if guessmode:
#                iguess = i - 1000
#                jguess = i + 1000
#            else:
#                iguess = 0
#                jguess = len(sens.Data)  - 1
                
            if sens.bool_sorted == False:
                sens.sortT()
                
            T = sens.get_list('time')
#            _ , i = genefun.find_nearest(T[iguess:jguess],epoch)
            _ , i = genefun.find_nearest(T,epoch)

            if (T[i] - epoch) > marge:
                continue
            try:
                out_mes_list.append(sens.Data[i])
#                guessmode = True
            except:
#                guessmode = False
                continue
        return out_mes_list
        
        
    def get_mesure_epoch2(self,epoch,marge=300):
        out_mes_list = []
        guessmode = False
        
        for sens in self.grp:

            if not sens.bool_sorted:
                sens.sortT()
                
            if not sens.T_uptodate:
                sens.updateT()
            
            posixepoch = geok.dt2posix(epoch)
            nearT , inearT = genefun.find_nearest(sens.T,posixepoch)


            if (nearT - posixepoch) > marge:
                continue

            try:
                out_mes_list.append(sens.Data[inearT])
#                guessmode = True
            except:
#                guessmode = False
                continue
        return out_mes_list
        
        
class Sensor_period(Sensor):
    def __init__(self,idin=np.nan):
        Sensor.__init__(self,idin)

    @property
    def start(self):
        return np.min([e.time for e in self.data])
    
    @property
    def end(self):
        return np.max([e.time for e in self.data])
    

def read_ssp_file(pathin,idin = 0): 
    """
    Read a file
    POSIX TIME SV DEPTH
    In a sensor object
    """
    print('file',pathin)
    if idin == 0:
        idin = int(pathin.split('id')[1].split('.')[0])
    
    Sens = Sensor(idin)
    print(Sens.id)
    M = np.loadtxt(pathin)
    p_mean = np.nanmean(M[:,2])
    Sens.depth = p_mean
        
    for i in range(M.shape[0]):
        t = geok.posix2dt(M[i,0])
        ss = M[i,1]
        p = M[i,2]
        Mes = Mesure()
        Mes.time = t
        Mes.sndspd = ss
        Mes.depth = p
        Mes.bool_valid = True
        Mes.sensor = idin
        Sens.append_data(Mes)
        
    return Sens
    
def read_multi_ssp_file(pathin,wildcardin):
    fil_lis = sorted(glob.glob(os.path.join(pathin,wildcardin)))
    Sensgrp = SensorGroup()
    for fil in fil_lis:
        sens = read_ssp_file(fil)
        Sensgrp.append_sensor(sens)
    return Sensgrp
        
def decimate_ssp_file(filein,start,end):   
    outputprefix = 'decssp_'
    fildir = os.path.dirname(filein)
    filnam = os.path.basename(filein)
    ID = int(filnam.split('_')[-1])
    M = np.loadtxt(filein)
       
    start_posix = geok.dt2posix(start)
    end_posix = geok.dt2posix(end)
    
    bol = (start_posix <= M[:,0]) * (M[:,0] < end_posix)
        
    outname = outputprefix + start.strftime("%Y%m%d") +'_' \
    + end.strftime("%Y%m%d") + '_' + str(ID)
    
    outpath = os.path.join(fildir,outname)
    np.savetxt(outpath,M[bol])
    return None

def plot_list_of_mesures(listin,fig=0):
    P   = [m.depth for m in listin]
    SSP = [m.sndspd for m in listin]
    
    fig = genefun.get_figure(fig)
    ax = plt.gca()
    
    ax.plot(SSP,-np.array(P)) 
        

    

def search_period_sensor(sensor_in,start,end,indstrt=0):
    """
    Make a new sensor object windowed from a Raw Sensor object
    """
    out_sens_period = Sensor()
    if not sensor_in.bool_sorted:
        sensor_in.sortT()
        
    datalist = sensor_in.Data
    datalist = datalist[indstrt:]        

    for i,data in enumerate(datalist):
        if data.time >= end:
            break
        elif start <= data.time < end:
            out_sens_period.append_data(data)
        else:
            continue
        
    out_sens_period.depth = sensor_in.depth # c'est pas normal de devoir fair ca mais au moins on est sur
    return out_sens_period , i


def search_multi_period_sensor(sensor_in,deltatime):
    """
    from a raw sensor (with all data)
    make several sensors with the data of the period inside
    """
        
    start = sensor_in.start
    end   = sensor_in.end
          
    if not sensor_in.bool_sorted:
        sensor_in.sortT()
        
    curstrtdate = start
    indstart = 0

    outsens_lis = []
    while curstrtdate < end:
        curenddate = curstrtdate + deltatime
        tempsens , inext = search_period_sensor(sensor_in,curstrtdate,curenddate,indstart)
        outsens_lis.append(tempsens)
        curstrtdate = curenddate
        indstart = inext
    return outsens_lis


def plot_TC_frontend(sensorin,close=True,path='',minmax=[],
                     minmax_general=[],plotmod=2,outputext = ('.png',)):
    """
    plot mode :
    
    1 = old style => the time window scale in blue + general sensor scale in green
    2 = general sensor scale in green + whole spatial temporal scale in blue
    """
    import matplotlib.dates as mdates

    if sensorin.Data  == []:
        return None

    plt.ioff()
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    T,Z,C = sensorin.get_TZC()

    y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    ax1.yaxis.set_major_formatter(y_formatter)
    ax2.yaxis.set_major_formatter(y_formatter)
    fig.suptitle(sensorin.title_plot_str())
    fig.set_size_inches(16.53,11.69)  
    stattxt = ['{:12.6f}'.format((e)) for e in sensorin.stats()]
    plt.text( 0.00, -0.18 , stattxt , transform = ax1.transAxes)

    ax1.plot(T,C,'.b')
    ax2.plot(T,C,'.g')

    if plotmod == 1:
        ax1.set_ylabel('daily scale (blue)')
        ax2.set_ylabel('sensor maxi. amplitude scale  (green)')
        if minmax != []:
            ax2.set_ylim(minmax)
    elif plotmod == 2:
        ax1.set_ylabel('general scale (blue)')
        ax2.set_ylabel('sensor maxi. amplitude scale  (green)')
        if minmax != []:
            ax2.set_ylim(minmax)
        if minmax_general != []:
            ax1.set_ylim(minmax_general)        
    
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%Y/%m/%d'))
    fig.autofmt_xdate()

    if path != '':
        for ext in outputext:
            fig.savefig(os.path.join(path,sensorin.filename_str() + ext )) 
    if close:
        plt.close(fig.number)

    return None


def plot_TC_frontend_OLD(sensorin,close=True,path='',minmax=[]):
    import matplotlib.dates as mdates

    if sensorin.Data  == []:
        return None

    plt.ioff()
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    T,Z,C = sensorin.get_TZC()
    ax1.plot(T,C,'.b')
    ax2.plot(T,C,'.g')
    y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    ax1.yaxis.set_major_formatter(y_formatter)
    ax2.yaxis.set_major_formatter(y_formatter)
    fig.suptitle(sensorin.title_plot_str())
    fig.set_size_inches(16.53,11.69)  
    stattxt = ['{:12.6f}'.format((e)) for e in sensorin.stats()]
    plt.text( 0.00, -0.18 , stattxt , transform = ax1.transAxes)

    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%Y/%m/%d'))
    fig.autofmt_xdate()
    ax1.set_ylabel('daily scale (blue)')
    ax2.set_ylabel('sensor maxi. amplitude scale  (green)')

    if minmax != []:
        ax2.set_ylim(minmax)
    if path != '':
        fig.savefig(os.path.join(path,sensorin.filename_str() + '.png' )) 
        fig.savefig(os.path.join(path,sensorin.filename_str() + '.pdf' )) 
    if close:
        plt.close(fig.number)

    return None

def ZC_from_mesure_lis(mes_lis_in):
    SSP = np.array([m.sndspd for m in mes_lis_in])
    Z   = np.array([m.depth  for m in mes_lis_in])
    return Z , SSP

  #_____                   _       ______        _               _   
 #|_   _|                 | |     / / __ \      | |             | |  
   #| |  _ __  _ __  _   _| |_   / / |  | |_   _| |_ _ __  _   _| |_ 
   #| | | '_ \| '_ \| | | | __| / /| |  | | | | | __| '_ \| | | | __|
  #_| |_| | | | |_) | |_| | |_ / / | |__| | |_| | |_| |_) | |_| | |_ 
 #|_____|_| |_| .__/ \__,_|\__/_/   \____/ \__,_|\__| .__/ \__,_|\__|
             #| |                                   | |              
             #|_|                                   |_|              
#

#
#def load_raw_move_files_in_sensgrp_obj(filelistin):
#    M_lis = []
#    sens_grp = mcls.SensorGroup()
#    for fil in fillist:
#        
#        bool_part = False
#        
#        fil_name   = os.path.basename(fil)
#        id_campagn = int(fil_name.split('_')[2])
#        
#        if 'PART' in fil:
#            bool_part = True
#            continue
#              
#        d = Dataset(fil)
#        
#        TIME = [dt.datetime(1950,1,1) + dt.timedelta(days=t) for t in d.variables['TIME'][:]]
#        
#        SALINITY = np.squeeze(np.array(d.variables['PSAL'][:]))
#        PRESURE  = np.squeeze(np.array(d.variables['PRES'][:]))
#        TEMP     = np.squeeze(np.array(d.variables['TEMP'][:]))
#        if not bool_part:
#            DEPTH = np.squeeze(np.array(d.variables['DEPTH'][:]))
#        else:
#            DEPTH = np.array(d.variables['DEPTH'][:])
#        
#        if not bool_part:
#            for di,d in enumerate(DEPTH):
#                for ti,t in enumerate(TIME):
#                    M = mcls.Mesure()
#                    M.depth = d
#                    M.pres  = PRESURE[ti,di] 
#                    M.temp  = TEMP[ti,di]
#                    M.sali  = SALINITY[ti,di] 
#                    M.time = t
#                    
#                    M.campagn = id_campagn
#                    M.sensor = di+1
#                    
#                    sens_grp(M.sensor).append_data(M)
#    #                sens_grp.append_mes_in_rigth_sens(M)
#                  
#    return sens_grp
#
#
#
#
#
#def write_sensgrp_obj_in_files(sens_grp):
#    for sens in sens_grp.grp:
#        print sens.id
#        sens.check()
#        sens.calc_sndspd()
#        
#        filobj = open(os.path.join(path,'data_'+str(sens.id)),'w')
#        filobj2 = open(os.path.join(path,'ssp_'+str(sens.id)),'w')
#        
#        for m in sens.Data:
#            final_str = ''
#            for a in ['time' , 'sensor' ,  'depth' , 'pres' , 'sali' , 'sndspd' , 'bool_valid']:
#                final_str = final_str + str(getattr(m,a)) + ' '
#            filobj.write(final_str + '\n')
#            if m.bool_valid:
#                filobj2.write(str(geok.dt2posix(m.time)) + ' ' + str(m.sndspd) + ' ' + str(m.depth) +   '\n')
#
#        filobj.close()
#        filobj2.close()
#    return None




