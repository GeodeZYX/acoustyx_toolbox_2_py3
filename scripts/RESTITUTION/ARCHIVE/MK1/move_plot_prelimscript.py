# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 14:07:17 2015

@author: psakicki
"""

import glob
import numpy as np
import matplotlib.pyplot as plt
import geodetik as geok
from scipy.fftpack import fft
import scipy

import MOVEclass as mcls

path= '/home/psakicki/THESE/bd_CTD/NEW_DATA_2015/MOVE/MOVE3_clean/ssp*'
path_fig = '/home/psakicki/THESE/bd_CTD/NEW_DATA_2015/MOVE/MOVE3_clean/'
fil_lis = glob.glob(path)


mcls.read_multi_ssp_file(path_fig,'ssp*')




#fil = fil_lis[0]
#M = np.loadtxt(fil)
#T = M[:,0]
#
#tfin = np.min(T) + 31557600 
#
#Myear = M[M[:,0] < tfin]
#
#Mdiff = np.diff(Myear[:,0])
#
#
#y = fft(Mdiff)
#
#plt.plot(np.abs(y))


mean_lis = []
std_lis = []
prof_lis = []
prof_mean_lis = []

f = plt.figure()
for i,fil in enumerate(fil_lis):
    f.clf()
    M = np.loadtxt(fil)
    ax = plt.gca()
    ax.ticklabel_format(useOffset=False)
    prof = list(set(list(M[:,2])))
    prof_lis.append(prof)
    prof_mean_lis.append(np.mean(prof))
    str_prof = ''
    for p in prof:
        str_prof = str_prof + str(p) + ' ' 
#    ax.plot(geok.posix2dt(M[:,0]),M[:,1],'.')
    f.suptitle('sensor ' + str(i+1) + 'prof : ' + str_prof)
    f.savefig(path_fig + 'plot_ssp' + str(i) + '.pdf')
    mean_lis.append(np.nanmean(M[:,1]))
    std_lis.append(np.std(M[:,1]))
#    plt.close(f)
    
mean_lis = [m for (p,m) in sorted(zip(prof_mean_lis,mean_lis))]
std_lis = [m for (p,m) in sorted(zip(prof_mean_lis,std_lis))]
prof_mean_lis = sorted(prof_mean_lis)
    
plt.plot(mean_lis,-np.array(prof_mean_lis))
plt.errorbar(mean_lis,-np.array(prof_mean_lis) , xerr=np.array(std_lis)*1)