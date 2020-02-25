# -*- coding: utf-8 -*-
"""
Created on Sat Aug 27 11:10:44 2016

@author: psakicki
"""

from megalib import *

from netCDF4 import Dataset

path = "/home/psakicki/THESE/bd_CTD/NEW_DATA_2015/MOVE/MOVE3/OTHERS/*CURRENT*nc"
FL = glob.glob(path)[-1:]
#FL = glob.glob(path)
Tjjulstk = []

Kstk = []


for f in FL:
    fh = Dataset(f, mode='r')
    Depth = fh['DEPTH'][:]
    T = geok.jjulCNES2dt(fh['TIME'][:])
    Tjjul = (fh['TIME'][:])
    K = fh['CSPD'][:]
    K2 = K[:,1]
    K1 = K[:,0]
    D = fh['CDIR'][:]
    D2 = D[:,1]
    D1 = D[:,0]
    
    Depth1 =  str(np.round(Depth[0])) + 'm'
    Depth2 =  str(np.round(Depth[1])) + 'm'
    



#T = np.concatenate(Tjjulstk)
#K = np.concatenate(Kstk)
#
#T = np.array(geok.jjulCNES2dt(T))


T1,K1 = gf.sort_binom_list(T,K1,1)
T2,K2 = gf.sort_binom_list(T,K2,1)

T1,D1 = gf.sort_binom_list(T,D1,1)
T2,D2 = gf.sort_binom_list(T,D2,1)


plt.plot(T1,K1,label=Depth1)
plt.plot(T2,K2,label=Depth2)
plt.figure()

plt.plot(T1,D1,label=Depth1)
plt.plot(T2,D2,label=Depth2)


# Filter requirements.
order = 6
fs     = 1 # sample rate, Hz
cutoff = 0.00125 # desired cutoff frequency of the filter, Hz


K1f = butter_lowpass_filter(K1, cutoff, fs,order)
K2f = butter_lowpass_filter(K2, cutoff, fs,order)
D1f = butter_lowpass_filter(D1, cutoff, fs,order)
D2f = butter_lowpass_filter(D2, cutoff, fs,order)

outplotdir = "/home/psakicki/Documents/manuscrit/FIG/oceano/"

fig = plt.figure()
plt.plot(T1,K1f,label=Depth1)
plt.plot(T2,K2f,label=Depth2)
plt.ylabel('Current speed (m/s)')
plt.legend(title='Depth')
fig.autofmt_xdate()
plt.savefig(outplotdir + 'kouran_spd.pdf')
fig = plt.figure()
plt.plot(T1,D1f,label=Depth1)
aaa = plt.plot(T2,D2f,label=Depth2)
plt.ylabel('Current direction (°)')
plt.legend(title='Depth')
fig.autofmt_xdate()
plt.savefig(outplotdir + 'kouran_dir.pdf')
fig.autofmt_xdate()

        
    
