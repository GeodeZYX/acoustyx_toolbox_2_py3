# -*- coding: utf-8 -*-
"""
Created on Sun Aug 16 10:43:27 2015

@author: pierre
"""

from megalib import *
from natsort import natsorted, ns

start = datetime.datetime(2007, 1, 4, 19, 20)
start_posix = geok.dt2posix(start)
len_sspt = 84600

sspfillis = glob.glob('/media/pierre/D34F-24E3/plots/SSPs/*ssp.dat')

# lecture brute
#sspdic = acls.sspfiles_list2sspdic(sspfillis)
# lecture d'un dico
sspdic = genefun.pickle_loader('/media/pierre/D34F-24E3/plots/SSPs/sspdic')

I , Zuniq = acls.sspdic2InterpoSSPT(sspdic,start_posix,len_sspt)

# un Z manuel eventuellement utile
Z = np.arange(31.5,  4000 , 1)

C = acls.SSPT_from_Interpo(I,Zuniq,23)

plt.clf()
plt.plot(C,-Zuniq,'+')









