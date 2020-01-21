#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 17 22:40:49 2017

@author: adminuser
"""

from megalib import *

export_path = '/home/adminuser/aaa_FOURBI/'

export_path = '/home/adminuser/aaa_FOURBI/'
Gpath = export_path + 'Gradtensor' + '.npy' 
Zpath = export_path + 'Ztensor' + '.npy'
Xpath = export_path + 'Xtensor' + '.npy'
Ypath = export_path + 'Ytensor' + '.npy'
Cpath = export_path + 'Ctensor' + '.npy'

C = np.load(Cpath)
X = np.load(Xpath)
Y = np.load(Ypath)
Z = np.load(Zpath)
G = np.load(Gpath)

for i in (0,500,1000,2500,4000,5999):
    if i == 5999:
        j = i +1
    else:
        j = i
    plt.plot( Y[0,:,i] ,G[:,0,i] , label=str(j) + 'm')
    
plt.ylabel('sound speed (m/s)')
plt.xlabel('distance to the center along y component (m)')

plt.legend(title='Depth')

