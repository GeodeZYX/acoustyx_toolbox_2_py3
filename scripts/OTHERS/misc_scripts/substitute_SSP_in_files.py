# -*- coding: utf-8 -*-
"""
Created on Tue Sep  6 18:32:47 2016

@author: psakicki
"""

from megalib import *

#['2007','2009','2010'],['Median','Pire','Meilleur']


Cpath = "/home/psakicki/THESE/RENDU/1609_SIMU_FINAL/FINALb_good_traingle/VariTempor/basic/FINALb_VariTempor_basic_1000_R50_nois1-0.0_vit1__date20070424_grand-carree_/FINALb_VariTempor_basic_1000_R50_nois1-0.0_vit1__date20070424_grand-carree__.C.dat"
Zpath = "/home/psakicki/THESE/RENDU/1609_SIMU_FINAL/FINALb_good_traingle/VariTempor/basic/FINALb_VariTempor_basic_1000_R50_nois1-0.0_vit1__date20070424_grand-carree_/FINALb_VariTempor_basic_1000_R50_nois1-0.0_vit1__date20070424_grand-carree__.Z.dat"

Cpath = "/home/psakicki/THESE/RENDU/1609_SIMU_FINAL/FINALb_good_traingle/VariTempor/basic/FINALb_VariTempor_basic_1000_R50_nois1-0.0_vit1__date20070424_grand-carree_/FINALb_VariTempor_basic_1000_R50_nois1-0.0_vit1__date20070424_grand-carree__.C.dat"
Zpath = "/home/psakicki/THESE/RENDU/1609_SIMU_FINAL/FINALb_good_traingle/VariTempor/basic/FINALb_VariTempor_basic_1000_R50_nois1-0.0_vit1__date20070424_grand-carree_/FINALb_VariTempor_basic_1000_R50_nois1-0.0_vit1__date20070424_grand-carree__.Z.dat"

Cpath="/home/psakicki/THESE/RENDU/1609_SIMU_FINAL/FINALb_good_traingle/VariTempor/basic/FINALb_VariTempor_basic_1000_R50_nois1-0.0_vit1__date20100225_petit-carree_/FINALb_VariTempor_basic_1000_R50_nois1-0.0_vit1__date20100225_petit-carree__.C.dat"
Zpath="/home/psakicki/THESE/RENDU/1609_SIMU_FINAL/FINALb_good_traingle/VariTempor/basic/FINALb_VariTempor_basic_1000_R50_nois1-0.0_vit1__date20100225_petit-carree_/FINALb_VariTempor_basic_1000_R50_nois1-0.0_vit1__date20100225_petit-carree__.Z.dat"

path = '/home/psakicki/THESE/RENDU/1609_SIMU_FINAL/FINALb_good_traingle/VariTempor/basic_mode'
path = '/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/working'
pathsuffix = '/*basic*'

for expdir in glob.glob(path + pathsuffix ):
    shutil.copy(Zpath , glob.glob(expdir + '/*.Z.dat')[0])
    shutil.copy(Cpath , glob.glob(expdir + '/*.C.dat')[0])





