#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 15:11:38 2020

@author: psakicki
"""

#### Import star style
from geodezyx import *                   # Import the GeodeZYX modules
from geodezyx.externlib import *         # Import the external modules
from geodezyx.megalib.megalib import *   # Import the legacy modules names

p="/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/01_/PAMELI_BREST_vJupyter_2_4transducers/PAMELI_BREST_vJupyter_2_4transducers.B.dat"

M = np.loadtxt(p)

DF = pd.DataFrame()

DF['date'] = np.zeros(np.multiply(*M.shape),dtype=int)
DF['ID_BEA_A'] = np.zeros(np.multiply(*M.shape),dtype=int)
DF['ID_BEA_B'] = np.zeros(np.multiply(*M.shape),dtype=int)
DF['dist'] = np.zeros(np.multiply(*M.shape),dtype=int)

idf = 0
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        DF.loc[idf,'date']     = pd.NaT
        DF.loc[idf,'ID_BEA_A'] = i+1
        DF.loc[idf,'ID_BEA_B'] = j+1
        DF.loc[idf,'dist']     = M[i,j]
        
        idf +=1
        

DF.to_csv("/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/01_/PAMELI_BREST_vJupyter_2_4transducers/PAMELI_BREST_vJupyter_2_4transducers_BASELINE.csv")
    
        