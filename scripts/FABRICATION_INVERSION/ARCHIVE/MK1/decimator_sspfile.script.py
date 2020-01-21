# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 17:19:01 2015

@author: psakicki
"""

import MOVEclass as mcls

filein='/home/psakicki/THESE/bd_CTD/NEW_DATA_2015/MOVE/MOVE3_clean/ssp*'
start = dt.datetime(2010,1,1)
end = dt.datetime(2011,1,1)

for f in glob.glob(filein):
    mcls.decimate_ssp_file(f,start,end)