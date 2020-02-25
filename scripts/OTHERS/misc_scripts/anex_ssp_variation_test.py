# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 14:10:32 2016

@author: psakicki
"""

from megalib import  *
import mayavi as mlab
from mayavi.mlab import *

T = 14.5
P = 820
S = 38.5

SSPlist = []

for dp in np.arange(-20,20,1):
    for ds in np.arange(-1,1,0.1):
        for dt in np.arange(-3,3,1):
            s = S+ds
            t = T+dt
            p = P+dp
            SSP = ssp.soundspeed(s,t,p,'del_grosso')
            SSPlist.append(np.array([s,t,p,SSP]))
            
SSParr = np.vstack(SSPlist)

ssparr = SSParr
plot3d(ssparr[:,0],ssparr[:,1],ssparr[:,2],ssparr[:,3])




            
            
            
            

            





