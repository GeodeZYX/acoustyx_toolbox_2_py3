# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 15:35:34 2015

@author: psakicki
"""

import acouclass as acls
import raytrace as rt

Zbazik = np.arange(0,5000,1)
Cbazik1 = 1500 + 20 * np.cos(Zbazik  * (0.001))
Cbazik2 = 1500 + 20 * np.sin(Zbazik  * (0.001))

SSPbazik1 = acls.SSP(Zbazik,Cbazik1,e=0)
SSPbazik2 = acls.SSP(Zbazik,Cbazik2,e=10)

zm,cm = rt.munk(6000)

SSPmunk=acls.SSP(zm,cm,e=0)



SSPbazik1.plot()
SSPbazik2.plot()
SSPmunk.plot()

#ssf1 = acls.make_SSF3D_from_SSP(SSPbazik1,xmin,xmax,ymin,ymax,xstep,ystep)
#ssf2 = acls.make_SSF3D_from_SSP(SSPbazik2,xmin,xmax,ymin,ymax,xstep,ystep)
ssf_munk = acls.make_SSF3D_from_SSP(SSPmunk,-2500,2500,-2500,2500)
ssf_munk = acls.make_SSF3D_from_SSP(SSPmunk,-2500,2500,-2500,2500)

ssf_munk.plot()
