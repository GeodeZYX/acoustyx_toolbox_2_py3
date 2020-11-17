#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  3 13:15:56 2020

@author: psakicki
"""

#### Import star style
from geodezyx import *                   # Import the GeodeZYX modules
from geodezyx.externlib import *         # Import the external modules
from geodezyx.megalib.megalib import *   # Import the legacy modules names


pikpath = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/05_c_cumulative_window/05_c.pik"

DF = utils.pickle_loader(pikpath)

DFgrp = DF.groupby(["BEA","W_DV"])

bea_col_dict = dict()

bea_col_dict[1] = "blue"
bea_col_dict[2] = "green"
bea_col_dict[3] = "orange"

fig,Ax = plt.subplots(2,2)


for igrp,grp in DFgrp:
    
    grp[["Napri","Eapri","Dapri"]] = grp.iloc[9][["Napri","Eapri","Dapri"]].values
        
    if np.isinf(igrp[1]):
        prefixcol = "light"
        sufixlab = "no DV" 
    else:
        prefixcol = "dark"
        sufixlab = "with DV" 
        
    colnam = "xkcd:" + prefixcol + " " + bea_col_dict[igrp[0]]
    
    #NEDapri = grp[["Napri","Eapri","Dapri"]].values
    NEDapri = grp.iloc[9][["N","E","D"]].values.astype(float)
    NED     = grp[["N","E","D"]].values
    
    dNED = NED - NEDapri
    dNED_dist = np.linalg.norm(dNED[:,:2] ,axis=1)
    
    dNEDplt = np.column_stack((np.abs(dNED),dNED_dist))
    
    print(dNED)
    
    Subtitle = ["North","East","Down","2D Distance"]
    
    for i,(dned,ax) in enumerate(zip(dNEDplt.T,Ax.flatten())):
        lab = "beacon " + str(igrp[0]) + ", " + sufixlab
        ax.plot(grp.win_len_twtt,dned,".-",color=colnam,label=lab)
        ax.set_title(Subtitle[i])
        if i in (2,3):
            ax.set_xlabel("size of the cumulative window")
        if i in (0,2):
            ax.set_ylabel("diff. w.r.t. complete solution (m)")
        plt.legend()

fig.set_size_inches(11.69,8.27)
plt.tight_layout()

