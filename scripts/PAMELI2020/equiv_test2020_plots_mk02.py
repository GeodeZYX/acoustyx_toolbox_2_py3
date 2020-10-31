#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 27 15:18:56 2020

@author: psakicki
"""

#### Import star style
from geodezyx import *                   # Import the GeodeZYX modules
from geodezyx.externlib import *         # Import the external modules
from geodezyx.megalib.megalib import *   # Import the legacy modules names



pikpath = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/03_g_weight_range_good_2/03_g.pik"
mode = "simple_range"

pikpath = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/03_j_win_weight_range_3/03_j.pik"
mode = "win_range"






DFall = utils.pickle_loader(pikpath)

Xbea_apri_all_init_dict = {1: np.array([-7.1673863 , 24.29822317, 40.07]),
                           2: np.array([-7.5009646 ,-26.21878226, 40.66]),
                           3: np.array([22.46625595, 2.22251482 , 39.68])}

import matplotlib.lines as mlines
import matplotlib.pyplot as plt

sym_list1 = ['1','2','3','4']
sym_list2 = ['v','^','<','>']

for idbea in range(1,4):

    DFbea = DFall[DFall.BEA == idbea]    
    DFgrp = DFbea.groupby(["NAME","C_esti","W_DV"])
    

    
    fig,ax = plt.subplots()
    ax.axis("equal")

    MarkerLegend = []
    nam_marker = mlines.Line2D([], [], 
                               color="red",
                               marker="*",
                               linestyle='None',
                               markersize=2, 
                               label="iXBlue position")   
    MarkerLegend.append(nam_marker)
    

    ################ ONLY FOR NAME COLORS ################
    Names = DFbea.NAME.unique()
    nam_col_dict = dict()

    for inam , nam in enumerate(Names):
        if inam < 3:
            nam_col_dict[nam] = "C" + str(inam)
        else:
            nam_col_dict[nam] = "C" + str(inam+4)
            
        nam_marker = mlines.Line2D([], [], 
                                  color=nam_col_dict[nam],
                                  marker='s',
                                  linestyle='None',
                                  markersize=2, 
                                  label=nam)   
        MarkerLegend.append(nam_marker)
    ################ ONLY FOR NAME COLORS ################
        
    
    ################ ONLY FOR PARAMETER SYMBOLS ################    
    IW_DVmap = list(DFbea.W_DV.unique())
    ICESTmap = list(DFbea.C_esti.unique())


    para_sym_dict = dict()

    for iparam , (cest,wdv) in enumerate(itertools.product(ICESTmap,IW_DVmap)):
        
        
        if cest:
            sym_list = sym_list2
        else:
            sym_list = sym_list1
            
        sym = sym_list[int(np.mod(iparam,len(sym_list)))]
        para_sym_dict[(cest,wdv)] = sym

        nam_marker = mlines.Line2D([], [], 
                                   color="k",
                                   marker=sym,
                                   linestyle='None',
                                   markersize=2, 
                                   label="Cest" + str(cest) + " Wdv" + str(wdv))   
        MarkerLegend.append(nam_marker)
    ################ ONLY FOR PARAMETER SYMBOLS ################    


    for irow ,row in DFbea.iterrows():
        
        sym = para_sym_dict[(row.C_esti,row.W_DV)]
        col = nam_col_dict[row.NAME]

        if  mode == "simple_range"  and     row.C_esti: ## Triangle
            linewidths = 0.5
            facco = col
            edgco = col
        elif mode == "simple_range" and not row.C_esti: ## Tri
            linewidths = 2
            facco = col
            edgco = col
        elif mode == "win_range" and     row.C_esti:
            linewidths = 0.5 
            facco = "none"
            edgco = col
        elif mode == "win_range" and not row.C_esti:
            linewidths = 0.5 
            facco = col
            edgco = col       
            
        ax.scatter(row.E,row.N,s=100,marker=sym,facecolors=facco,
                   edgecolors=edgco,
                   linewidths=linewidths)
        if not mode == "win_range":
            ax.errorbar(row.E,row.N,yerr=row.sE,xerr=row.sN,marker=sym,
                        c=nam_col_dict[row.NAME],elinewidth=1)        

    Xbea_ixblue = Xbea_apri_all_init_dict[idbea]
    ax.scatter(Xbea_ixblue[1],Xbea_ixblue[0],s=100,marker="*",c="red")

    if idbea == 2:
        loc = 4
    else:
        loc = 0
        
    ax.legend(handles=MarkerLegend,loc=loc)
    #ax.set_aspect('equal','box')
    ax.set_xlabel("East (m)")
    ax.set_ylabel("North (m)")
    ax.set_title("Beacon " + str(idbea) )
    ax.axis('equal')
    fig.set_size_inches(8., 8.)
    fig.tight_layout()
        
