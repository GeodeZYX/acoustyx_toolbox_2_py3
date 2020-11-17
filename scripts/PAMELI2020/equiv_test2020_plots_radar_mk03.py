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


def core_plot(DFbea_in,mode_in,figax_in=None):
        if figax_in:
            fig,ax = figax_in
        else:
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
        Names = DFbea_in.NAME.unique()
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
        IW_DVmap = list(sorted(DFbea_in.W_DV.unique()))
        ICESTmap = list(DFbea_in.C_esti.unique())
            
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
    
    
        for irow ,row in DFbea_in.iterrows():
            
            sym = para_sym_dict[(row.C_esti,row.W_DV)]
            col = nam_col_dict[row.NAME]
    
            if  mode_in == "simple_range"  and     row.C_esti: ## Triangle
                linewidths = 0.5
                facco = col
                edgco = col
            elif mode_in == "simple_range" and not row.C_esti: ## Tri
                linewidths = 2
                facco = col
                edgco = col
            elif mode_in == "win_range" and     row.C_esti:
                linewidths = 0.5 
                facco = "none"
                edgco = col
            elif mode_in == "win_range" and not row.C_esti:
                linewidths = 0.5 
                facco = col
                edgco = col       
                
            ax.scatter(row.E,row.N,s=100,marker=sym,facecolors=facco,
                       edgecolors=edgco,
                       linewidths=linewidths)
            if not mode_in == "win_range":
                ax.errorbar(row.E,row.N,xerr=row.sE,yerr=row.sN,marker=sym,
                            c=nam_col_dict[row.NAME],elinewidth=1)        
    
        Xbea_ixblue = Xbea_apri_all_init_dict[idbea]
        ax.scatter(Xbea_ixblue[1],Xbea_ixblue[0],s=100,marker="*",c="red")
    
        # if idbea == 2:
        #     loc = 4
        # else:
        #     loc = 0
            
        ax.legend(handles=MarkerLegend,loc=0,prop={'size': 10  })
        
        xlim = np.round(np.mean(ax.get_xlim()),2)
        ylim = np.round(np.mean(ax.get_ylim()),2)
        
        siz = .12
        ax.set_xlim((xlim - siz,xlim + siz))
        ax.set_ylim((ylim - siz,ylim + siz))
        
        ax.set_xlabel("East (m)")
        ax.set_ylabel("North (m)")
        ax.set_title("Beacon " + str(idbea) )
        #ax.set_aspect('equal','box')
        ax.axis('equal')
        fig.set_size_inches(8., 8.)
        fig.tight_layout()
        
        return fig,ax


#####################################################################################################################################################
#####################################################################################################################################################

# pikpath = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/03_j_win_weight_range_3/03_j.pik"
# mode = "win_range" ### BAD



# pikpath = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/04_b/04_b.pik"
# mode = "simple_range" ### BAD
# filterr = 1


pikpath = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/04_b/PAM_BST_v403a_m22so_d11_ITRF14_RTKLIB/PAM_.pik"
pikpath = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/04_b/PAM_BST_v401a_m21ma_d01_ITRF14_RTKLIB/PAM_.pik"
pikpath = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/04_b/PAM_BST_v402a_m21so_d10_ITRF14_RTKLIB/PAM_.pik"
mode = "simple_range"
filterr = 3 ####### DO NOT FORGET TO CHANGE THE NAME MANUALLY


######### Window mode
pikpath = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/03_l_win_weight_range_5_no_overlap/03_l.pik"
mode = "win_range"
filterr = 1

pikpath = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/03_g_weight_range_good_2/03_g.pik"
mode = "simple_range"
filterr = 2

pikpath = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/03_g_weight_range_good_2/03_g.pik"
mode = "simple_range"
filterr = 1


DFall = utils.pickle_loader(pikpath)

Xbea_apri_all_init_dict = {1: np.array([-7.1673863 , 24.29822317, 40.07]),
                           2: np.array([-7.5009646 ,-26.21878226, 40.66]),
                           3: np.array([22.46625595,  2.22251482, 39.68])}

import matplotlib.lines as mlines
import matplotlib.pyplot as plt

sym_list1 = ['1','2','3','4']
sym_list2 = ['v','^','<','>']

for idbea in range(1,4):

    DFbea_orig = DFall[DFall.BEA == idbea]
    DFbea = DFbea_orig.copy()
    
    
    ########################### FILTER ZONE ###########################
    ############    
    if filterr == 1:
        DFbea   = DFbea[DFbea.NAME == 'RGF93_RTKLIB']
        outname = "radar_plot_RGF93_RTKLIB_BEA"  + str(idbea)
    ############
    elif filterr == 2:
        DFbea   = DFbea[(DFbea.C_esti) & (DFbea.W_DV == 0.001)]
        outname = "radar_plot_best_param"
    elif filterr == 3:
        DFbea   = DFbea[(DFbea.BEA == idbea)]
        outname = "radar_plot_manips2soir_BEA"  + str(idbea)
    else:
        outname = "radar_plot"

        
    ########################### FILTER ZONE ###########################
    
    if mode == "win_range":
        outname += "_win_range"
        DFgrp = DFbea.groupby(["NAME","C_esti","W_DV"])
        
        DFwin_std = DFgrp.std().reset_index()
        DFwin_std = DFwin_std.drop(["C","sN","sE","sD","sC"],axis=1)
        DFwin_std = DFwin_std.rename({"N":"sN","E":"sE","D":"sD"},axis=1)
        
        DFwin_mean = DFgrp.mean().reset_index()
        DFwin_mean = DFwin_mean[["N","E","D"]]
        
        DFwin = pd.concat([DFwin_std,DFwin_mean],axis=1)
        
        DFcore_in = DFbea
        mode_core_in = mode
        figax = core_plot(DFcore_in,mode_in=mode_core_in)
        
        DFcore_in = DFwin
        mode_core_in = "simple_range"
        figax = core_plot(DFcore_in,mode_in=mode_core_in,figax_in=figax) 
    else:
        DFcore_in = DFbea
        mode_core_in = mode
        figax = core_plot(DFcore_in,mode_in=mode_core_in)

    outtype = (".png",".pdf",".figpik")
    figoutpath = utils.figure_saver(figax[0],
                       os.path.dirname(pikpath),
                       outname,
                       outtype)
    
    figpik = utils.pickle_loader(figoutpath[-1])
                                

    
        
