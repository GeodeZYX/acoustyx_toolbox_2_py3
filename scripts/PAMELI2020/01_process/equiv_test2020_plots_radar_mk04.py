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


pixbluecoord = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/ixblue_coords.pik"
DFixblue = utils.pickle_loader(pixbluecoord)
DFixblue.loc[DFixblue.meth == "A","meth"] = "ACW"
DFixblue.loc[DFixblue.meth == "C","meth"] = "CW"




pikpath = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/05_b_eq_a_plus_Up_from_postproc/05_b.pik"
mode = "simple_range"
filterr = 0

pikpath = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/03_g_weight_range_good_2/03_g.pik"
mode = "simple_range"
filterr = 2



pikpath = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/03_l_win_weight_range_5_no_overlap/03_l.pik"
mode = "win_range"
filterr = 1


pikpath = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/03_g_weight_range_good_2/03_g.pik"
mode = "simple_range"
filterr = 10

pikpath = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/03_l_win_weight_range_5_no_overlap/03_l.pik"
mode = "win_range2"
filterr = 1
pikpath_bis = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/03_g_weight_range_good_2/03_g.pik"
filterr_bis = 10



pikpath = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/05_f_day12_cw_acw_3/05_f.pik"
mode = "boxin"
filterr = 5


pikpath = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/06_a_StaticAbove/06_a.pik"
mode = "simple_range"
filterr = 4


pikpath = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/04_b_StaticCenter/PAM_BST_v402a_m21so_d10_ITRF14_RTKLIB/PAM_.pik"
# pikpath = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/04_b_StaticCenter/PAM_BST_v401a_m21ma_d01_ITRF14_RTKLIB/PAM_.pik"
# pikpath = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/04_b_StaticCenter/PAM_BST_v403a_m22so_d11_ITRF14_RTKLIB/PAM_.pik"
mode = "simple_range"
filterr = 3

DFall = utils.pickle_loader(pikpath)
DFall_orig = DFall.copy()

Xbea_apri_iXBlue_dict = {1: np.array([-7.1673863 , 24.29822317, 40.07]),
                         2: np.array([-7.5009646 ,-26.21878226, 40.66]),
                         3: np.array([22.46625595,  2.22251482, 39.68])}

Xbea_apri_calc_LSQ_best_dict = {1: np.array([  -7.09260362,   24.39713696,   37.80365583]),
                                2: np.array([  -7.36952229,  -26.22316284,   38.19201817]),
                                3: np.array([  22.61862272,    2.203012,     37.69790289])}


np.array([[-7.22298422, 24.29822341, 91.99175988],
          [-7.1673863 , 24.29822317, 91.99175981],
          [-7.44537355, 24.84208254, 91.99176222]])

np.array([[ -7.55656224, -26.28058458,  92.58176811],
          [ -7.5009646 , -26.21878226,  92.58176779],
          [ -7.53802982, -26.23114284,  92.58176788]])

np.array([[22.5033212 ,  2.12363201, 91.60174968],
          [22.54038651,  2.19779409, 91.60174984],
          [22.83690879,  2.22251467, 91.6017509 ]])

import matplotlib.lines as mlines
import matplotlib.pyplot as plt

sym_list1 = ['1','2','3','4']
sym_list2 = ['v','^','<','>']

color_list       = ["orange","green","red","blue"]
color_dark_list  = ["xkcd:dark "  + e for e in color_list]
color_light_list = ["xkcd:light " + e for e in color_list]

color_dark_list  = ["C0","C1","C2","C3"]
color_light_list = ["C4","C5","C6","C7"]

sym_boxin_dict = dict()
sym_boxin_dict["CW"]  = ">"
sym_boxin_dict["ACW"] = "<"
sym_boxin_dict["ALL"] = "o"


def core_plot(DFbea_in,mode_in,figax_in=None,refpoint = None,DFixbluein=None):
        if figax_in:
            fig,ax = figax_in
        else:
            fig,ax = plt.subplots()
            ax.axis("equal")

        MarkerLegend = []

        print("refpoint",refpoint)

        
        if not refpoint:
            if not utils.is_iterable(refpoint):
                refpoint = (refpoint,)
                
        activate_DFixblue1 = False

        for rpt in refpoint:
            print("refpoint2",rpt)
            if rpt == "ixblue":
                Xbea_ref = Xbea_apri_iXBlue_dict[idbea]
                col_refpt = "blue"
                label_refpt = "iXBlue position"
    
            elif rpt == "calc_LSQ_best":
                Xbea_ref = Xbea_apri_calc_LSQ_best_dict[idbea]
                col_refpt = "red"
                label_refpt = "LSQ-determined position"
                
            elif "DFixblue1" == rpt:
                activate_DFixblue1 = True
                continue
            
            ax.scatter(Xbea_ref[1],Xbea_ref[0],s=100,marker="*",c=col_refpt)

            nam_marker = mlines.Line2D([], [], 
                                       color=col_refpt,
                                       marker="*",
                                       linestyle='None',
                                       markersize=3, 
                                       label=label_refpt)   
            MarkerLegend.append(nam_marker)
            
        print("activate_DFixblue1",activate_DFixblue1)

        if activate_DFixblue1:
            for irowixb,rowixb in DFixbluein.iterrows():
                col = "C" + str(rowixb.day-1)
                if rowixb.meth == "CW":
                    sym = "4"
                else:
                    sym = "3"
                    
                ax.scatter(rowixb.E,rowixb.N,s=200,marker=sym,c=col)
    
        
                    
                    
                    
                

    
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
            if not activate_DFixblue1:
                MarkerLegend.append(nam_marker)
        ################ ONLY FOR NAME COLORS ################
            
        
        ################ ONLY FOR PARAMETER SYMBOLS ################   
        IW_DVmap = list(sorted(DFbea_in.W_DV.unique()))
        ICESTmap = list(DFbea_in.C_esti.unique())
            
        para_sym_dict = dict()
        para_col_dict = dict()

    
        for iparam , (cest,wdv) in enumerate(itertools.product(ICESTmap,IW_DVmap)):
            if cest:
                sym_list  = sym_list2
                col_list1 = color_light_list
            else:
                sym_list  = sym_list1
                col_list1 = color_dark_list
                
            sym = sym_list[int(np.mod(iparam,len(sym_list)))]
            para_sym_dict[(cest,wdv)] = sym

            col = col_list1[int(np.mod(iparam,len(sym_list)))]
            para_col_dict[(cest,wdv)] = col
    
            nam_marker = mlines.Line2D([], [], 
                                       color="k",
                                       marker=sym,
                                       linestyle='None',
                                       markersize=2, 
                                       label="Cest" + str(cest) + " Wdv" + str(wdv))   
            
            if not activate_DFixblue1:
                MarkerLegend.append(nam_marker)
            
        ################ ONLY FOR PARAMETER SYMBOLS ################    


        ################ ONLY FOR PARAMETER SYMBOLS ################  
        
        if mode == "boxin":    
            
            marker = mlines.Line2D([], [], 
                                       color="k",
                                       marker=">",
                                       linestyle='None',
                                       markersize=2, 
                                       label="clockwise (CW)")   
            MarkerLegend.append(marker)            

            
            marker = mlines.Line2D([], [], 
                                       color="k",
                                       marker="<",
                                       linestyle='None',
                                       markersize=2, 
                                       label="anticlockwise (ACW)")   
            MarkerLegend.append(marker) 
            
            
            marker = mlines.Line2D([], [], 
                                       color="k",
                                       marker="o",
                                       linestyle='None',
                                       markersize=2, 
                                       label="CW+ACW")   
            MarkerLegend.append(marker)    

            marker = mlines.Line2D([], [], 
                                       color="C0",
                                       marker="s",
                                       linestyle='None',
                                       markersize=2, 
                                       label="Day 1")   
            MarkerLegend.append(marker)            

            marker = mlines.Line2D([], [], 
                                       color="C1",
                                       marker="s",
                                       linestyle='None',
                                       markersize=2, 
                                       label="Day 2")   
            MarkerLegend.append(marker)  
            
            
            nam_marker = mlines.Line2D([], [], 
                                       color="k",
                                       marker='4',
                                       linestyle='None',
                                       markersize=3,
                                       label="iXBlue position CW",
                                       linewidth=4)  
            MarkerLegend.append(nam_marker)
            
            nam_marker = mlines.Line2D([], [], 
                                       color="k",
                                       marker='3',
                                       linestyle='None',
                                       markersize=3,
                                       label="iXBlue position ACW",
                                       linewidth=4)  
            MarkerLegend.append(nam_marker)
                
    
    
    
    
        for irow ,row in DFbea_in.iterrows():
            
            symini = para_sym_dict[(row.C_esti,row.W_DV)]
            colini = nam_col_dict[row.NAME]
    
            if  mode_in == "simple_range"  and     row.C_esti: ## Triangle
                linewidths = 0.5
                sym = symini
                facco = colini
                edgco = colini
                plot_error_bar = True
            elif mode_in == "simple_range" and not row.C_esti: ## Tripods
                linewidths = 2
                sym = symini
                facco = colini
                edgco = colini
                plot_error_bar = True

            elif mode_in == "win_range" and     row.C_esti:
                linewidths = 0.5 
                sym = symini
                facco = "none"
                edgco = colini
                plot_error_bar = False
                
            elif mode_in == "win_range" and not row.C_esti:
                linewidths = 0.5 
                sym = symini
                facco = colini
                edgco = colini
                plot_error_bar = False
                
            elif mode_in == "win_range2_uniq_wins":
                linewidths = 0.5
                edgco = para_col_dict[(row.C_esti,row.W_DV)]
                sym = symini
                if sym in sym_list1:
                    linewidths = 0.5
                    facco = edgco
                else:
                    linewidths = 0.5 
                    facco = "none"
                plot_error_bar = False

            elif mode_in == "win_range2_mean":
                facco = para_col_dict[(row.C_esti,row.W_DV)]
                edgco = para_col_dict[(row.C_esti,row.W_DV)]
                sym = symini
                if sym in sym_list1:
                    linewidths = 2.
                else:
                    linewidths = 0.5 
                plot_error_bar = True

            elif mode_in == "win_range2_lsq":
                linewidths = 0.5 
                facco = para_col_dict[(row.C_esti,row.W_DV)]
                edgco = para_col_dict[(row.C_esti,row.W_DV)]
                sym = "*"
                plot_error_bar = True

            elif mode_in == "boxin":
                linewidths = 0.5 
                facco = "C" + str(int(row.day)-1)
                edgco = facco
                sym = sym_boxin_dict[row.boxin]
                plot_error_bar = True
                
            ax.scatter(row.E,row.N,s=100,marker=sym,facecolors=facco,
                       edgecolors=edgco,
                       linewidths=linewidths)
            
            if plot_error_bar:
                ax.errorbar(row.E,row.N,xerr=row.sE,yerr=row.sN,marker=sym,
                            c=facco,elinewidth=1)  
            
    
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

def core_fiter(DFbea,filterr,idbea):
    
    DFbea = DFbea[DFbea.BEA == idbea]

    ########################### FILTER ZONE ###########################
    ############    
    print(filterr)
    DFixblue_use = None
    
    if filterr == 1:
        DFbeaout   = DFbea[DFbea.NAME == 'RGF93_RTKLIB_GPSonly']
        outname = "radar_plot_RGF93_RTKLIB_GPSonly_BEA"  + str(idbea)
        refpt = "calc_LSQ_best"
        refpt = ("ixblue","calc_LSQ_best")
        
    elif filterr == 10:
        DFbeaout   = DFbea[DFbea.NAME == 'RGF93_RTKLIB_GPSonly']
        outname = "radar_plot_RGF93_RTKLIB_GPSonly_BEA"  + str(idbea)
        refpt = "calc_LSQ_best"
        refpt = ("ixblue",)
        
    elif filterr == 2:
        DFbeaout   = DFbea[(DFbea.C_esti) & (DFbea.W_DV == 0.001)]
        outname = "radar_plot_best_param"
        refpt = "ixblue"
        refpt = "calc_LSQ_best"
        refpt = ("ixblue",)

    elif filterr == 3:
        DFbeaout   = DFbea[(DFbea.BEA == idbea)]
        p = re.compile("m[0-9]{2}..")
        manip = p.search(pikpath).group(0)
        outname = "radar_plot_" + manip + "_BEA"  + str(idbea)
        refpt = "ixblue"
        refpt = ("ixblue","calc_LSQ_best")

    elif filterr == 4:
        DFbeaout = DFbea[(DFbea.BEA == idbea)]
        outname = "radar_plot_" + "_BEA"  + str(idbea)
        refpt = ("ixblue","calc_LSQ_best")

    elif filterr == 5:
        print("A")
        DFixblue_use = DFixblue.copy()
        DFixblue_use = DFixblue_use[DFixblue_use.trait == "M"]
        DFixblue_use = DFixblue_use[DFixblue_use.meth.str.contains("CW")]
        DFixblue_use = DFixblue_use[DFixblue_use.bea == idbea]
        print("C",DFixblue_use)

        print("B")
        
        DFbeaout = DFbea[(DFbea.BEA == idbea)]
        outname = "radar_plot_boxin_ACW_" + "_BEA"  + str(idbea)
        refpt = ("DFixblue1",)

    else:
        DFbeaout = DFbea
        outname = "radar_plot"
        refpt = "ixblue"
        refpt = ("ixblue","calc_LSQ_best")
        
    return DFbeaout,outname,refpt,DFixblue_use
        
    
#####################################################################################################################################################

# pikpath = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/03_j_win_weight_range_3/03_j.pik"
# mode = "win_range" ### BAD



# pikpath = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/04_b/04_b.pik"
# mode = "simple_range" ### BAD
# filterr = 1


for idbea in range(1,4):

    DFbea = DFall.copy()
    DFbea,outname,refpt,DFixblue_use = core_fiter(DFbea,filterr,idbea)

    #################pre processing of the mean ##################
    
    DFgrp = DFbea.groupby(["NAME","C_esti","W_DV"])
        
    DFwin_std = DFgrp.std().reset_index()
    DFwin_std = DFwin_std.drop(["C","sN","sE","sD","sC"],axis=1)
    DFwin_std = DFwin_std.rename({"N":"sN","E":"sE","D":"sD"},axis=1)
    
    DFwin_mean = DFgrp.mean().reset_index()
    DFwin_mean = DFwin_mean[["N","E","D"]]
    
    DFwin = pd.concat([DFwin_std,DFwin_mean],axis=1)
    
    ########################### FILTER ZONE ###########################
    
    if mode == "win_range":
        outname += "_win_range"
        
        DFcore_in = DFbea
        mode_core_in = mode
        figax = core_plot(DFcore_in,mode_in=mode_core_in,refpoint=refpt)
        
        DFcore_in = DFwin
        mode_core_in = "simple_range"
        figax = core_plot(DFcore_in,mode_in=mode_core_in,
                          figax_in=figax,refpoint=refpt) 
        
    elif mode == "win_range2":
        DFcore_in = DFbea
        figax = core_plot(DFcore_in,mode_in="win_range2_uniq_wins",refpoint=refpt)
        DFcore_in = DFwin
        figax = core_plot(DFcore_in,mode_in="win_range2_mean",refpoint=refpt,figax_in=figax)
        
        DFall_bis = utils.pickle_loader(pikpath_bis)
        DFall_bis,_,_ = core_fiter(DFall_bis,10,idbea)
        figax = core_plot(DFall_bis,mode_in="win_range2_lsq",figax_in=figax)
        
    elif mode == "boxin":
        DFcore_in = DFbea
        mode_core_in = mode
        figax = core_plot(DFcore_in,mode_in=mode_core_in,
                          refpoint=refpt,
                          DFixbluein=DFixblue_use)
        
    else:
        DFcore_in = DFbea
        mode_core_in = mode
        figax = core_plot(DFcore_in,mode_in=mode_core_in,refpoint=refpt)

    outtype = (".png",".pdf",".figpik")
    figoutpath = utils.figure_saver(figax[0],
                       os.path.dirname(pikpath),
                       outname,
                       outtype)

    

    #figpik = utils.pickle_loader(figoutpath[-1])
                                

    
        
