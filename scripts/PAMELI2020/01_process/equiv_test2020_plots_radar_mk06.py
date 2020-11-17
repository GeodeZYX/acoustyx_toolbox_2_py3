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



import matplotlib.lines as mlines
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerLine2D, HandlerTuple 

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
mode = "win_range2"
filterr = 1
pikpath_bis = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/03_g_weight_range_good_2/03_g.pik"
filterr_bis = 10


pikpath = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/03_l_win_weight_range_5_no_overlap/03_l.pik"
mode = "win_range2"
filterr = 1


pikpath = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/03_g_weight_range_good_2/03_g.pik"
mode = "simple_range"
filterr = 10

######    07_ new circle modes / lever arms DelphINS  #####################################

#### single circle - slice
pikpath = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/07_b_slice_circle_2/07_b.pik"
mode = "win_range2"
filterr = 1
pikpath_bis = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/07_a_weight_range_1/07_a.pik"
filterr_bis = 10

#### single circle - chunck
pikpath = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/07_c_chunk_circle/07_c.pik"
mode = "win_range2"
filterr = 1
pikpath_bis = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/07_a_weight_range_1/07_a.pik"
filterr_bis = 10

#### single circle
pikpath = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/07_d_unik_circle_2/07_d.pik"
mode = "win_range2"
filterr = 1
pikpath_bis = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/07_a_weight_range_1/07_a.pik"
filterr_bis = 10

#### 5 arbitry group circle (03_g)
pikpath = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/03_l_win_weight_range_5_no_overlap/03_l.pik"
mode = "win_range2"
filterr = 1
pikpath_bis = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/03_g_weight_range_good_2/03_g.pik"
filterr_bis = 10





######    03_  with new circle modes   ########################################

#### single circle - slice
pikpath = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/03_o_slice_circle/03_o.pik"
mode = "win_range2"
filterr = 1
pikpath_bis = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/03_g_weight_range_good_2/03_g.pik"
filterr_bis = 10

#### single circle
pikpath = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/03_m_unik_circle/03_m.pik"
mode = "win_range2"
filterr = 1
pikpath_bis = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/03_g_weight_range_good_2/03_g.pik"
filterr_bis = 10

#### single circle - chunck - THE BEST ONE 
pikpath = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/03_n_chunk_circle/03_n.pik"
mode = "win_range2"
filterr = 1
pikpath_bis = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/03_g_weight_range_good_2/03_g.pik"
filterr_bis = 10


#################### CLOCKWISE ANTI CLOCKWISE
pikpath = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/05_f_day12_cw_acw_3/05_f.pik"
mode = "boxin"
filterr = 5




DFall = utils.pickle_loader(pikpath)
DFall_orig = DFall.copy()


#########################################################################################################

############# STATIC
if 0:
    pikpath1 = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/04_b_StaticCenter/PAM_BST_v402a_m21so_d10_ITRF14_RTKLIB/PAM_.pik"
    pikpath2 = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/04_b_StaticCenter/PAM_BST_v401a_m21ma_d01_ITRF14_RTKLIB/PAM_.pik"
    #pikpath3 = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/04_b_StaticCenter/PAM_BST_v403a_m22so_d11_ITRF14_RTKLIB/PAM_.pik"
    DF1 = utils.pickle_loader(pikpath1)
    DF2 = utils.pickle_loader(pikpath2)
    DF1.NAME = "BaAM"
    DF2.NAME = "BaPM"
    filterr = 3
    
    pikpath = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/06_a_StaticAbove/06_a.pik"
    DF3 = utils.pickle_loader(pikpath)
    DF3.NAME = "ABOV"
    DFall_static = pd.concat([DF1,DF2,DF3])
    filterr = 4
    
    DFall = DFall_static
    
    utils.pickle_saver(DFall_static,"/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/","06_static_n_center")
    mode = "static"
#########################################################################################################

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

#IW_DVmap = list(sorted(DFall.W_DV.unique()))
#ICESTmap = list(DFall.C_esti.unique())

IW_DVmap = [0.001, 0.01, 0.1, np.inf]
ICESTmap = [False,True]

col_static_dict = dict()
col_static_dict["BaAM"] = "C0"
col_static_dict["BaPM"] = "C2"
col_static_dict["ABOV"] = "C5"

col_static_long_short_name_dict = dict()
col_static_long_short_name_dict["BaAM"] = "Static slanted (AM)"
col_static_long_short_name_dict["BaPM"] = "Static slanted (PM)"
col_static_long_short_name_dict["ABOV"] = "Static above"    

col_static_mapper_name_dict = dict()
col_static_mapper_name_dict["B"] = "BaAM"
col_static_mapper_name_dict["B"] = "BaPM"
col_static_mapper_name_dict["S"] = "ABOV"    

######################### BLOC FOR WIN RANGE 2, NOT SIMPLE RANGE #########################'
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

Names = DFall.NAME.unique()
nam_col_dict = dict()

for inam , nam in enumerate(Names):     
    if inam < 3:
        nam_col_dict[nam] = "C" + str(inam)
    else:
        nam_col_dict[nam] = "C" + str(inam+4)
        

        # MarkerLegend = []

        # print("refpoint",refpoint)


def core_ref_pt(figax_in,DFbea_in,DFixbluein,mode_ref_pt):
    
    if mode_ref_pt == "boxin":
        ## IXBLUE  REF PT MEDIAN ##############################################
        DFixblue2 = DFixbluein[(DFixbluein.bea.isin(DFbea_in.BEA) ) & (DFixbluein.trait == "M") & (DFixbluein.meth.str.contains("CW"))]
        #print(DFixblue2)

        for irowixb,rowixb in DFixblue2.iterrows():
            col = "C" + str(rowixb.day-1)
            if rowixb.meth == "CW":
                sym = "4"
            else:
                sym = "3"
            figax_in[1].scatter(float(rowixb.E),float(rowixb.N),s=200,marker=sym,c=col)
        ## IXBLUE REF PT LSQ ##################################################
        DFixblue3 = DFixbluein[(DFixbluein.bea.isin(DFbea_in.BEA) ) & (DFixbluein.trait == "L") & (DFixbluein.meth.str.contains("CW"))]
        #print(DFixblue2)

        for irowixb,rowixb in DFixblue3.iterrows():
            col = "C" + str(rowixb.day-1)
            if rowixb.meth == "CW":
                sym = "+"
            else:
                sym = "x"
            figax_in[1].scatter(float(rowixb.E),float(rowixb.N),s=200,marker=sym,c=col)

    if mode_ref_pt == "static":
        DFixblue2 = DFixbluein[(DFixbluein.bea.isin(DFbea_in.BEA) ) & (DFixbluein.trait == "M") & (DFixbluein.meth.isin(("S","B"))) & (DFixbluein.day == 2) ]
        print(DFixblue2)
        
        for irowixb,rowixb in DFixblue2.iterrows():
            col = col_static_dict[col_static_mapper_name_dict[rowixb.meth]]
            figax_in[1].scatter(float(rowixb.E),float(rowixb.N),s=200,marker="*",c=col)



                
def core_legend(figax_in,DFbea_in):
    
        MarkerLegend = []
        LabelLegend  = []
        with_extern_label = False
        
        lba_onoff = lambda x: "I" if x == True else "0"
    
        ###################################################################################################################
        if mode in ("simple_range","win_range"):
            print("WARN: simple_range uses wrong dict !!!")
            Names = DFbea_in.NAME.unique()        
            for inam , nam in enumerate(Names):     
                nam_marker = mlines.Line2D([], [], 
                                          color=nam_col_dict[nam],
                                          marker='s',
                                          linestyle='None',
                                          markersize=2, 
                                          label=nam)  

                MarkerLegend.append(nam_marker)
    
            for iparam , (cest,wdv) in enumerate(itertools.product(ICESTmap,IW_DVmap)):                    
                sym = para_sym_dict[(cest,wdv)]
                col = col_list1[int(np.mod(iparam,len(sym_list)))]
        
                nam_marker = mlines.Line2D([], [], 
                                            color="k",
                                            marker=sym,
                                            linestyle='None',
                                            markersize=2, 
                                            label="Cest" + str(cest) + " Wdv" + str(wdv))   
                
                MarkerLegend.append(nam_marker)
                
        elif mode == "win_range2":
            for iparam , (cest,wdv) in DFbea[["C_esti","W_DV"]].drop_duplicates().iterrows():                    
                sym = para_sym_dict[(cest,wdv)]
                col = para_col_dict[(cest,wdv)]
                
                if np.isinf(wdv):
                    wdvstr = "dv not used"
                else:
                    if wdv == 0.01:
                        wdvstr = "ςdv loose"
                    else:
                        wdvstr = "ςdv constr."
                
                lab = "C estim. " + lba_onoff(cest) + ", " + str(wdvstr) 
                                
                marker  = mlines.Line2D([], [], color=col,marker=sym,linestyle='None',markersize=2,label=lab)   
                marker2 = mlines.Line2D([], [], color=col,marker="*",linestyle='None',markersize=2,label=lab)   

                MarkerLegend.append((marker,marker2))
                LabelLegend.append(lab)
            marker = mlines.Line2D([], [], color="k",marker="*",linestyle='None',markersize=2,label="complete-period LSQ posi.")   
            MarkerLegend.append(marker)                 
            LabelLegend.append("complete-period LSQ posi.")
            with_extern_label = True
        
        ###################################################################################################################
        elif mode == "boxin":  
            marker = mlines.Line2D([], [], color="k",marker=">",linestyle='None',markersize=2,label="clockwise (CW)")   
            MarkerLegend.append(marker)            
            marker = mlines.Line2D([], [],color="k",marker="<",linestyle='None',markersize=2,label="anticlockwise (ACW)")   
            MarkerLegend.append(marker) 
            marker = mlines.Line2D([], [],color="k",marker="o",linestyle='None',markersize=2,label="CW+ACW")   
            MarkerLegend.append(marker)
            marker = mlines.Line2D([], [],color="C0",marker="s",linestyle='None',markersize=2,label="Day 1")   
            MarkerLegend.append(marker)            
            marker = mlines.Line2D([], [],color="C1",marker="s",linestyle='None',markersize=2,label="Day 2")   
            MarkerLegend.append(marker)  
            nam_marker = mlines.Line2D([], [],color="k",marker='4',linestyle='None',markersize=3,label="iXBlue position CW",linewidth=4)  
            MarkerLegend.append(nam_marker)
            nam_marker = mlines.Line2D([], [],color="k",marker='3',linestyle='None',markersize=3,label="iXBlue position ACW",linewidth=4)  
            MarkerLegend.append(nam_marker)
            nam_marker = mlines.Line2D([], [],color="k",marker='+',linestyle='None',markersize=3,label="iXBlue position CW",linewidth=4)  
            MarkerLegend.append(nam_marker)
            nam_marker = mlines.Line2D([], [],color="k",marker='x',linestyle='None',markersize=3,label="iXBlue position ACW",linewidth=4)  
            MarkerLegend.append(nam_marker)
                        
            
        ###################################################################################################################
        elif mode == "static":  
            for iparam , (cest,wdv) in DFbea[["C_esti","W_DV"]].drop_duplicates().iterrows():                    
                sym = para_sym_dict[(cest,wdv)]
                col = para_col_dict[(cest,wdv)]
                                
                if np.isinf(wdv):
                    wdvstr = "dv not used"
                else:
                    if wdv == 0.01:
                        wdvstr = "ςdv loose"
                    else:
                        wdvstr = "ςdv constr."
                        
                lab = "C estim. " + lba_onoff(cest) + ", " + str(wdvstr) 
                                
                marker  = mlines.Line2D([], [], color="k",marker=sym,linestyle='None',markersize=2,label=lab)  
                MarkerLegend.append(marker)
                
            for staticday in col_static_dict.keys():
                col = col_static_dict[staticday]
                sym = "s"
                lab = col_static_long_short_name_dict[staticday]
                
                marker  = mlines.Line2D([], [], color=col,marker=sym,linestyle='None',markersize=2,label=lab)
                MarkerLegend.append(marker)

            nam_marker = mlines.Line2D([], [],color="k",marker='*',linestyle='None',markersize=3,label="iXBlue position",linewidth=4)  
            MarkerLegend.append(nam_marker)
                
                

                        
            
        ##########################################################
        print("Marker",len(MarkerLegend))
        if with_extern_label:
            figax_in[1].legend(handles=MarkerLegend,labels=LabelLegend,loc=0,prop={'size': 10  },handler_map={tuple: HandlerTuple(ndivide=None)})
        else:
            figax_in[1].legend(handles=MarkerLegend,loc=0,prop={'size': 10  },handler_map={tuple: HandlerTuple(ndivide=None)})
            

            
            
        
                
    

def core_plot(DFbea_in,mode_in,figax_in=None,refpoint = None,DFixbluein=None):
        if figax_in:
            fig,ax = figax_in
        else:
            fig,ax = plt.subplots()
            ax.axis("equal")
    
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
                
            elif mode_in == "static":
                linewidths = 0.5 
                facco = col_static_dict[row.NAME]
                edgco = facco
                sym   = para_sym_dict[(row.C_esti,row.W_DV)]
                
                if sym in sym_list1:
                    linewidths = 2.
                else:
                    linewidths = 0.5 
                    
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
            
        # ax.legend(handles=MarkerLegend,loc=0,prop={'size': 10  })
        
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
        plt.xticks(rotation=45, ha='right')
        fig.set_size_inches(8., 8.)
        #fig.tight_layout()
        
        return fig,ax


#####################################################################################################################################################

def core_fiter(DFbea,filterr,idbea,filter_bea=True,filter_dv=True):

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
        DFixblue_use = DFixblue.copy()
        DFixblue_use = DFixblue_use[DFixblue_use.trait == "M"]
        DFixblue_use = DFixblue_use[DFixblue_use.meth.str.contains("CW")]
        DFixblue_use = DFixblue_use[DFixblue_use.bea == idbea]
        
        DFbeaout = DFbea[(DFbea.BEA == idbea)]
        outname = "radar_plot_boxin_ACW_" + "_BEA"  + str(idbea)
        refpt = ("DFixblue1",)

    else:
        DFbeaout = DFbea
        outname = "radar_plot"
        refpt = "ixblue"
        refpt = ("ixblue","calc_LSQ_best")
        
    if filter_bea:
        DFbeaout = DFbeaout[DFbeaout.BEA == idbea]
    
    ##### FILer of ONE component
    if filter_dv:
        DFbeaout = DFbeaout[DFbeaout.W_DV != 0.1]
        
    return DFbeaout,outname,refpt,DFixblue_use
        
    
#####################################################################################################################################################

extra_size_dict = dict()


FigStk = []

for iidbea,idbea in enumerate(range(1,4)):

    DFbea = DFall.copy()
    DFbea,outname,refpt,DFixblue_use = core_fiter(DFbea,filterr,idbea)

    #################pre processing of the mean ##################
    
    DFgrp = DFbea.groupby(["NAME","C_esti","W_DV"])
        
    DFwin_std = DFgrp.std().reset_index()
    DFwin_std = DFwin_std.drop(["C","sN","sE","sD","sC"],axis=1)
    DFwin_std = DFwin_std.rename({"N":"sN","E":"sE","D":"sD"},axis=1)
    
    print("STD",DFwin_std.to_string())
    
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
        
        core_legend(figax,DFbea_in)
        
    elif mode == "win_range2":
        DFcore_in = DFbea
        figax = core_plot(DFcore_in,mode_in="win_range2_uniq_wins",refpoint=refpt)
        DFcore_in = DFwin
        figax = core_plot(DFcore_in,mode_in="win_range2_mean",refpoint=refpt,figax_in=figax)
        
        DFall_bis = utils.pickle_loader(pikpath_bis)
        DFall_bis,_,_,_ = core_fiter(DFall_bis,10,idbea)
        figax = core_plot(DFall_bis,mode_in="win_range2_lsq",figax_in=figax)
        
        core_legend(figax,DFbea)

        #                       center     dE dN
        # extra_size_dict[1] = ((24.5,-7.1)   ,(.1,.1))
        # extra_size_dict[2] = ((-26.29,-7.38),(.1,.1))
        # extra_size_dict[3] = ((2.23+0.01,22.58)  ,(.125,.125))


        extra_size_dict[1] = ((24.4,-7.095)   ,(.08,.08))
        extra_size_dict[2] = ((-26.29,-7.38),(.08,.08))
        extra_size_dict[3] = ((2.23+0.01,22.58)  ,(.08,.08))

        
    elif mode == "boxin":
        DFcore_in = DFbea
        mode_core_in = mode
        figax = core_plot(DFcore_in,mode_in=mode_core_in,
                          refpoint=refpt,
                          DFixbluein=DFixblue_use)
        
        core_legend(figax,DFcore_in)
        core_ref_pt(figax,DFbea,DFixblue,mode)

    elif mode == "static":
        DFcore_in = DFbea
        mode_core_in = mode
        figax = core_plot(DFcore_in,mode_in=mode_core_in,
                          refpoint=refpt,
                          DFixbluein=DFixblue_use)
        core_legend(figax,DFcore_in)
        core_ref_pt(figax,DFbea,DFixblue,mode)
        
        #                       center     dE dN
        extra_size_dict[1] = ((24.5,-7.38),(.6,.6))
        extra_size_dict[2] = ((-26.2,-7.2),(.6,.6))
        extra_size_dict[3] = ((2.1,22.7),(.6,.6))

    elif mode == "simple_range":
        DFcore_in = DFbea
        mode_core_in = mode
        figax = core_plot(DFcore_in,mode_in=mode_core_in,refpoint=refpt)
        core_ref_pt(figax,DFcore_in,DFixblue,mode)

        
    else:
        DFcore_in = DFbea
        mode_core_in = mode
        figax = core_plot(DFcore_in,mode_in=mode_core_in,refpoint=refpt)
        core_legend(figax,DFcore_in,DFixblue,mode)

    from matplotlib_scalebar.scalebar import ScaleBar
    #from matplotlib_scalebar.scalebar import
    scalebar = ScaleBar(1, 'm',location='lower right',box_alpha=0,length_fraction=0.2) # 1 pixel = 0.2 feet
    

    figax[1].add_artist(scalebar)
    plt.show()

    
    outtype = (".png",".pdf",".figpik")
    figoutpath = utils.figure_saver(figax[0],
                       os.path.dirname(pikpath),
                       outname,
                       outtype)
    


    FigStk.append(figax)
    
    
# for idbea in range(0,3):
    
    print(idbea)

    figax = FigStk[iidbea]

    if len(extra_size_dict) > 0:
        
        #plt.tight_layout()
        axx = figax[1]
        ((xcentr,ycentr),(xstp,ystp)) = extra_size_dict[iidbea+1]
        
        xlimm = (xcentr-xstp,xcentr+xstp)
        ylimm = (ycentr-ystp,ycentr+ystp)
    
        #axx.set_xbound(xlimm)
        #axx.set_ybound(ylimm)    
    
        axx.set_aspect(1)
        axx.autoscale(False)
    
        axx.set_xlim(xlimm)
        axx.set_ylim(ylimm)   
        
        from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)
        
        #axx.xaxis.set_major_locator(MultipleLocator(0.02))
        #axx.yaxis.set_major_locator(MultipleLocator(0.02))
        
        
        # print("****************************************")        
        # print(xlimm)
        # print(axx.get_xlim())
        # print(np.array(xlimm) - np.array(axx.get_xlim()))
        # print("----------------------------------------")        
        # print(ylimm)
        # print(axx.get_ylim())
        # print(np.array(ylimm) - np.array(axx.get_ylim()))
        # print("****************************************")    
       
    
        outtype = (".png",".pdf",".figpik")
    
        figoutpath = utils.figure_saver(figax[0],
                            os.path.dirname(pikpath),
                            outname + "_zoom",
                            outtype)   
                
        print("****************************************")        
        print(xlimm)
        print(axx.get_xlim())
        print(np.array(xlimm) - np.array(axx.get_xlim()))
        print("----------------------------------------")        
        print(ylimm)
        print(axx.get_ylim())
        print(np.array(ylimm) - np.array(axx.get_ylim()))
        print("****************************************")    
       
            
        
        #figpik = utils.pickle_loader(figoutpath[-1])
                                    

    
        
