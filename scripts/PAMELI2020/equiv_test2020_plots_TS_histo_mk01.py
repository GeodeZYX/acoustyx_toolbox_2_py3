#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 15:40:05 2020

@author: psakicki
"""

#### Import star style
from geodezyx import *                   # Import the GeodeZYX modules
from geodezyx.externlib import *         # Import the external modules
from geodezyx.megalib.megalib import *   # Import the legacy modules names


plt.ioff()

# p="/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/02_/PAM_BST_v221_m2507-1109_m5_d02_bea2_RGF93_RTKLIB_GPSonly/log/20201021_135322_PAM_BST_v221_m2507-1109_m5_d02_bea2_RGF93_RTKLIB_GPSonly_0.pik"
# p="/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/02_/PAM_BST_v211_m2507-1158_m3_d04_bea1_RGF93_RTKLIB_GPSonly/log/20201021_213501_PAM_BST_v211_m2507-1158_m3_d04_bea1_RGF93_RTKLIB_GPSonly_0.pik"
# p="/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/02_/PAM_BST_v231_m2507-1245_m4_d06_bea3_RGF93_RTKLIB_GPSonly/log/20201021_221405_PAM_BST_v231_m2507-1245_m4_d06_bea3_RGF93_RTKLIB_GPSonly_0.pik"

pbig = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/03_c_weight_range"
pbig = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/03_a_preliminary_runs"
pbig = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/03_l_win_weight_range_5_no_overlap"
pbig = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/05_b_eq_a_plus_Up_from_postproc"
pbig = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/03_g_weight_range_good_2"
pbig = "/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/05_a_reboot03_n_RealTimeUpCorrectedmv"

Pbig = utils.find_recursive(pbig,"*log")
Pbig = list(np.unique([os.path.dirname(e) for e in Pbig]))

Weight = [10**float(n) for n in np.arange(-6,7)]

for plog in Pbig:
    
    P = utils.find_recursive(str(plog),"*.pik")
    
    Lines = []
    
    Xnew_stk = []
    for ip,p in enumerate(P):
        
        #print(p)
        
        plt.close("all")

        # p="/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/02_/PAM_BST_v221_m2507-1109_m5_d02_bea2_RGF93_RTKLIB_GPSonly/log/20201021_135322_PAM_BST_v221_m2507-1109_m5_d02_bea2_RGF93_RTKLIB_GPSonly_0.pik"
        # p="/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/02_/PAM_BST_v211_m2507-1158_m3_d04_bea1_RGF93_RTKLIB_GPSonly/log/20201021_213501_PAM_BST_v211_m2507-1158_m3_d04_bea1_RGF93_RTKLIB_GPSonly_0.pik"
        # p="/home/psakicki/GFZ_WORK/PROJECTS_OTHERS/2001_PAMELi_GNSSA/02_PREPROCESSING/02_/PAM_BST_v231_m2507-1245_m4_d06_bea3_RGF93_RTKLIB_GPSonly/log/20201021_221405_PAM_BST_v231_m2507-1245_m4_d06_bea3_RGF93_RTKLIB_GPSonly_0.pik"

        DictIteraStore = utils.pickle_loader(p)
        DFbig = DictIteraStore[list(DictIteraStore.keys())[-1]]["DFbig"]
        
        if not 'E_TDC_emi' in DFbig.columns:
            ixbluemode = True
        else:
            ixbluemode = False
            
        
        pdir = p[:-4] 
        utils.create_dir(pdir)
        suffix_out = os.path.basename(pdir) + "_"
        outtype = (".png",".pdf",".figpik")
        
        
        # ###################### PLOT FOR TRAJECTORY
        
        if not ixbluemode:

            ############ Trajectory EN
            fig,ax = plt.subplots()
            ax.plot(DFbig["E_AHD_emi"],DFbig["N_AHD_emi"],"+")
            ax.plot(DFbig["E_TDC_rec"],DFbig["N_TDC_rec"],"x")
            DFtdc1 = DFbig[DFbig.ID_TDC == 1]
            ax.plot(DFtdc1["E_TDC_rec"],DFtdc1["N_TDC_rec"],".")
            ax.axis("equal")
            
            ax.set_xlabel("East (m)")
            ax.set_ylabel("North (m)")
            
            plt.tight_layout()
            utils.figure_saver(fig,
                                pdir,
                                suffix_out + "traj_NE",
                                outtype )
            
            
            ############ Trajectory Time Series        
            figt,(axte,axtn) = plt.subplots(2,1)
            figt.suptitle("Trajectory TimeSeries")
            
            axte.plot(DFbig["date_emi"],DFbig["E_AHD_emi"],"+")
            axte.plot(DFbig["date_rec"],DFbig["E_TDC_emi"],"x")

            axtn.plot(DFbig["date_emi"],DFbig["N_AHD_emi"],"+")            
            axtn.plot(DFbig["date_rec"],DFbig["N_TDC_emi"],"x")
            
            
            axte.set_ylabel("East (m)")
            axtn.set_ylabel("North (m)")
            

            plt.tight_layout()
            utils.figure_saver(figt,
                                pdir,
                                suffix_out + "ts_traj_NE",
                                outtype )

            
            figttdc,(axttdce,axttdcn) = plt.subplots(2,1)
            figttdc.suptitle("Trajectory TimeSeries per transducer")
            
            
            for tdc in DFbig.ID_TDC.unique():
                DFtsc = DFbig[DFbig.ID_TDC == tdc]
                tdcstr = "TDC" + str(tdc)
                cols_posi_rec = ['N_' +tdcstr+ '_rec',
                                  'E_' +tdcstr+ '_rec',
                                  'D_' +tdcstr+ '_rec']
                axttdce.plot(DFtsc["date_rec"],DFtsc[cols_posi_rec[1]],"x")
                axttdcn.plot(DFtsc["date_rec"],DFtsc[cols_posi_rec[0]],"x")
                
                axttdce.set_ylabel("East (m)")
                axttdcn.set_ylabel("North (m)")
                
                
            plt.tight_layout()
            utils.figure_saver(figttdc,
                                pdir,
                                suffix_out + "ts_traj_NE2",
                                outtype )
            
            
            if "bea3" in p:
                raise Exception
        
            
            figup,axup = plt.subplots()
            axup.plot(DFbig["date_rec"],DFbig["D_AHD_emi"],"x",c='C0')
            figup.suptitle("Up + TWTT")
            
            axtwtt = axup.twinx()
            axtwtt.plot(DFbig.date_rec[DFbig.VALID],DFbig.TWTT_obs[DFbig.VALID],"+",c='C1')
            
            axup.set_ylabel("Up (m)")
            axtwtt.set_ylabel("TWTT (s)")
            
            plt.tight_layout()
            
            utils.figure_saver(figup,
                                pdir,
                                "ts_TWTT_Up",
                                outtype )
                
        ###################### PLOT FOR TRAJECTORY
        

            # if with_direction_vectors:
            #     fighistangl, axhistangl = plt.subplots()
            #     axhistangl.hist(Angle_residuals,100)
            #     fighistangl.suptitle("Histogram of Direction cosine angle")
        
        
        ###################### PLOT FOR RESIDUALS
        figres,axres   = plt.subplots()
        figres.suptitle("TWTT residuals")
        for idbea in DFbig.ID_BEA.unique():    
            DFbeaplt = DFbig[DFbig.ID_BEA == idbea] 
            axres.plot(DFbeaplt["date_rec"],DFbeaplt["B_TWTT"],"x")
            DFbeaplt_valid = DFbeaplt[DFbeaplt["VALID"]]
            axres.plot(DFbeaplt_valid["date_rec"],DFbeaplt_valid["B_TWTT"],"+")
            axres.set_ylabel("TWTT (s)")
            
        plt.tight_layout()
        utils.figure_saver(figres,
                           pdir,
                           suffix_out + "ts_TWTT_resid",
                           outtype )
        
        B_1stiter       = DictIteraStore[0]['DFbig']["B_TWTT"]
        B_1stiter_bool  = DictIteraStore[0]['DFbig']["VALID"]
        B_1stiter_valid = B_1stiter[B_1stiter_bool] 
        B_lastiter      = DFbeaplt_valid["B_TWTT"] 
        
        nbin = 50
        
        B_1stiter = B_1stiter * 1530 * 100
        B_lastiter = B_lastiter * 1530 * 100
        
        
        # ############# Double iteration iteration
        # fighist,axhist = plt.subplots()
        # fighist.suptitle("Histogram of the TWTT residuals")
        # axhist.hist(B_1stiter_valid,nbin,density=False,label="1st iter",color='C3')
        # xpdf,ypdf = utils.gaussian_for_plot(B_1stiter_valid,False,nbin)
        # axhist.plot(xpdf, ypdf,c='C1')
        
        # axhist.hist(B_lastiter,nbin,density=False,label="last iter",color='C0')
        # xpdf,ypdf = utils.gaussian_for_plot(B_lastiter,False,nbin)
        # axhist.plot(xpdf,ypdf,c='C1')
        
        # axhist.set_xlabel("TWTT (s)")
        # plt.tight_layout()
        
        # utils.figure_saver(fighist,
        #                    pdir,
        #                    suffix_out + "hist_TWTT_resid_1st_lastiter",
        #                    outtype )
        
        
        # ############# Last iteration
        # fighist_last,axhist_last = plt.subplots()
        # fighist_last.suptitle("Histogram of the TWTT residuals")
        # _ = axhist_last.hist(B_lastiter,nbin,density=False,
        #                      label="last iter",color='C0')
        # xpdf,ypdf = utils.gaussian_for_plot(B_lastiter,False,nbin)
        # axhist_last.plot(xpdf,ypdf,c='C1')
        
        # axhist_last.set_xlabel("TWTT (s)")
        # axhist_last.set_xlim((-0.00025,0.00025))
        # axhist_last.set_ylim((0,155))

        # plt.tight_layout()        
        
        # utils.figure_saver(fighist_last,
        #                    pdir,
        #                    suffix_out + "hist_TWTT_resid_lastiter",
        #                    outtype )
        
        ############# Barbapapa plots
        fighist_barba,axhist_barba = plt.subplots()
        fighist_barba.suptitle("Histogram of the TWTT residuals")
        _ = axhist_barba.hist(B_lastiter,
                              nbin,
                              density=False,
                              facecolor="None",
                              label = "full period",
                              edgecolor='C0',
                              histtype='stepfilled')
        
        xpdf,ypdf = utils.gaussian_for_plot(B_lastiter,False,nbin)
        axhist_barba.plot(xpdf,ypdf,c='C0')
        
        axhist_barba.set_xlabel("TWTT (cm)")
        axhist_barba.set_xlim((-40,40))
        axhist_barba.set_ylim((0,155))
        
        nchunk = 5
        B_chunked = utils.chunkIt(B_lastiter,5)
        
        for ichunk,B_chunk_i in enumerate(B_chunked):
            _ = axhist_barba.hist(B_chunk_i,nbin,density=False,
                                  facecolor="None",
                                  edgecolor='C' + str(ichunk +1),
                                  histtype='stepfilled')
            
            xpdf,ypdf = utils.gaussian_for_plot(B_chunk_i,False,nbin)
            axhist_barba.plot(xpdf,ypdf,
                              c='C' + str(ichunk +1),
                              label="period " + str(ichunk*20) + "%-" + str((1+ichunk)*20) + "%")
            
        plt.legend()
        plt.tight_layout()        
        
        utils.figure_saver(fighist_barba,
                           pdir,
                           suffix_out + "hist_TWTT_resid_miniwin",
                           outtype )