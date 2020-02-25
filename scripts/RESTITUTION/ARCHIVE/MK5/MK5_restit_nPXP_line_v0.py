# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 13:04:35 2015

@author: psakicki


v4 : en mode mono Z , on injecte le dZi comme observable
v5 : integration des modifs relatives a Geodesea
"""
from tempfile import mkdtemp
import os.path as path
from megalib import *


multiexp_mode = 0
batch_mode    = 0
force_mode    = 0
purge_mode    = 0

if platform.node() == 'calipso':
    gene_path = "/home/psakicki/THESE/DATA/1506_GEODESEA/WORKING/OPERA_DATA"
    gene_path = '/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/working'
    alternat_SSP_path='/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/exemple/SSP/SSP_NOAA_dep5781_20050827000000' 
elif platform.node() == 'psakicki-MS-16Y1':
    gene_path = '/media/psakicki/D32E-8E7E/CODES/acoustyx_toolbox_2/working'
    alternat_SSP_path='/media/psakicki/D32E-8E7E/CODES/acoustyx_toolbox_2/exemple/SSP/SSP_NOAA_dep5781_20050827000000'

if platform.node() == 'calipso':
    exp = "batch_3x500_x2000_y2000_nois1-1e-06_"
    exp = "tst1_3x50_x2000_y2000_nois1-1e-06__"
    exp = 'repris1_3x500_x2000_y2000_nois1-1e-06__'
    exp = 'repris1_3x50_x2000_y2000_nois1-1e-06__'
    exp = 'batch_3x10_x2000_y2000_nois1-1e-06_'
    exp = 'batc2_3x100_x0_y0_nois1-1e-06_'
    exp = "tst1_3x50_x2000_y2000_nois1-1e-06__"
    exp = 'batc2_3x100_x0_y0_nois1-1e-06_'
    exp = 'batc2_3x100_x2000_y2000_nois1-1e-06_'
    exp = 'batc2b_deriv_20000_R50_noisTrue-1e-06__'
    exp = 'batc2b_deriv_5000_R50_noisTrue-1e-06__'
    exp = 'batc2b_deriv_1500_R50_noisTrue-1e-06__'
    
    exp = 'batc2b_deriv_10000_R50_noisTrue-1e-06__'
    exp = 'batc2_3x1500_x0_y0_nois1-1e-06_'
    exp = "batc2_3x100_x2000_y2000_nois1-1e-06_"
    exp = "batc2_3x1500_x2000_y2000_nois1-1e-06_"

    exp = 'batc2c_deriv_decalZ_10000_R50_noisTrue-1e-06__'
    
    exp = 'batc2d_deriv_megadecalZ_1000_R50_noisTrue-1e-06__'
    exp = "batc2c_deriv_decalZ_1000_R50_noisTrue-1e-06__" 

    exp  = 'gsea_1'
    exp  = 'gsea_1b'
    exp  = 'IUEM_LAPTOP-3Beacontracking'
    
    exp = 'batc2c_deriv_decalZ_nopingnoise_5000_R50_noisFalse-0__'
    exp = 'batc2c_deriv_decalZ_noping_notraj_noise_5000_R50_noisFalse-0__'



elif platform.node() == 'psakicki-MS-16Y1':
    exp = "batch_3x50_x2000_y2000_nois1-1e-06_"
    exp = 'batc2c_deriv_decalZ_nopingnoise_5000_R50_noisFalse-0__'
    exp = 'batc2c_deriv_decalZ_noping_notraj_noise_5000_R50_noisFalse-0__'


if not multiexp_mode:
    exp_lis = [exp]
else:
    # == en manuel
    exp_lis = ["batc2_3x100_x2000_y2000_nois1-1e-06_",
               "batc2_3x500_x2000_y2000_nois1-1e-06_",
               "batc2_3x100_x0_y0_nois1-1e-06_",
               "batc2_3x500_x0_y0_nois1-1e-06_",
               "batc2_3x1500_x0_y0_nois1-1e-06_",
               "batc2_3x1500_x2000_y2000_nois1-1e-06_",
               "batc2_3x1000_x0_y0_nois1-1e-06_",
               "batc2_3x1000_x2000_y2000_nois1-1e-06_"]
               
    exp_lis = ['batc2b_deriv_1500_R50_noisTrue-1e-06__',
               'batc2b_deriv_5000_R50_noisTrue-1e-06__']

    # == avec une wildcard
    expprefix = 'tpo1'        
    expprefix = 'batc2_*x0_y0*'        
    expprefix = 'batc2_*'        
    expprefix = 'batc2c_*'        

    exp_lis = glob.glob(os.path.join(gene_path,expprefix + '*'))
    exp_lis = [os.path.basename(e) for e in exp_lis]
    
    
for exp in exp_lis:
    exp_path = os.path.join(gene_path,exp)
    bigdico  = acls.give_me_the_path(exp_path,exp)#,[1,3,5])

    if purge_mode:
        print("WARN : purge mode !!!")
        print("5sec for cancel")
        time.sleep(8)
        
        for ext in ("/*.log","/*.exp"):
            rmlis = glob.glob(exp_path + ext)
            for f in rmlis:
                print('INFO : remove :' , f)
                os.remove(f)
        
    if batch_mode:
        bool_kw_lis = [ 'with_alternat_SSP' ,
        'with_monoZ' ,
        'with_barycenter' ,
         ]
        booldic_var = genefun.boolean_dict(bool_kw_lis)
        booldic_lis = geok.kwargs_for_jacobian({'with_ssp_bilin' : 0 ,
                                                'with_decimate_SSP' : 0},
                                                booldic_var)
    else:
        booldic_lis = [dict()]  
        
    for booldic in booldic_lis:
        try:
            if batch_mode:
                genefun.eval_a_dict(booldic,globals())
                bool_4_exp = [str(int(eval(e))) for e in bool_kw_lis]
                bool_4_exp_str = ''.join(bool_4_exp)
            else:
                bool_4_exp_str = ''
                
            F = genefun.Tee_frontend(exp_path , exp , bool_4_exp_str )
      
            # ==============================================
            # INITIALISATION DES PARAMETRES
            # ==============================================
            
            Npxp = len(bigdico['N'])
            
            # PARAMETRES ANNEXES
            kriter = 10**7
            kriterold = 10**10
            iiter = 0
            h = 0 
            if platform.node() == 'calipso':
                nbproc   = 3
                iitermax = 4
            elif platform.node() == 'psakicki-MS-16Y1':
                nbproc   = 6
                iitermax = 5
                
            PXPold = np.array([9999999,9999999,9999999] * Npxp) 
            PXPnew_stk = []
            ObsASM_lis = []
            Xbato_lis  = []
            TTT_lis    = []

            PXP_lis = []
            PXPapri_lis = []
            fuv = 1
            cleaning_coef = 3
            expdic = collections.OrderedDict()
            iterdic = acls.get_iterdic_in_expdic(expdic,iiter)
            
            A , N , Ninv , P , B = [],[],[],[],[]
        
            with_V_4_P_reinject = 0
            with_time_window    = 0
   
            # PARAMETRES IMPORTANTS
            if not batch_mode:
                with_BL           = 0
                with_ssp_bilin    = 0
                with_monoZ        = 1
                with_alternat_SSP = 1
                with_barycenter   = 1
                with_decimate_SSP = 1
                
                with_ssp_bilin    = 0
                with_decimate_SSP = 1
                with_alternat_SSP = True
                with_barycenter   = False
                with_monoZ        = True
                with_BL           = False
                
                with_ssp_bilin    = 0
                with_decimate_SSP = 1
                with_alternat_SSP = 1
                with_barycenter   = 1
                with_monoZ        = 1
                with_BL           = 1
                with_dzcst        = 1
                with_zmaster      = 0

                with_ssp_bilin    = 0
                with_decimate_SSP = 0
                with_alternat_SSP = 0
                with_barycenter   = 0
                with_monoZ        = 0
                with_BL           = 0
                with_dzcst        = 0
                with_zmaster      = 0
                
                with_delta_extern = 1 # a laisser vrai !!!
                with_tempo_evol   = 0
                with_time_window  = 0
                with_monoc0       = 1
                with_c0_estim     = 1
                    
            if with_zmaster:
                with_barycenter = 1
    
            if with_BL:
                BL = bigdico['B']['d']
            
            if with_ssp_bilin:
                Z  = bigdico['2lin_Z']['d']
                C  = bigdico['2lin_C']['d'] 
            else:
                Z  = bigdico['Z']['d']
                C  = bigdico['C']['d'] 
            
            if with_alternat_SSP:
                alternat_SSP = np.loadtxt(alternat_SSP_path)
                Z = alternat_SSP[:,0]
                C = alternat_SSP[:,1]

            if 1:
                Zraw , Craw = Z , C
                Z,C = ssp.SSP_extrapolate(Z,C,3000,1)


            if with_decimate_SSP:
                print("decimating SSP")
                Zfull , Cfull = Z , C
                Z , C = ssp.SSP_light(Z,C)
                
            # ========================================================
            # DEFINITION DE LA FONCTION EMPIRIQUE
            # ========================================================


            # 1) estimation du dc en fonction de l'angle
            # ========================================================
            zref = 4000
            cccstk = []
            aaastk = []
            #for a in range(30):
            for a in np.arange(-50,50,1):
                ccc = ssp.SSP_mean(Z,C,a,zref)
                cccstk.append(ccc)
                aaastk.append(a)
            
            Poly = np.polyfit(aaastk, cccstk,12)
            
            fited = []
            for a in aaastk:
                #    fited.append(func(a,*popt))
                fited.append(np.polyval(Poly,a))
                
            #plt.clf()
            #plt.plot(aaastk, cccstk,"*")    
            #plt.plot(aaastk,fited,"-+")
            
            # 2) estimation du dc en fct de la prof
            # ========================================================
            cccstk = []
            zzzstk = []
            for z in np.arange(-200,200,10):
                ccc = ssp.SSP_mean(Z,C,0,zref + z)
                cccstk.append(ccc)
                zzzstk.append(z)
                
            k_dz , _ = geok.linear_regression(zzzstk,cccstk)
            
            # 3) fabrication de la fct proprement dite
            # ========================================================            
            sym_var_tup = sympy.symbols('xbato ybato zbato xpxp ypxp zpxp kc c0 tau t t0')
            xbato , ybato , zbato , xpxp , ypxp , zpxp , kc , c0 , tau , t , t0 = sym_var_tup
            sym_var_lis = list(sym_var_tup)
            
            expr_D    = sympy.sqrt((xbato - xpxp) ** 2 + (ybato - ypxp) ** 2 + (zbato - zpxp) ** 2)  
            expr_D2d  = sympy.sqrt((xbato - xpxp) ** 2 + (ybato - ypxp) ** 2)  
            
            expr_angl = sympy.asin((xpxp - zbato) / expr_D) * (180./np.pi)
            expr_angl = sympy.acos((zpxp - zbato) / expr_D) * (180./np.pi)
            
            #expr_C = c0 + popt[0] * sympy.exp(expr_angl * popt[1]) - k_dz * (zpxp - zref)
            
            #PTDMODIF B
            expr_C   = 1500
            expr_C   = np.polyval(Poly,41.47)
            expr_C   = c0
            Polymod  = np.array(Poly)
            Polymod[-1] = 0
            expr_C   = np.polyval(Poly,expr_angl) + k_dz * (zpxp - zref)
            expr_C   = c0 + np.polyval(Polymod,expr_angl) + k_dz * (zpxp - zref)
            
            #PTDMODIF B
            
            #PTDMODIF D
            expr_2   = expr_D / (kc * (t-t0) + c0)
            expr_1a  = expr_D / (c0)
            expr_1b  = expr_D / expr_C
            #PTDMODIF D
            
            #expr_3  = expr_D_sym / ((c0 + dc) 
            
            
            #PTDMODIF C
            if not with_tempo_evol:
                if with_delta_extern:
                    fctobs = expr_1a
            else:
                fctobs = expr_2
                
            fctobs = expr_1b
            #PTDMODIF C
            
            fctobs_lambda = sympy.lambdify(sym_var_lis,fctobs, "numpy")  
            
            fctobs_diff_list = []
            fctobs_diff_lambda_list = []
            
            for sym_var in sym_var_lis:
                fctobs_diff = sympy.diff(fctobs,sym_var)
                fctobs_diff_lambda = sympy.lambdify(sym_var_lis,fctobs_diff, "numpy") 
                fctobs_diff_list.append(fctobs_diff)  
                fctobs_diff_lambda_list.append(fctobs_diff_lambda)
            
            
            sigma_pxp_apri = 10
            kmeta_pxp_apri = 1 # le K normal c'est l'ID du PXP, le meta c'est pour tous
                               # le Kmeta sera un multiple de 10 par convention
            
            k_z = 42
            R_z = np.random.RandomState(k_z)
            err_z = R_z.randn(1) * sigma_pxp_apri
            i_pxp_master = 0
            
            MNfile = 'M' # M pour les SIMU, N pour GEODESEA
            
            
            for ipxp,Mpxp in bigdico[MNfile].items():
                Mdata  = Mpxp['d']
                Mcomm  = Mpxp['c']
                
                if MNfile == 'N':
                    Xbato = list(Mdata[:,1:4])
                elif MNfile == 'M':
                    Xbato = list(Mdata[:,:3])
                Xbato_lis.append(Xbato)

                if MNfile == 'N':
                    ObsASM_load = Mdata[:,-1]
                elif MNfile == 'M':
                    ObsASM_load = Mdata[:,3]
                ObsASM_lis.append(ObsASM_load)  
                
                TTT_lis.append(Mdata[:,0])
                
                PXP = acls.pxp_string_2_array(Mcomm['pxp_coords'])
                PXP_lis.append(PXP)
                nPXP = len(PXP_lis)
                Npxp = nPXP

                k_pxp_apri = ipxp + kmeta_pxp_apri
                R_pxp_apri = np.random.RandomState(k_pxp_apri)
                if with_monoZ:
                    PXPapri_mono = PXP + np.concatenate((R_pxp_apri.randn(3)[0:2] * sigma_pxp_apri , err_z))
                else:
                    if MNfile == 'N': 
                        PXPapri_mono = PXP
                    elif MNfile == 'M':
                        PXPapri_mono = PXP + R_pxp_apri.randn(3) * sigma_pxp_apri

                PXPapri_lis.append(PXPapri_mono)
                
            if not with_monoc0:
                c0_lis0 = [ssp.SSP_mean(Z,C,10,zpxp) for zpxp in [e[-1] for e  in PXPapri_lis]]
                c0_lis  = list(c0_lis0)
                kc_lis  = [0] * Npxp
            else:
                c0_lis0 = [ ssp.SSP_mean(Z,C,0,PXPapri_lis[-1][-1]) ]    
                c0_lis  = list(c0_lis0) 
    
                
            if with_monoZ:
                for i,pxp in enumerate(PXPapri_lis):
                    if i_pxp_master != i:
                        pxp[2] = PXPapri_lis[i_pxp_master][2] 
                
            PXPZ_lis    = np.array(PXP_lis)[:,-1]
            PXPdZ_lis_true = PXPZ_lis - PXPZ_lis[i_pxp_master] 
            noise4PXPdZ = np.random.RandomState(1789).randn(nPXP) * 0.001
            PXPdZ_lis   = noise4PXPdZ +  PXPdZ_lis_true  # <= par defaut c'est le PXP1 qui est maitre
            
            PXPtrue_arr = np.array(PXP_lis)
            shape_arr   = PXPtrue_arr.shape
            PXPapri0_arr = np.array(PXPapri_lis)
            PXPapri0_vct = genefun.vectorialize(PXPapri0_arr)
            PXPapri = PXPapri0_vct
            
            if with_barycenter:
                PXPbary0 = acls.barycenter_calc(PXPapri0_arr)
                PXPbary = PXPbary0
                dPXPapri0_arr = PXPapri0_arr - PXPbary0
                dPXPapri0_vct = genefun.vectorialize(dPXPapri0_arr)
                dPXPapri  = dPXPapri0_vct
                dPXPapri0 = dPXPapri0_vct
                dPXPapri_lis = list(dPXPapri0_arr)   
            
            if with_zmaster:
                PXPref0    = np.array(PXPbary0)
                PXPref0[2] = PXPapri_lis[i_pxp_master][2] 
                PXPref     = PXPref0
                dPXPapri0_arr = PXPapri0_arr - PXPref0
                dPXPapri0_vct = genefun.vectorialize(dPXPapri0_arr)
                dPXPapri0 = dPXPapri0_vct
                dPXPapri  = dPXPapri0_vct
                dPXPapri_lis = list(dPXPapri0_arr)  
                
            if with_time_window:
                ObsASM_lis_orig = list(ObsASM_lis)
                Xbato_lis_orig  = list(Xbato_lis)
                ObsASM_lis  = []
                Xbato_lis   = []   

                ssss = dt.datetime(2015,6,21,22)
                eeee = dt.datetime(2015,6,22,3)

                ssss = dt.datetime(2015,6,22,3)
                eeee = dt.datetime(2015,6,22,7)
                
                ssssp = geok.dt2posix(ssss)
                eeeep = geok.dt2posix(eeee)
                
                for ttt,xbato in zip(TTT_lis,Xbato_lis_orig):
                    _ , outxbato = geok.time_win_basic(ssssp,eeeep,ttt,xbato)
                    Xbato_lis.append(outxbato)
                for ttt,obsasm in zip(TTT_lis,ObsASM_lis_orig):
                    _ , outobsasm = geok.time_win_basic(ssssp,eeeep,ttt,obsasm)
                    ObsASM_lis.append(outobsasm)
                
            ObsASMgoodbool = []
            for obsasm in ObsASM_lis:
                ObsASMgoodbool.append( np.array(len(obsasm) * [True] ))
               
            A_stk = []
   
            # ===============================
            # BOUCLE D'INVERSION
            # ===============================
            print("START")
            print("===================== DEBUT =====================")
            start = genefun.get_timestamp(0)
            acls.print_n_dicwrit( "debut" , str(start) , iterdic , 1 )
            acls.print_n_dicwrit( "nom" , exp , iterdic , 1)
            acls.print_n_dicwrit( "batch_mode" , bool(batch_mode) , iterdic,1)
            acls.print_n_dicwrit( "force_mode" , bool(force_mode) ,  iterdic,1)         
            acls.print_n_dicwrit( "with_BL" , bool(with_BL) , iterdic,1)
            acls.print_n_dicwrit( "with_ssp_bilin" , bool(with_ssp_bilin) , iterdic,1)
            acls.print_n_dicwrit( "with_monoZ" , bool(with_monoZ), iterdic,1)
            acls.print_n_dicwrit( "with_alternat_SSP" , bool(with_alternat_SSP), iterdic,1)
            acls.print_n_dicwrit( "with_barycenter" , bool(with_barycenter), iterdic,1)
            acls.print_n_dicwrit( "with_decimate_SSP" , bool(with_decimate_SSP), iterdic,1)
            acls.print_n_dicwrit( "with_dzcst" , bool(with_dzcst), iterdic,1)
            acls.print_n_dicwrit( "with_zmaster" , bool(with_zmaster), iterdic,1)
            acls.print_n_dicwrit( "with_time_window" , bool(with_time_window), iterdic,1)
            acls.print_n_dicwrit( "with_V_4_P_reinject" , bool(with_V_4_P_reinject), iterdic,1)
            
            if with_time_window:
                acls.print_n_dicwrit( "time window start" , (ssss), iterdic,1)
                acls.print_n_dicwrit( "time window end  " , (eeee), iterdic,1)

            if batch_mode:
                acls.print_n_dicwrit( "params variables" , bool_kw_lis , iterdic,1)
                acls.print_n_dicwrit( "params var. de l'exp." , booldic , iterdic,1)
                
            print("")
            
            iiter = 1
            
#            while np.linalg.norm(PXPapri - PXPold) > 5 * 10**-5 and iiter <= iitermax:
            for iiter in range(iitermax):
                print("------------ iter No" , iiter , "------------")

                Jacobstk  = []
                Bstk      = []
                c0column  = []
                
            #    for ipxp , (Xbato , ObsASM , pxpapri , Tbato) in enumerate(zip(Xbato_lis ,
            #                                                         ObsASM_lis ,
            #                                                         PXPapri_lis,
            #                                                         Tbato_lis)):
                for ipxp , (Xbato , ObsASM , pxpapri ) in enumerate(zip(Xbato_lis ,
                                                                        ObsASM_lis ,
                                                                        PXPapri_lis)):
            
                    lines_stk = []
                    c_v_stk   = []
                
            #        for xbato , oasm , t in zip(Xbato , ObsASM , Tbato):
                    for xbato_v , oasm  in zip(Xbato , ObsASM):
            
                        # xbato xpxp kc c0 tau t
                        if not with_monoc0:
                            if not with_tempo_evol:
                                argstup = (xbato_v[0],xbato_v[1],xbato_v[2],pxpapri[0],pxpapri[1],
                                           pxpapri[2],0,c0_lis0[ipxp],oasm,0,0)
                            else:
                                argstup = (xbato_v[0],xbato_v[1],xbato_v[2],pxpapri[0],pxpapri[1],
                                           pxpapri[2],kc_lis[ipxp],c0_lis[ipxp],oasm,t,Tbato[0])
                        else:
                            ang = np.arccos((pxpapri[2] - xbato_v[2]) / np.linalg.norm(xbato_v - pxpapri)) * (180./np.pi)
                            c_precalc = ssp.SSP_mean(Z,C,ang,pxpapri[2])
                            
                            # PTDMODIF A  
                            c_v = c_precalc
                            c_v = c0_lis[0]
                            # PTDMODIF A
                            c_v_stk.append(c_v)
                            
                            argstup = (xbato_v[0],xbato_v[1],xbato_v[2],
                                       pxpapri[0],pxpapri[1],pxpapri[2],
                                       0,c_v,oasm,0,0) 
                                       
                        if not with_delta_extern:
                            Bstk.append(fctobs_lambda(*argstup))
                        else:
                            Bstk.append(oasm - fctobs_lambda(*argstup))
                
                        xpxpdiff = fctobs_diff_lambda_list[3](*argstup)
                        ypxpdiff = fctobs_diff_lambda_list[4](*argstup)
                        zpxpdiff = fctobs_diff_lambda_list[5](*argstup)  
            
                        kcdiff   = fctobs_diff_lambda_list[6](*argstup) 
                        c0diff   = fctobs_diff_lambda_list[7](*argstup) 
                        
                        if with_tempo_evol:
                            line = [ xpxpdiff , ypxpdiff , zpxpdiff , c0diff , kcdiff ]
                        else:
                            line = [ xpxpdiff , ypxpdiff , zpxpdiff , c0diff ]
                        if with_monoc0:
                            line = [ xpxpdiff , ypxpdiff , zpxpdiff ]
                            c0column.append( c0diff )
                        
                        lines_stk.append(line)
                    
                    Jacobstk.append(np.vstack(lines_stk))
                    
                A    = scipy.linalg.block_diag(*Jacobstk)
                
                if with_c0_estim:
                    if with_monoc0:
                        A = np.column_stack((A,np.array(c0column)))
                B    = np.array(Bstk)   
                
                Ninv = scipy.linalg.inv((A.T).dot(A))
                dX   = Ninv.dot(A.T).dot(B)
            
                if not with_monoc0:
                    if not with_tempo_evol:    
                        dXreshape = dX.reshape((Npxp,4))
                    else:
                        dXreshape = dX.reshape((Npxp,5))
                    if with_c0_estim:
                        c0_lis    = np.array(c0_lis) + dXreshape[:,3]
            
                else:
                    if not with_c0_estim:
                        dXreshape = dX.reshape((Npxp,3))
                    else:
                        dXreshape = dX[:-1].reshape((Npxp,3))
                        c0_lis    = np.array(c0_lis) + dX[-1]
            
                PXPapri_lis = list(np.array(np.array(PXPapri_lis) + dXreshape[:,0:3]))
                if with_tempo_evol:
                    kc_lis = np.array(kc_lis) + dXreshape[:,4]
            
                print("interm PXP" , PXPapri_lis)
                print("interm barycenter" , geok.barycenter(PXPapri_lis))    
            
                V = B - A.dot(dX)
                fuv = geok.fuv_calc(V,A,1/((10 **-6) **2),1)
                Sigma = np.sqrt((np.diag(Ninv) * fuv))
                
#                plt.figure()
#                plt.hist(V,100)
#                plt.figure()
#                plt.plot(c_v_stk)      
            
                acls.print_n_dicwrit('f.u.v.' , fuv , iterdic , 1)
                acls.print_n_dicwrit('dX' , dX , iterdic )
                acls.print_n_dicwrit('som. dX' , np.sum(dX) , iterdic , 1)
                 
                PXPnew_arr = np.array(PXPapri_lis)
                PXPnew     = list(PXPnew_arr)
               
                acls.print_n_dicwrit( 'nouvelles coords.' , PXPnew_arr , iterdic )
                if not batch_mode:
                    PXPnew_stk.append(PXPnew)
                with np.errstate(invalid='ignore'):
                    sigmas = np.sqrt(np.diag(Ninv) * fuv)
                acls.print_n_dicwrit( 'sigmas' , sigmas , iterdic )
                
                barybrut = acls.barycenter_calc(PXPnew_arr.reshape(shape_arr))            

                acls.print_n_dicwrit( "ecart 3D a la postion vraie en distance" , 
                                np.linalg.norm( PXPnew_arr - PXPtrue_arr,axis=1) , iterdic)
                acls.print_n_dicwrit( "ecart 2D a la postion vraie en distance" , 
                                np.linalg.norm((PXPnew_arr - PXPtrue_arr)[:,:2],axis=1) , iterdic)

                acls.print_n_dicwrit( "ecart a la postion vraie en coordonnees"  , 
                                PXPnew_arr - PXPtrue_arr , iterdic ) 
                
                acls.print_n_dicwrit( "BLs Vraies " , acls.print_BL(PXPtrue_arr) , iterdic)
                
                if with_BL:   
                    acls.print_n_dicwrit( "BLs Observees" , acls.print_BL(ObsBL,0) , iterdic)                

                acls.print_n_dicwrit("BL des coords a prioris('Mod')",            
                                     acls.print_BL(PXPapri0_arr),iterdic)
                    
                acls.print_n_dicwrit("BLs News Ã  l'issue de l'inversion",
                                     acls.print_BL(PXPnew_arr),iterdic)
                
                acls.print_n_dicwrit("Barycentre 'brut' : Sum Xpxp / Npxp",
                                     barybrut , iterdic)

                baryvrai = acls.barycenter_calc(PXPtrue_arr)            
                acls.print_n_dicwrit("ecart 3D bary brut/vrai en distance" , 
                np.linalg.norm(barybrut - baryvrai) , iterdic)

                acls.print_n_dicwrit("ecart 2D bary brut/vrai en distance" ,
                np.linalg.norm(barybrut[0:2] - baryvrai[0:2]) , iterdic ) 
                
                acls.print_n_dicwrit("ecart au bary brut/vrai en coords.",
                barybrut - baryvrai , iterdic ) 
                
                if with_barycenter:
                    acls.print_n_dicwrit("Barycentre estime dans le MC",
                                         PXPbary , iterdic )

                    acls.print_n_dicwrit("Barycentre des dXs dans le MC",                
                    acls.barycenter_calc(dPXPnew.reshape(shape_arr)) , iterdic)
            
                    acls.print_n_dicwrit("ecart 3D bary MC/vrai en distance", 
                    np.linalg.norm(PXPbary - baryvrai) , iterdic)
                    
                    acls.print_n_dicwrit("ecart 2D bary MC/vrai en distance", 
                    np.linalg.norm(PXPbary[0:2] - baryvrai[0:2]) , iterdic)
                
                    acls.print_n_dicwrit("ecart au bary MC/bary vrai en coords.",
                    PXPbary - baryvrai , iterdic )
                
                if with_zmaster:
                    acls.print_n_dicwrit("ecart Ã  la balise de rÃ©ference",
                                         PXPnew_arr[:,2] - PXPnew_arr[i_pxp_master,2],
                                         iterdic)
                    acls.print_n_dicwrit("ecart en input Ã  la balise de rÃ©ference",
                                         PXPdZ_lis,
                                         iterdic)                                         
                    acls.print_n_dicwrit("ecart vrai Ã  la balise de rÃ©ference",
                                         PXPdZ_lis_true,
                                         iterdic)   

                iiter=iiter+1           

            print("===================== FINAL =====================")
            end = genefun.get_timestamp(0)
            acls.print_n_dicwrit("fin"   , str(end) , iterdic , 1 )
            acls.print_n_dicwrit("duree" , str(end - start) , iterdic , 1 )
            acls.print_n_dicwrit("nom" , exp , iterdic , 1)
            acls.print_n_dicwrit("batch_mode" , bool(batch_mode) , iterdic , 1)
            acls.print_n_dicwrit("force_mode" , bool(force_mode) , iterdic , 1)
            acls.print_n_dicwrit("with_BL" , bool(with_BL) , iterdic , 1)
            acls.print_n_dicwrit("with_ssp_bilin" , bool(with_ssp_bilin) , iterdic , 1)
            acls.print_n_dicwrit("with_monoZ" , bool(with_monoZ) , iterdic , 1)
            acls.print_n_dicwrit("with_alternat_SSP" , bool(with_alternat_SSP) , iterdic , 1)
            acls.print_n_dicwrit( "with_barycenter" , bool(with_barycenter), iterdic,1)
            acls.print_n_dicwrit( "with_decimate_SSP" , bool(with_decimate_SSP), iterdic,1)
            acls.print_n_dicwrit( "with_dzcst" , bool(with_dzcst), iterdic,1)
            acls.print_n_dicwrit( "with_zmaster" , bool(with_zmaster), iterdic,1)
            acls.print_n_dicwrit( "with_time_window" , bool(with_time_window), iterdic,1)
            acls.print_n_dicwrit( "with_V_4_P_reinject" , bool(with_V_4_P_reinject), iterdic,1)
            if with_time_window:
                acls.print_n_dicwrit( "time window start" , (ssss), iterdic,1)
                acls.print_n_dicwrit( "time window end  " , (eeee), iterdic,1)
                


            if batch_mode:
                acls.print_n_dicwrit( "params variables" , bool_kw_lis , iterdic,1)
                acls.print_n_dicwrit( "params var. de l'exp." , booldic , iterdic,1)
            acls.print_n_dicwrit("Nb pings     " , len(Xbato) , iterdic , 1)
            acls.print_n_dicwrit("Nb PXPs      " , Npxp , iterdic , 1)
            acls.print_n_dicwrit("Taille Jacob." , A.shape , iterdic , 1)
            acls.print_n_dicwrit("Taille Resid." , V.shape , iterdic , 1)
            acls.print_n_dicwrit("Nb iter.     " , iiter  , iterdic , 1)
            
            F.stop()
            
            expdic[-1] = expdic[max(expdic.keys())]
            
            genefun.save_obj_as_file(expdic , exp_path , exp, suffix=bool_4_exp_str)
            
            print("STOP")

        except Exception as e:
            if force_mode:
                print("ERR : exception ... mais on continue !")     
                print(e)
                F.stop()
                continue
            else:
                print("ERR : exception ... on arrÃªte tout !")
                print(e)
                F.stop()
                raise
    try:      
        out_sum_fil = acls.exp_summary(exp_path)
#        acls.plot_cntr_from_expdics(exp_path,exp_path,exp)
#        acls.plot_dist_from_expdics(exp_path,exp_path,exp)
    except:
        pass


#Vnorma      = np.array(V / P.diagonal())
#Vnormclean  = Vnorma[np.abs(Vnorma) < 0.0001 * np.std(Vnorma)]
#
#plt.clf()
#n,bins,patches = plt.hist(Vnormclean ,1000,normed=1)
#gauss = scipy.stats.norm.pdf(bins,np.mean(Vnormclean),np.std(Vnormclean))
#plt.plot(bins,gauss)
#
#plt.clf()
#n,bins,patches = plt.hist(V ,1000,normed=1)
#gauss = scipy.stats.norm.pdf(bins,np.mean(V),np.std(V))
#plt.plot(bins,gauss)
#
#protopath = os.path.join( exp_path , gf.get_timestamp() + '_' + exp )
#plt.savefig(protopath + '.png')
#np.savetxt(protopath + '.V' , V)
#
#geok.chi2_test_lsq(V,A,P)

if 1:
    plt.clf()
    plt.plot(zip(*Xbato)[0] , zip(*Xbato)[1])
    plt.scatter(PXPnew_arr[:,0] , PXPnew_arr[:,1],s=200)
    for (pxp,(i,j)) in zip(PXPnew_arr,itertools.combinations(list(range(nPXP)),2)):
        print(i,j)
        
    for ipxp,pxp in enumerate(PXPnew_arr):
        if with_barycenter:
            kkk = 3
        else:
            kkk = 0
            
        xe,ye,_,_ = geok.error_ellipse(pxp[0],pxp[1], np.sqrt(Ninv[ kkk + ipxp * 3 , kkk + ipxp * 3]) , 
                                       np.sqrt(Ninv[kkk + ipxp * 3 +1 ,kkk + ipxp*3 +1]) , 
                                       Ninv[kkk + ipxp * 3,kkk + ipxp*3 +1]  , scale= 1000)
        plt.plot(xe,ye)
        plt.axis('equal')
    
    


        
        print("final PXP               " , PXPapri_lis)
        print("final barycenter        " , geok.barycenter(PXPapri_lis))
        print("diff true barycenter    " , geok.barycenter(PXPapri_lis) - geok.barycenter(PXP_lis))
        print("diff true bary 3D (norm)" , np.linalg.norm(geok.barycenter(PXPapri_lis) - geok.barycenter(PXP_lis)))
        print("diff true bary 2D (norm)" , np.linalg.norm(geok.barycenter(PXPapri_lis)[0:2] - geok.barycenter(PXP_lis)[0:2]))
        
        
            
            