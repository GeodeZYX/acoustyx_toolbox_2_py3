# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 13:04:35 2015

@author: psakicki

v4 : en mode mono Z , on injecte le dZi comme observable
"""
from tempfile import mkdtemp
import os.path as path
from megalib import *

multiexp_mode = 1
batch_mode    = 0
force_mode    = 0
purge_mode    = 1

if platform.node() == 'calipso':
    gene_path = '/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/working'
    gene_path = "/home/psakicki/THESE/DATA/1506_GEODESEA/WORKING/OPERA_DATA"
    alternat_SSP_path = '/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/exemple/SSP/SSP_NOAA_dep5781_20050827000000'

    gene_path = '/home/psakicki/THESE/DATA/1506_GEODESEA/WORKING/OPERA2//'
    gene_path = '/home/psakicki/THESE/DATA/1506_GEODESEA/WORKING/OPERA3//'

elif platform.node() == 'pierre-MS-16Y1':
    gene_path         = '/home/pierre/Documents/CODES/acoustyx_toolbox_2/working'
    alternat_SSP_path = '/home/pierre/Documents/CODES/acoustyx_toolbox_2/exemple/SSP/SSP_NOAA_dep5781_20050827000000'

elif platform.node() == 'diamant':
    gene_path = '/home/psakicki/GEODESEA_DATA/OPERA1'
    gene_path = '/home/psakicki/GEODESEA_DATA/OPERA2'
    gene_path = '/home/psakicki/GEODESEA_DATA/'


MNfile = 'O' # M pour les SIMU, N pour GEODESEA


if platform.node() == 'calipso':
    exp  = 'IUEM_LAPTOP-3Beacontracking'
else:
    exp  = 'IUEM_LAPTOP-3Beacontracking'


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


    exp_lis = ['OPERA1/IUEM_LAPTOP-3Beacontracking',
               'OPERA2/IUEM_LAPTOP-3Beacontracking']    

    exp_lis = ['IUEM_LAPTOP-3Beacontracking']    
    
for exp in exp_lis:
    expini = exp
    exp_path = os.path.join(gene_path,expini)
    exp = os.path.basename(expini)
    bigdico  = acls.give_me_the_path(exp_path,exp)#,[1,3,5])
    
    if 'OPERA2' in exp_path or 'OPERA3' in exp_path:
        MNfile = 'O'
    else:
        MNfile = 'N'

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


        # bool_kw_lis if you want to apply a Jack Knife
        if 1:
            bool_kw_lis = [ 'with_V_4_P_reinject'  ,
            'with_jackknife' ,
            'with_invert_jk' ]
            
            booldic_var = genefun.boolean_dict(bool_kw_lis)
            booldic_lis = geok.kwargs_for_jacobian({'with_ssp_bilin' : 0 ,
                                                    'with_decimate_SSP' : 0,               
                                                    'with_alternat_SSP' : 0,
                                                    'with_barycenter'   : 0,
                                                    'with_monoZ'        : 0,
                                                    'with_BL'           : 0,
                                                    'with_dzcst'        : 0,
                                                    'with_zmaster'      : 0,
                                                    'with_time_window'    : 0,
                                                    'with_specific_timwin_bool' : 0},
                                                    booldic_var)

        else:
            # bool_kw_lis if you want to apply a Time WIndows
            bool_kw_lis = [ 'with_V_4_P_reinject'  ,
            'with_time_window' ,
            'with_specific_timwin_bool' ]

            booldic_var = genefun.boolean_dict(bool_kw_lis)
            booldic_lis = geok.kwargs_for_jacobian({'with_ssp_bilin'    : 0 ,
                                                    'with_decimate_SSP' : 0,               
                                                    'with_alternat_SSP' : 0,
                                                    'with_barycenter'   : 0,
                                                    'with_monoZ'        : 0,
                                                    'with_BL'           : 0,
                                                    'with_dzcst'        : 0,
                                                    'with_zmaster'      : 0,
                                                    'with_jackknife'    : 0,
                                                    'with_invert_jk'    : 0 },
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
                
            if 'V_4_P_reinject' in list(globals().keys()):
                del V_4_P_reinject
                
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
                iitermax = 1
            elif platform.node() == 'psakicki-MS-16Y1':
                nbproc   = 6
                iitermax = 5
            elif platform.node() == 'diamant':
                nbproc   = 3
                iitermax = 3
                
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
                with_decimate_SSP = 1
                with_alternat_SSP = 0
                with_barycenter   = 0
                with_monoZ        = 0
                with_BL           = 0
                with_dzcst        = 0
                with_zmaster      = 0
                
                with_jackknife    = 0
                with_invert_jk    = 1
                
                with_V_4_P_reinject = 0
                with_time_window    = 0
                with_specific_timwin_bool = 0
                
                with_specific_timwin_bool =  0
                with_ssp_bilin = 0
                with_invert_jk = False
                with_time_window = 0
                with_V_4_P_reinject = 0
                with_decimate_SSP = 0
                with_alternat_SSP = 0
                with_zmaster = 0
                with_barycenter = 0
                with_monoZ = 0
                with_jackknife = 0 
                with_BL = 0
                with_dzcst = 0
   
                    
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
              
            # for jackknife
            keep_ratio        = 0.5
            jk_rand_seed      = 1452
            
            sigma_pxp_apri = 10
            kmeta_pxp_apri = 1 # le K normal c'est l'ID du PXP, le meta c'est pour tous
                               # le Kmeta sera un multiple de 10 par convention
            
            k_z = 42
            R_z = np.random.RandomState(k_z)
            err_z = R_z.randn(1) * sigma_pxp_apri
            i_pxp_master = 0
            

            for ipxp,Mpxp in bigdico[MNfile].items():
                Mdata  = Mpxp['d']
                Mcomm  = Mpxp['c']
                
                if MNfile == 'N':
                    Xbato = list(Mdata[:,1:4])
                elif MNfile == 'M':
                    Xbato = list(Mdata[:,:3])
                elif MNfile == 'O':
                    Xbato_f = list(Mdata[:,1:4])
                    Xbato_b = list(Mdata[:,12:15])
                    Xbato = list(zip(Xbato_f , Xbato_b))
                    
                    
                Xbato_lis.append(Xbato)

                if MNfile == 'N':
                    ObsASM_load = Mdata[:,-1]
                elif MNfile == 'M':
                    ObsASM_load = Mdata[:,3]
                elif MNfile == 'O':
                    ObsASM_load = Mdata[:,-1]

                ObsASM_lis.append(ObsASM_load)  
                if MNfile == 'N':
                    TTT_lis.append(Mdata[:,0])
                elif MNfile == 'O':
                    TTT_lis.append((Mdata[:,0] , Mdata[:,11]))
                    
                PXP = acls.pxp_string_2_array(Mcomm['pxp_coords'])
                PXP_lis.append(PXP)
                nPXP = len(PXP_lis)

                k_pxp_apri = ipxp + kmeta_pxp_apri
                R_pxp_apri = np.random.RandomState(k_pxp_apri)
                if with_monoZ:
                    PXPapri_mono = PXP + np.concatenate((R_pxp_apri.randn(3)[0:2] * sigma_pxp_apri , err_z))
                else:
                    if MNfile in ('N','O'): 
                        PXPapri_mono = PXP
                    elif MNfile == 'M':
                        PXPapri_mono = PXP + R_pxp_apri.randn(3) * sigma_pxp_apri

                PXPapri_lis.append(PXPapri_mono)
                
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

                # Dates en UTC+2
                ssss = dt.datetime(2015,6,22,3)
                eeee = dt.datetime(2015,6,22,7)

                ssss = dt.datetime(2015,6,21,22)
                eeee = dt.datetime(2015,6,22,3)
                
                # Dates en Pacific Time
                ssss = dt.datetime(2015,6,21,16)
                eeee = dt.datetime(2015,6,21,22)

                ssss = dt.datetime(2015,6,21,12)    
                eeee = dt.datetime(2015,6,21,16)
                
                # Dates en Zulu :)
                if with_specific_timwin_bool:
                    ssss = dt.datetime(2015,6,21,20)
                    eeee = dt.datetime(2015,6,22,1) 
                else:            
                    ssss = dt.datetime(2015,6,22,1)
                    eeee = dt.datetime(2015,6,22,5)

               
                
                ssssp = geok.dt2posix(ssss)
                eeeep = geok.dt2posix(eeee)
                
                for ttt,xbato in zip(TTT_lis,Xbato_lis_orig):
                    if MNfile != 'O':
                        _ , outxbato  = geok.time_win_basic(ssssp,eeeep,ttt,xbato)
                    else:
                        _ , outxbato  = geok.time_win_basic(ssssp,eeeep,ttt[-1],xbato)
                        outxbato = [tuple(e) for e in outxbato]
                    Xbato_lis.append(outxbato)
                for ttt,obsasm in zip(TTT_lis,ObsASM_lis_orig):
                    if MNfile != 'O':
                        _ , outobsasm = geok.time_win_basic(ssssp,eeeep,ttt,obsasm)
                    else:
                        _ , outobsasm = geok.time_win_basic(ssssp,eeeep,ttt[-1],obsasm)
                        
                    ObsASM_lis.append(outobsasm)

                
            ObsASMgoodbool = []
            for obsasm in ObsASM_lis:
                if not with_jackknife:
                    ObsASMgoodbool.append( np.array(len(obsasm) * [True] )) 
                else:
                    RState = np.random.RandomState(jk_rand_seed)
                    rand_vals = RState.rand(len(obsasm))
                    bool_arr  = rand_vals < keep_ratio 
                    if with_invert_jk:
                        bool_arr = np.logical_not(bool_arr)
                    ObsASMgoodbool.append( bool_arr )
               
            A_stk = []
   
            # ===============================
            # BOUCLE D'INVERSION
            # ===============================
            print("START")
            print("===================== DEBUT =====================")
            start = genefun.get_timestamp(0)
            acls.print_n_dicwrit( "debut" , str(start) , iterdic , 1 )
            acls.print_n_dicwrit( "nom" , exp , iterdic , 1)
            acls.print_n_dicwrit( "path" , exp_path , iterdic , 1)

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
            acls.print_n_dicwrit( "with_jackknife" , bool(with_jackknife), iterdic,1)

            if with_time_window:
                acls.print_n_dicwrit( "time window start" , (ssss), iterdic,1)
                acls.print_n_dicwrit( "time window end  " , (eeee), iterdic,1)
                

            if batch_mode:
                acls.print_n_dicwrit( "params variables" , bool_kw_lis , iterdic,1)
                acls.print_n_dicwrit( "params var. de l'exp." , booldic , iterdic,1)
                
            if with_jackknife:
                acls.print_n_dicwrit( "jackknife inverted" , with_invert_jk , iterdic,1)
                acls.print_n_dicwrit( "jackknife random seed" , jk_rand_seed , iterdic,1)
                acls.print_n_dicwrit( "keep_ratio" , keep_ratio , iterdic,1)

            print("")
            
            iiter = 1
            
#            while np.linalg.norm(PXPapri - PXPold) > 5 * 10**-5 and iiter <= iitermax:
            while np.linalg.norm(kriter) > 1 * 10**-5 and iiter <= iitermax:
            #while np.linalg.norm(kriter - kriterold) > 10**-4:
                print("------------ iter No" , iiter , "------------")
                iterdic = acls.get_iterdic_in_expdic(expdic,iiter)
            
                # Partie ASM
                if with_dzcst:
                    dz_cst = PXPdZ_lis
                else:
                    dz_cst = None
                    
                ObsASM , ModASM  = acls.vectorialize_ASM_multi(PXPapri_lis,
                                                               ObsASM_lis,
                                                               Xbato_lis,Z,C,
                                                               nbprocs=nbproc,
                                                               dz_cst=dz_cst,
                                                               ObsASMBoolinp_lis=ObsASMgoodbool)
                B_ASM = ObsASM - ModASM
     
                
                if with_zmaster:
                    JacobASM = acls.jacob_ASM((PXPref,dPXPapri_lis), ObsASM_lis,
                                              Xbato_lis,Z,C,h,
                                              nbprocs=nbproc,monoZ=with_monoZ,
                                              accur=1,dz_cst=dz_cst,
                                              ObsASMBoolinp_lis=ObsASMgoodbool)
                elif not with_barycenter:
                    JacobASM = acls.jacob_ASM(PXPapri_lis,ObsASM_lis,Xbato_lis,Z,C,h,
                                                  nbprocs=nbproc,monoZ=with_monoZ,
                                                  accur=1,dz_cst=dz_cst,
                                                  ObsASMBoolinp_lis=ObsASMgoodbool)
                else:
                    JacobASM = acls.jacob_ASM((PXPbary,dPXPapri_lis),ObsASM_lis,Xbato_lis,Z,C,h,
                                                  nbprocs=nbproc,monoZ=with_monoZ,
                                                  accur=1,dz_cst=dz_cst,
                                                  ObsASMBoolinp_lis=ObsASMgoodbool)  
                
                # les baselines sont injectées comme observables => mode contraint et non mode fixé
                #Partie BL
                if with_BL:
                    if with_barycenter:
                        ObsBL,ModBL = acls.vectorialize_BL(BL,dPXPapri_lis,dz_cst)
                    else:
                        ObsBL,ModBL = acls.vectorialize_BL(BL,PXPapri_lis,dz_cst)
                    B_BL = ObsBL - ModBL
                    JacobBL =  acls.jacob_BL(PXPapri_lis,with_monoZ,with_barycenter,dz_cst)
                    
                if with_zmaster:
                    ObsZ = PXPdZ_lis
                    ModZ = np.array([dpxp[-1] for dpxp in dPXPapri_lis])
                    B_Z  = ObsZ - ModZ
                    
                    nz = len(PXPdZ_lis)
                    JacobZ = np.hstack((np.zeros((nz,JacobBL.shape[1] - nz)) , np.diag(np.ones((nz,nz))[:,0])))
                    JacobZ = np.zeros(JacobZ.shape)
                    for i,lin in enumerate(JacobZ):
                        lin[3+ i*3+2] = 1  
                        
                # Partie Commune

                print("fin du remplissage, debut du redesign")
                A = JacobASM
                B = B_ASM
                #P = gf.diagonalize(0.001,np.max(A.shape))
                if not with_V_4_P_reinject or not "V_4_P_reinject" in list(globals().keys()):
                    Ptmp = [1 / ((10**-3)**2)] * np.max(A.shape)
                else:
                    Ptmp = 1 / ((np.array(V_4_P_reinject) * 1) ** 2)
                
                if with_BL:
                    #_,_,Ptmp = geok.weight_mat([10,1],[len(ObsASM),len(ObsBL)],
                    #                        fuv,sparsediag=True)
                    #Ptmp = gf.diagonalize(0.001,len(ObsBL))
                    Ptmp = Ptmp + [0.001] * len(ObsBL)

                    #P = scipy.linalg.block_diag(P,Ptmp)
                    A = np.vstack((A,JacobBL))
                    B = np.hstack((B,B_BL))
                    
                if with_zmaster:
                    #Ptmp = gf.diagonalize(10**0,nz)
                    Ptmp = Ptmp + [10**0] *  nz

                    #P = scipy.linalg.block_diag(P,Ptmp)
                    A = np.vstack((A,JacobZ))
                    B = np.hstack((B,B_Z))
                    
                P = sprs.diags(Ptmp,0)
                
                del Ptmp
                print("fin des poids ")

#                if not batch_mode:        
#                    A_stk.append(A)
                
                # la contrainte de barycentre est injectée suivant la methode fixée d'Helmert
                if with_barycenter:
                    G = acls.constraints_matrix(len(PXP_lis),with_barycenter,with_monoZ)
                    n_lagrangian = 3


                # ==== A partir de là on fat les conv memap / sparse
#                A = np.matrix(A)
#                B = np.matrix(B).T
#                P = np.matrix(P)
                 
                A = genefun.memmap_from_array(A)
                B = genefun.memmap_from_array(B)
                #P = genefun.memmap_from_array(P)
#                N = A.T * P * A
                
                
                print("conversion en sparse ")
                Asprs = sprs.csc_matrix(A)
                Bsprs = sprs.csc_matrix(B).T
                Psprs = sprs.csc_matrix(P)
                print("fin de la conversion en sparse ")
                
                Nsprs = (Asprs.T).dot(Psprs).dot(Asprs)
                N     = Nsprs.toarray()
            
#                N = np.dot(np.dot(A.T , P ) , A)
                
                # Liquidation des matrices inutiles
                JacobASM , ObsASM = [] , []

                # debut contraintes
                if with_barycenter:
                    Gsiz = G.shape[0]
                    O = np.zeros((G.shape[0],G.shape[0]))
                    N = np.vstack((N,G))    
                    N = np.hstack((N,np.hstack((G,O)).T))
                # fin contraintes       
                Ninv = scipy.linalg.inv(N) 
#                AtPB = A.T * P * B
#                AtPB = np.dot(np.dot(A.T , P ) , B)  
                AtPBsprs = (Asprs.T).dot(Psprs).dot(Bsprs)
                AtPB     = AtPBsprs.toarray()
                
                Corel = np.corrcoef(Ninv)
                if with_barycenter:
                    AtPBbis = np.vstack((AtPB,np.matrix(np.zeros(Gsiz)).T))
                    dX = np.dot( Ninv , AtPBbis )
                else:
                    dX = np.dot( Ninv , AtPB )

            #    c,resid,rank,sigma = scipy.linalg.lstsq(A,B)
                # ==== fin de la section des matrix  
                
                dX = np.array(dX).squeeze()
                dXorig = dX
                if with_barycenter:
                    if not with_monoZ:        
                        dX = dX[:-3]
                    else:
                        dX = dX[:-2]
                
                if with_monoZ:
                    dX = np.insert(dX,np.arange(2,len(dX)-2,2),dX[-1]) # manip sioux pour ajouter les dZ communs à tous les PXP 
                if not with_barycenter:
                    PXPold  = PXPapri
                    PXPnew  = PXPapri + dX
                    PXPapri = PXPnew
                    PXPapri_lis = list(PXPnew.reshape(shape_arr))
                else:
                    PXPbaryold = np.array(PXPbary)
                    PXPbarynew = np.array(PXPbary) + dX[:3]
                    PXPbary    = np.array(PXPbarynew) 

                    dPXPnew      = dPXPapri + dX[n_lagrangian:]
                    dPXPapri     = dPXPnew
                    dPXPapri_lis = list(dPXPnew.reshape(shape_arr))    

                    PXPbarynew_tmp    = np.array(PXPbarynew)
                    PXPbarynew_tmp[2] = PXPbaryold[2]
                    PXPnew = dPXPnew + np.tile(PXPbarynew_tmp,Npxp)
                    PXPapri_lis = list(PXPnew.reshape(shape_arr))
                

#                ObsASM_V , ModASM4_V  = acls.vectorialize_ASM_multi(PXPapri_lis,
#                                                               ObsASM_lis,
#                                                               Xbato,Z,C,
#                                                               nbprocs=nbproc,
#                                                               dz_cst=dz_cst)
#                                                               
#                                                               
#                Vnaif = ObsASM_V - ModASM4_V
#                
#                if with_BL:
#                    if with_barycenter:
#                        ObsBL_V,ModBL_V = acls.vectorialize_BL(BL,dPXPapri_lis)
#                    else:
#                        ObsBL_V,ModBL_V = acls.vectorialize_BL(BL,PXPapri_lis)
#                        
#                    V_BL = ObsBL_V - ModBL_V
#                    Vnaif = np.concatenate((Vnaif,V_BL))
                    

                Vrigour = B - A.dot(dXorig[:A.shape[1]])
                
                V = Vrigour

                if cleaning_coef != 0 and not ObsASMgoodbool is None:
                    
                    lenObsASM = [np.sum(o) for o in ObsASMgoodbool]
                    V_per_ObsASM = gf.sublistsIt(V,lenObsASM,True)

                    ObsASMgoodbool_tmp = []
                    
                    if with_V_4_P_reinject:
                        V_4_P_reinject = []
                    for iii  , (boolping  , Voasm) in enumerate(zip(ObsASMgoodbool , V_per_ObsASM )):
                        print("ping valides AVANT nettoyage",np.sum(boolping),'for PXP',iii)
                        indices_true  = np.where(boolping)[0]
                        indices_false = np.where(np.logical_not(boolping))[0]
                                                
                        actualised_bool = boolping[indices_true] * (np.abs(np.array(Voasm)) < cleaning_coef * np.std(Voasm))
                        boolping_new = np.array(boolping)
                        boolping_new[indices_true] = actualised_bool
                            
                        ObsASMgoodbool_tmp.append(boolping_new)
                        if with_V_4_P_reinject:
                            V_4_P_reinject = V_4_P_reinject + list(Voasm[actualised_bool])
                        print("ping valides APRES nettoyage",np.sum(boolping_new),'for PXP',iii)
                       
                    ObsASMgoodbool = ObsASMgoodbool_tmp
                    if with_V_4_P_reinject:
                        V_4_P_reinject = np.array(V_4_P_reinject)
                    
                if not with_BL:
                    fuv = geok.fuv_calc(V,A,P)
                else:
                    fuv = geok.fuv_calc(V,A,P)
                    
                acls.print_n_dicwrit('f.u.v.' , fuv , iterdic , 1)
                acls.print_n_dicwrit('dX' , dX , iterdic )
                acls.print_n_dicwrit('som. dX' , np.sum(dX) , iterdic , 1)
                 
                PXPnew_arr = PXPnew.reshape(shape_arr)
               
                acls.print_n_dicwrit( 'nouvelles coords.' , PXPnew_arr , iterdic )
                if not batch_mode:
                    PXPnew_stk.append(PXPnew)
                with np.errstate(invalid='ignore'):
                    sigmas = np.sqrt(np.diag(Ninv) * fuv)
                acls.print_n_dicwrit( 'sigmas' , sigmas , iterdic )
                
                barybrut = acls.barycenter_calc(PXPnew.reshape(shape_arr))            
                barybrut_old = acls.barycenter_calc(PXPold.reshape(shape_arr))   
                
                kriterold = kriter
                kriter = np.linalg.norm(barybrut_old - barybrut)
                acls.print_n_dicwrit( "kriter/kriterold" , [ kriter , kriterold ] , iterdic)
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
                    
                acls.print_n_dicwrit("BLs News à l'issue de l'inversion",
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
                    acls.print_n_dicwrit("ecart à la balise de réference",
                                         PXPnew_arr[:,2] - PXPnew_arr[i_pxp_master,2],
                                         iterdic)
                    acls.print_n_dicwrit("ecart en input à la balise de réference",
                                         PXPdZ_lis,
                                         iterdic)                                         
                    acls.print_n_dicwrit("ecart vrai à la balise de réference",
                                         PXPdZ_lis_true,
                                         iterdic)   

                iiter=iiter+1           

            print("===================== FINAL =====================")
            end = genefun.get_timestamp(0)
            acls.print_n_dicwrit("fin"   , str(end) , iterdic , 1 )
            acls.print_n_dicwrit("duree" , str(end - start) , iterdic , 1 )
            acls.print_n_dicwrit("nom" , exp , iterdic , 1)
            acls.print_n_dicwrit( "path" , exp_path , iterdic , 1)
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
            acls.print_n_dicwrit( "with_jackknife" , bool(with_jackknife), iterdic,1)

            if with_time_window:
                acls.print_n_dicwrit( "time window start" , (ssss), iterdic,1)
                acls.print_n_dicwrit( "time window end  " , (eeee), iterdic,1)
                


            if batch_mode:
                acls.print_n_dicwrit( "params variables" , bool_kw_lis , iterdic,1)
                acls.print_n_dicwrit( "params var. de l'exp." , booldic , iterdic,1)
                
            if with_jackknife:
                acls.print_n_dicwrit( "jackknife inverted" , with_invert_jk , iterdic,1)
                acls.print_n_dicwrit( "jackknife random seed" , jk_rand_seed , iterdic,1)
                acls.print_n_dicwrit( "keep_ratio" , keep_ratio , iterdic,1)

                
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
                print("ERR : exception ... on arrête tout !")
                print(e)
                F.stop()
                raise
    try:      
        out_sum_fil = acls.exp_summary(exp_path)
#        acls.plot_cntr_from_expdics(exp_path,exp_path,exp)
#        acls.plot_dist_from_expdics(exp_path,exp_path,exp)
    except:
        pass


Vnorma      = np.array(V / P.diagonal())
Vnormclean  = Vnorma[np.abs(Vnorma) < 3 * np.std(Vnorma)]

plt.clf()
n,bins,patches = plt.hist(Vnormclean ,100,normed=1)
gauss = scipy.stats.norm.pdf(bins,np.mean(Vnormclean),np.std(Vnormclean))
plt.plot(bins,gauss)

plt.clf()
n,bins,patches = plt.hist(V ,100,normed=1)
gauss = scipy.stats.norm.pdf(bins,np.mean(V),np.std(V))
plt.plot(bins,gauss)

protopath = os.path.join( exp_path , gf.get_timestamp() + '_' + exp )
plt.savefig(protopath + '.png')
np.savetxt(protopath + '.V' , V)

geok.chi2_test_lsq(V,A,P)

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
    
    

print("MK3_restit 4 gsea")