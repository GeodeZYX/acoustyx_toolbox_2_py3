# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 13:04:35 2015

@author: psakicki
"""
from tempfile import mkdtemp
import os.path as path
from megalib import *

multiexp_mode = 0
batch_mode    = 1
force_mode    = 1
purge_mode    = 1

if platform.node() == 'calipso':
    gene_path = '/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/working'
    alternat_SSP_path='/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/exemple/SSP/SSP_NOAA_dep5781_20050827000000'
  
elif platform.node() == 'pierre-MS-16Y1':
    gene_path = '/home/pierre/Documents/CODES/acoustyx_toolbox_2/working'
    alternat_SSP_path='/home/pierre/Documents/CODES/acoustyx_toolbox_2/exemple/SSP/SSP_NOAA_dep5781_20050827000000'

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
    
    exp = "batc2c_deriv_decalZ_1000_R50_noisTrue-1e-06__"  
    exp = 'batc2b_deriv_1500_R50_noisTrue-1e-06__'
    expdir_suffix = 'OLDSTYLE'
    



elif platform.node() == 'pierre-MS-16Y1':
    exp = "batch_3x50_x2000_y2000_nois1-1e-06_"

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

    exp_lis = glob.glob(os.path.join(gene_path,expprefix + '*'))
    exp_lis = [os.path.basename(e) for e in exp_lis]
    
for exp in exp_lis:
    exp_path = os.path.join(gene_path,exp+expdir_suffix)
    bigdico  = acls.give_me_the_path(exp_path,exp) #,[2,3,4])
    
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
        bool_kw_lis = ['with_alternat_SSP' ,
        'with_monoZ',
        'with_barycenter',
        'with_BL']
        booldic_var = genefun.boolean_dict(bool_kw_lis)
        booldic_lis = geok.kwargs_for_jacobian({'with_ssp_bilin' : 0 , 'with_decimate_SSP' : 1},booldic_var)
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
            
            Npxp = len(bigdico['M'])
            
            # PARAMETRES ANNEXES
            kriter = 10**7
            kriterold = 10**10
            iiter = 0
            h = 0 
            if platform.node() == 'calipso':
                nbproc = 4
                iitermax = 4
            elif platform.node() == 'pierre-MS-16Y1':
                nbproc = 6
                iitermax = 4

            PXPold = np.array([9999999,9999999,9999999] * Npxp) 
            PXPnew_stk = []
            ObsASM_lis  = []
            PXP_lis = []
            PXPapri_lis = []
            fuv = 1
            expdic = collections.OrderedDict()
            iterdic = acls.get_iterdic_in_expdic(expdic,iiter)
            
            A , N , Ninv , P , B = [],[],[],[],[]
                 
            # PARAMETRES IMPORTANTS
            if not batch_mode:
                with_BL = 0
                with_ssp_bilin = 0
                with_monoZ = 1
                with_alternat_SSP = 1
                with_barycenter   = 1
                with_decimate_SSP = 1
                
                with_ssp_bilin = 0
                with_decimate_SSP= 1
                with_alternat_SSP= True
                with_barycenter= False
                with_monoZ= True
                with_BL= False

                with_ssp_bilin    = 0
                with_decimate_SSP = 1
                with_alternat_SSP = 1
                with_barycenter   = 1
                with_monoZ        = 1
                with_BL           = 1

                    
            if with_BL:
                BL = bigdico['B']['d']
            
            if with_ssp_bilin:
                Z = bigdico['2lin_Z']['d']
                C = bigdico['2lin_C']['d'] 
            else:
                Z      = bigdico['Z']['d']
                C      = bigdico['C']['d']
            
            if with_alternat_SSP:
                alternat_SSP = np.loadtxt(alternat_SSP_path)
                Z = alternat_SSP[:,0]
                C = alternat_SSP[:,1]
            
            if with_decimate_SSP:
                print("decimating SSP")
                Zfull , Cfull = Z , C
                Z , C = ssp.SSP_light(Z,C)
            
            sigma_pxp_apri = 10
            kmeta_pxp_apri = 1 # le K normal c'est l'ID du PXP, le meta c'est pour tous
                               # le Kmeta sera un multiple de 10 par convention
            
            k_z = 42
            R_z = np.random.RandomState(k_z)
            err_z = R_z.randn(1) * sigma_pxp_apri
                    
            for ipxp,Mpxp in bigdico['M'].items():
                Mdata  = Mpxp['d']
                Mcomm  = Mpxp['c']
                
                Xbato  = list(Mdata[:,:3])
                
                ObsASM = Mdata[:,3]
                ObsASM_lis.append(ObsASM)
                
                PXP = acls.pxp_string_2_array(Mcomm['pxp_coords'])
                PXP_lis.append(PXP)
                k_pxp_apri = ipxp + kmeta_pxp_apri
                R_pxp_apri = np.random.RandomState(k_pxp_apri)
                if with_monoZ:
                    PXPapri_mono = PXP + np.concatenate((R_pxp_apri.randn(3)[0:2] * sigma_pxp_apri , err_z))
                else:
                    PXPapri_mono = PXP + R_pxp_apri.randn(3) * sigma_pxp_apri
                PXPapri_lis.append(PXPapri_mono)
            
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

            if batch_mode:
                acls.print_n_dicwrit( "params variables" , bool_kw_lis , iterdic,1)
                acls.print_n_dicwrit( "params var. de l'exp." , booldic , iterdic,1)
                
            print("")
            
            iiter = 1
#            while np.linalg.norm(PXPapri - PXPold) > 5 * 10**-5 and iiter <= iitermax:
            while np.linalg.norm(kriter) > 1 * 10**-5 and iiter <= iitermax:
            #while np.linalg.norm(kriter - kriterold) > 10**-4:
                print("------------ iter No" , iiter , "------------")
                iterdic = acls.get_iterdic_in_expdic(expdic,iiter)
            
                # Partie ASM
                ObsASM , ModASM  = acls.vectorialize_ASM_multi(PXPapri_lis,ObsASM_lis,
                                                               Xbato,Z,C,nbprocs=nbproc)
                B_ASM = ObsASM - ModASM
                if not with_barycenter:
                    JacobASM = acls.jacob_ASM(PXPapri_lis,ObsASM_lis,Xbato,Z,C,h,
                                                  nbprocs=nbproc,monoZ=with_monoZ,
                                                  accur=1)
                else:
                    JacobASM = acls.jacob_ASM((PXPbary,dPXPapri_lis),ObsASM_lis,Xbato,Z,C,h,
                                                  nbprocs=nbproc,monoZ=with_monoZ,
                                                  accur=1)        
                                   
                #Partie BL
                if with_BL:
                    if with_barycenter:
                        ObsBL,ModBL = acls.vectorialize_BL(BL,dPXPapri_lis)
                    else:
                        ObsBL,ModBL = acls.vectorialize_BL(BL,PXPapri_lis)
                    B_BL = ObsBL - ModBL
                    JacobBL =  acls.jacob_BL(PXPapri_lis,with_monoZ,with_barycenter)
            
                # Partie Commune
                if with_BL:
                    _,_,P = geok.weight_mat([10,1],[len(ObsASM),len(ObsBL)],
                                            fuv,sparsediag=True)
                    A = np.vstack((JacobASM,JacobBL))
                    B = np.hstack((B_ASM,B_BL))      
                else:
                    A = JacobASM
                    B = B_ASM
                    P = sprs.diags(np.ones(len(B)),0)

#                if not batch_mode:        
#                    A_stk.append(A)
                
                # la contrainte 
                if with_barycenter:
                    G = acls.constraints_matrix(len(PXP_lis),with_barycenter,with_monoZ)
                
                # ==== A partir de là on fat les conv memap / sparse
#                A = np.matrix(A)
#                B = np.matrix(B).T
#                P = np.matrix(P)
                 
                A = genefun.memmap_from_array(A)
                B = genefun.memmap_from_array(B)
#                P = genefun.memmap_from_array(P)
#                N = A.T * P * A
                
                Asprs = sprs.csc_matrix(A)
                Bsprs = sprs.csc_matrix(B).T
                Psprs = sprs.csc_matrix(P)
                
                Nsprs = (Asprs.T).dot(Psprs).dot(Asprs)
                N     = Nsprs.toarray()
            
#                N = np.dot(np.dot(A.T , P ) , A)
                
                # LIquidation des matrices inutiles
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

                    dPXPnew      = dPXPapri + dX[3:]
                    dPXPapri     = dPXPnew
                    dPXPapri_lis = list(dPXPnew.reshape(shape_arr))    

                    PXPbarynew_tmp    = np.array(PXPbarynew)
                    PXPbarynew_tmp[2] = PXPbaryold[2]
                    PXPnew = dPXPnew + np.tile(PXPbarynew_tmp,Npxp)
                    PXPapri_lis = list(PXPnew.reshape(shape_arr))
                
                ObsASM_V , ModASM4_V = acls.vectorialize_ASM_multi(PXPapri_lis,ObsASM_lis,
                                                                   Xbato,Z,C,nbprocs=nbproc)
                V = ObsASM_V - ModASM4_V
                
                if with_BL:
                    if with_barycenter:
                        ObsBL_V,ModBL_V = acls.vectorialize_BL(BL,dPXPapri_lis)
                    else:
                        ObsBL_V,ModBL_V = acls.vectorialize_BL(BL,PXPapri_lis)
                        
                    V_BL = ObsBL_V - ModBL_V
                    V = np.concatenate((V,V_BL))
               
                if not with_BL:
                    fuv = geok.fuv_calc(V,A)
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
                    acls.print_n_dicwrit( 'sigmas' , np.sqrt(np.diag(Ninv) * fuv) , iterdic )
                
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