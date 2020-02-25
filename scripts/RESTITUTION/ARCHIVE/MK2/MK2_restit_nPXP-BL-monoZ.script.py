# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 13:04:35 2015

@author: psakicki
"""

from megalib import *

bool_kw_lis = ['with_BL' ,
'with_monoZ' ,
'with_alternat_SSP' ]

booldic_var = genefun.boolean_dict(bool_kw_lis)
booldic_lis = geok.kwargs_for_jacobian({'with_ssp_bilin' : 0},booldic_var)

for a in [1]:
#for booldic in booldic_lis:
#    try:
    if 1:
#        genefun.eval_a_dict(booldic,globals())

        if platform.node() == 'calipso':
            gene_path = '/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/working'
        elif platform.node() == 'pierre-MS-16Y1':
            gene_path = '/home/pierre/Documents/CODES/acoustyx_toolbox_2/working'

        exp = 'tst_3x300_x2000_y2000_nois+True_'
        exp = 'tst_3x300_x2000_y2000_nois+0_'
        exp = 'msiexp3_300_104'
        exp = 'monoZtst_3x100_x2000_y2000_nois+0_'
        exp = 'monoZtst_3x500_x2000_y2000_nois+0_'
        
        exp = "new_3x500_x2000_y2000_nois+0_"
        exp = "new_3x500_x2000_y2000_nois0-2e-06_"
        exp = "new_3x500_x2000_y2000_nois0-0.0002_"
        exp = "new_3x1500_x2000_y2000_nois1-1e-06_"
        exp = "new_3x500_x2000_y2000_nois1-1e-06_"
        exp = "new_3x100_x2000_y2000_nois1-1e-06_"
    
        exp = "batch_3x500_x2000_y2000_nois1-1e-06_"  
        
        exp = "batch_3x50_x2000_y2000_nois1-1e-06_"

    
        exp_path = os.path.join(gene_path,exp)
    
#        bool_4_exp = [str(int(eval(e))) for e in bool_kw_lis]
#        bool_4_exp_str = ''.join(bool_4_exp)
        bool_4_exp_str = ''
        F = genefun.Tee_frontend(exp_path , exp,'monZ' + bool_4_exp_str)
    
        bigdico  = acls.give_me_the_path(exp_path,exp) #,[2,3,4])
        
        #Alternative SSP
        alternat_SSP_path='/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/exemple/SSP/SSP_NOAA_dep5781_20050827000000'
        
        # ==============================================
        # INITIALISATION DES PARAMETRES
        # ==============================================
        
        Npxp = len(bigdico['M'])
        
        # PARAMETRES ANNEXES
        kriter = 10**7
        kriterold = 10**10
        iiter = 1
        iitermax = 2
        h = 0 
        nbproc = 2
        PXPold = np.array([9999999,9999999,9999999] * Npxp) 
        PXPnew_stk = []
        ObsASM_lis  = []
        PXP_lis = []
        PXPapri_lis = []
        fuv = 1
        
        # PARAMETRES IMPORTANTS
        with_BL = 0
        with_ssp_bilin = 0
        with_monoZ = 1
        with_alternat_SSP = 0
                
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
           
        A_stk = []
        
        # ===============================
        # BOUCLE D'INVERSION
        # ===============================
        
        print("************** Résumé (debut) *************")
        print("nom" , exp)
        print("with_BL" , bool(with_BL))
        print("with_ssp_bilin" , bool(with_ssp_bilin))
        print("with_monoZ" , bool(with_monoZ))
        print("with_alternat_SSP" , bool(with_alternat_SSP))
        print("")
        
        while np.linalg.norm(PXPapri - PXPold) > 5 * 10**-5 and iiter < iitermax:
        #while np.linalg.norm(kriter - kriterold) > 10**-4:
            print("============ iter No" , iiter , '============')
        
            # Partie ASM
            ObsASM , ModASM  = acls.vectorialize_ASM_multi(PXPapri_lis,ObsASM_lis,Xbato,Z,C,nbprocs=nbproc)
            B_ASM = ObsASM - ModASM
            JacobASM = acls.jacob_ASM_old(PXPapri_lis,ObsASM_lis,Xbato,Z,C,h,
                                          nbprocs=nbproc,monoZ=with_monoZ,
                                          accur=1)
        
            #Partie BL
            if with_BL:
                ObsBL,ModBL = acls.vectorialize_BL(BL,PXPapri_lis)
                B_BL = ObsBL - ModBL
                JacobBL =  acls.jacob_BL(PXPapri_lis,with_monoZ)
        
            # Partie Commune
            if with_BL:
                K,Q,P = geok.weight_mat([10,1],[len(ObsASM),len(ObsBL)],fuv)
                A = np.vstack((JacobASM,JacobBL))
                B = np.hstack((B_ASM,B_BL))      
            else:
                A = JacobASM
                B = B_ASM
                P = np.eye(len(B))
        
            A_stk.append(A)
            
            # ==== A partir de là on bosse ac des matrix
            A = np.matrix(A)
            B = np.matrix(B).T
            P = np.matrix(P)
        
            N = A.T * P * A
            
            Ninv = scipy.linalg.inv(N) 
            Corel = np.corrcoef(Ninv)
        
            dX = Ninv * A.T * P * B
            
        #    c,resid,rank,sigma = scipy.linalg.lstsq(A,B)
            # ==== fin de la section des matrix  
            
            dX = np.array(dX).squeeze()
            dXorig = dX
            if with_monoZ:
                dX = np.insert(dX,np.arange(2,len(dX)-2,2),dX[-1]) # manip sioux pour ajouter les dZ communs à tous les PXP 
            PXPold  = PXPapri  
            PXPnew  = PXPapri + dX
            
            PXPapri = PXPnew
            PXPapri_lis = list(PXPnew.reshape(shape_arr))
        
            ObsASM_V , ModASM4_V = acls.vectorialize_ASM_multi(PXPapri_lis,ObsASM_lis,Xbato,Z,C,nbprocs=nbproc)
            V = ObsASM_V - ModASM4_V
            if with_BL:
                ObsBL_V,ModBL_V = acls.vectorialize_BL(BL,PXPapri_lis)
                V_BL = ObsBL_V - ModBL_V
                V = np.concatenate((V,V_BL))
           
            if not with_BL:
                fuv = geok.fuv_calc(V,A)
            else:
                fuv = geok.fuv_calc(V,A,P)
            print('f.u.v. :' , fuv)
            
            PXPnew_arr = PXPnew.reshape(shape_arr)
            print('new :',PXPnew_arr)
            PXPnew_stk.append(PXPnew)
            print('sigmas : ' , np.sqrt(np.diag(Ninv) * fuv)) 
            
            kriterold = kriter
            kriter = np.linalg.norm(PXPnew - PXPold)
            print("kriter , kriterold", kriter , kriterold)
            print("ecart à la postion vraie en distance : ") 
            print(np.linalg.norm( PXPnew_arr - PXPtrue_arr,axis=1))
            print("ecart à la postion vraie en coordonnées : ") 
            print(PXPnew_arr - PXPtrue_arr)

            print(" BLs Vraies ")
            acls.print_BL(PXPtrue_arr)
            
            if with_BL:   
                print("BLs Observées ")
                acls.print_BL(ObsBL,0)
            
            print(" BL des coords a prioris('Mod') ")
            acls.print_BL(PXPapri0_arr)
            
            print(" BLs News à l'issue de l'inversion")
            acls.print_BL(PXPnew_arr)
            
            print("Barycentre")
            print(acls.barycenter_calc(PXPnew.reshape(shape_arr)))

            iiter=iiter+1
            
            
        print("************** Résumé (fin) *************")
        print("nom" , exp)
        print("with_BL" , bool(with_BL))
        print("with_ssp_bilin" , bool(with_ssp_bilin))
        print("with_monoZ" , bool(with_monoZ))
        print("with_alternat_SSP" , bool(with_alternat_SSP))
        print("Nb pings     " , len(Xbato))
        print("Nb PXPs      " , len(PXPapri_lis))
        print("Taille Jacob." , A.shape)
        print("Taille Résid." , V.shape)
        
        F.stop()
#    except:
#        continue


    
        
        
    
    