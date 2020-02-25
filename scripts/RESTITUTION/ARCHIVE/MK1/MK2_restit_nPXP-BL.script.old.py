# -*- coding: utf-8 -*-
"""
Created on Fri Jul 31 11:12:30 2015

@author: psakicki
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 13:04:35 2015

@author: psakicki
"""

from megalib import *

gene_path = "/home/psakicki/THESE/DATA/1506_GEODESEA/WORKING/OPERA_DATA"
gene_path = '/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox/working'
gene_path = '/home/pierre/Documents/CODES/acoustyx_toolbox//working'
gene_path = '/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox/working'


exp  = 'gsea_1'
exp  = 'compar_SSPorig_SSPbilin_nonoise'
exp  = 'compar_SSPorig_SSPbilin_noised'
exp  = 'test6'
exp  = 'test8_1500'


exp = "msiexp1"
exp = 'msiexp2_300'

exp = 'test_impred_std'
exp  = 'test3'
exp = 'msiexp3_300_104'


exp_path = os.path.join(gene_path,exp)
bigdico  = acls.give_me_the_path(exp_path,exp)#,[2,3,4])

#F = genefun.Tee_frontend(exp_path , exp)

# ==============================================
# INITIALISATION DES PARAMETRES
# ==============================================

Npxp = len(bigdico['M'])

# PARAMETRES ANNEXES
kriter = 10**7
kriterold = 10**10
iiter = 1
iitermax = 5
h = 0 
nbproc = 4
PXPold = np.array([9999999,9999999,9999999] * Npxp) 
PXPnew_stk = []
ObsASM_lis  = []
PXP_lis = []
PXPapri_lis = []

# PARAMETRES IMPORTANTS
with_BL = 1
with_ssp_bilin = 0

fuv = 1

if with_BL:
    BL = bigdico['B']['d']

if with_ssp_bilin:
    Z = bigdico['2lin_Z']['d']
    C = bigdico['2lin_C']['d'] 
else:
    Z      = bigdico['Z']['d']
    C      = bigdico['C']['d']

sigma_pxp_apri = 10
kmeta_pxp_apri = 1 # le K normal c'est l'ID du PXP, le meta c'est pour tous
                   # le Kmeta sera un multiple de 10 par convention
for ipxp,Mpxp in bigdico['M'].items():
    Mdata  = Mpxp['d']
    Mcomm  = Mpxp['c']
    
    Xbato  = list(Mdata[:,:3])
    
    ObsASM = Mdata[:,3]
    ObsASM_lis.append(ObsASM)
    
    PXP = PXP    = acls.pxp_string_2_array(Mcomm['pxp_coords'])
    PXP_lis.append(PXP)
    k_pxp_apri = ipxp + kmeta_pxp_apri
    R_pxp_apri = np.random.RandomState(k_pxp_apri)
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

while np.linalg.norm(PXPapri - PXPold) > 5* 10**-3 and iiter < iitermax:
#while np.linalg.norm(kriter - kriterold) > 10**-4:
    print("============ iter No" , iiter , '============')

    # Partie ASM
    ObsASM , ModASM  = acls.vectorialize_ASM_multi(PXPapri_lis,ObsASM_lis,Xbato,Z,C,nbprocs=nbproc)
    B_ASM = ObsASM - ModASM
    JacobASM = acls.jacob_ASM(PXPapri_lis,ObsASM_lis,Xbato,Z,C,h,nbprocs=nbproc)

    #Partie BL
    if with_BL:
        ObsBL,ModBL = acls.vectorialize_BL(BL,PXPapri_lis)
        B_BL = ObsBL - ModBL
        JacobBL =  acls.jacob_BL(PXPapri_lis)

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
    
    # ==== A partir de lÃ  on bosse ac des matrix
    A = np.matrix(A)
    B = np.matrix(B).T
    P = np.matrix(P)

    N = A.T * P * A
    Ninv = scipy.linalg.inv(N) 

    dX = Ninv * A.T * P * B
    
#    c,resid,rank,sigma = scipy.linalg.lstsq(A,B)
    # ==== fin de la section des matrix  
    
    dX = np.array(dX).squeeze()
    PXPold  = PXPapri    
    PXPnew  = PXPapri + dX
    PXPapri = PXPnew
    PXPapri_lis = list(PXPnew.reshape(shape_arr))


    ObsASM_V , ModASM4_V = acls.vectorialize_ASM_multi(PXPapri_lis,ObsASM_lis,Xbato,Z,C,nbprocs=nbproc)
    V = ObsASM_V - ModASM4_V
    fuv = geok.fuv_calc(V,A)
    print('f.u.v. :' , fuv)
    
    PXPnew_arr = PXPnew.reshape(shape_arr)
    print('new :',PXPnew_arr)
    print(PXPnew_stk.append(PXPnew))
    print('sigmas : ' , np.sqrt(np.diag(N)))
    print('sigmas 2 : ' , np.sqrt(np.diag(Ninv) * fuv)) 
    


    
    kriterold = kriter
    kriter = np.linalg.norm(PXPnew - PXPold)
    print("kriter , kriterold", kriter , kriterold)
    print("ecart Ã  la postion vraie : ", end=' ')
    print(np.linalg.norm( PXPnew_arr - PXPtrue_arr,axis=1))
    iiter=iiter+1

print(" BLs Vraies ")
for cpl in itertools.combinations(PXPtrue_arr,2):
    print(np.linalg.norm(cpl[0] - cpl[1])) 

if with_BL:   
    print("BLs ObservÃ©es ")
    for bl in ObsBL:
        print(bl)

print(" BL des coords a prioris('Mod') ")
for cpl in itertools.combinations(PXPapri0_arr,2):
    print(np.linalg.norm(cpl[0] - cpl[1]))

print(" BLs News Ã  l'issue de l'inversion")
for cpl in itertools.combinations(PXPnew_arr,2):
    print(np.linalg.norm(cpl[0] - cpl[1]))

print("Barycentre")
print(acls.barycenter_calc(PXPnew.reshape(shape_arr)))

#F.close()