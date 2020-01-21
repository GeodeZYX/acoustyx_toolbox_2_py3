# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 13:04:35 2015

@author: psakicki
"""

from megalib import *

if platform.node() == 'calipso':
    gene_path = "/home/psakicki/THESE/DATA/1506_GEODESEA/WORKING/OPERA_DATA"
    gene_path = '/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox/working'
    gene_path = '/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox/working'
elif platform.node() == 'msi':
    gene_path = '/home/pierre/Documents/CODES/acoustyx_toolbox//working'


exp = "tst_3x300_x2000_y2000_nois+True_"

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
iitermax = 12
h = 0 
nbproc = 4
PXPold = np.array([9999999,9999999,9999999] * Npxp) 
PXPnew_stk = []
ObsASM_lis  = []
PXP_lis = []
PXPapri_lis = []

# PARAMETRES IMPORTANTS
with_BL = 0
with_dtepoc = 1
with_ssp_bilin = 0
with_mono_epoc = 1

fuv = 1

if with_BL:
    BL = bigdico['B']['d']

if with_ssp_bilin:
    Z = bigdico['2lin_Z']['d']
    C = bigdico['2lin_C']['d'] 
else:
    Z = bigdico['Z']['d']
    C = bigdico['C']['d']

sigma_pxp_apri = 10**-3
kmeta_pxp_apri = 1 # le K normal c'est l'ID du PXP, le meta c'est pour tous
                   # le Kmeta sera un multiple de 10 par convention
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
    PXPapri_mono = PXP + R_pxp_apri.randn(3) * sigma_pxp_apri
    PXPapri_lis.append(PXPapri_mono)

Nepoc = len(Xbato)

PXPtrue_arr = np.array(PXP_lis)
shape_arr   = PXPtrue_arr.shape
# 1st apriori parameter, won't gonna change
PXPapri0_arr = np.array(PXPapri_lis)
PXPapri0_vct = genefun.vectorialize(PXPapri0_arr)
# apriori vector which gonna change at each iter
PXPapri = PXPapri0_vct 
# apriori for dt epoc
if not with_mono_epoc:
    dt_apri0 = [0.] * Nepoc
    dt_apri0 = Mpxp['d'][:,-1]
else:
    dt_apri0 = [ 0.00001 ]
    dt_apri0 = [ Mpxp['d'][0,-1] + 2*10**-6 ]
    
dt_apri = dt_apri0

# Apriori vector which gonna change at every iteration
Xapri = np.array(list(PXPapri0_vct)  + list(dt_apri0))
Xapri0 = Xapri

A_stk = []

# ===============================
# BOUCLE D'INVERSION
# ===============================

while np.linalg.norm(PXPapri - PXPold) > 5 * 10**-5 and iiter < iitermax:
    #while np.linalg.norm(kriter - kriterold) > 10**-4:
    print("============ iter No" , iiter , '============')

    # Partie ASM
#    ObsASM , ModASM  = acls.vectorialize_ASM_multi(PXPapri_lis,ObsASM_lis,
#                                                   Xbato,Z,C,nbprocs=nbproc,
#                                                   dtepoch_lis=dt_apri0 + 10**-4 * np.random.randn(len(dt_apri0)))

    ObsASM , ModASM  = acls.vectorialize_ASM_multi(PXPapri_lis,ObsASM_lis,
                                                   Xbato,Z,C,nbprocs=nbproc,
                                                   dtepoch_lis=dt_apri)

    B_ASM = ObsASM - ModASM
    
    ObsASM2 , ModASM2  = acls.vectorialize_ASM_multi(PXPapri_lis,ObsASM_lis,
                                               Xbato,Z,C,nbprocs=nbproc)
    B_ASM2 = ObsASM2 - ModASM2
    
    print("B_ASM 1&2" , np.sum(B_ASM) , np.sum(B_ASM2))

                                               
    JacobASM = acls.jacob_ASM(PXPapri_lis,ObsASM_lis,Xbato,Z,C,h,nbprocs=nbproc)

    #Partie BL
    if with_BL:
        ObsBL,ModBL = acls.vectorialize_BL(BL,PXPapri_lis)
        B_BL = ObsBL - ModBL
        JacobBL =  acls.jacob_BL(PXPapri_lis)
        
    #partie dtepoc
    if with_dtepoc:
        if not with_mono_epoc:
            lineS = []
            for i in range(Nepoc):
                lin = Nepoc * Npxp * [0]
                for j in range(Npxp):
                    lin[i + j*Nepoc] = 1
                lineS.append(lin)
        else:
            lineS = [[1] * Nepoc * Npxp]
            
        JacobDT = genefun.array_from_lists(*lineS)
        
    # Partie Commune
    A = JacobASM
    B = B_ASM
    P = np.eye(len(B))

    if with_BL:
        K,Q,P = geok.weight_mat([10,1],[len(ObsASM),len(ObsBL)],fuv)
        A = np.vstack((A,JacobBL))
        B = np.hstack((B,B_BL))
    
    if with_dtepoc:
        A = np.hstack((A,JacobDT))

    A_stk.append(A)
    
    # ==== A partir de là on bosse ac des matrix
    A = np.matrix(A)
    B = np.matrix(B).T
    P = np.matrix(P)

    N = A.T * P * A
#    Ninv = np.linalg.inv(N) 
    Ninv = scipy.linalg.inv(N) 

    dX = Ninv * A.T * P * B

#    c,resid,rank,sigma = scipy.linalg.lstsq(A,B)
    # ==== fin de la section des matrix  
    
    dX = np.array(dX).squeeze()
    Xold = Xapri
    Xnew = Xapri + dX
    Xapri = Xnew
    
    
    PXPold  = PXPapri    
    PXPnew  = Xnew[:np.product(PXPtrue_arr.shape)]
    PXPapri = PXPnew
    PXPapri_lis = list(PXPnew.reshape(shape_arr))
    
    dt_apri = Xnew[np.product(PXPtrue_arr.shape):]

    ObsASM_V , ModASM4_V = acls.vectorialize_ASM_multi(PXPapri_lis,ObsASM_lis,
                                                       Xbato,Z,C,nbprocs=nbproc,
                                                       dtepoch_lis = dt_apri)
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
    print(PXPnew_stk.append(PXPnew))
    print('sigmas 2 : ' , np.sqrt(np.diag(Ninv) * fuv)[:np.product(PXPtrue_arr.shape)])
    
    kriterold = kriter
    kriter = np.linalg.norm(PXPnew - PXPold)
    print("kriter , kriterold", kriter , kriterold)
    print("ecart à la postion vraie : ", end=' ')
    print(np.linalg.norm( PXPnew_arr - PXPtrue_arr,axis=1))
    
    print("Barycentre")
    print(acls.barycenter_calc(PXPnew.reshape(shape_arr)))

    iiter=iiter+1

print(" Résumé ")
print("Nb pings     " , Nepoc)
print("Nb PXPs      " , Npxp)
print("Taille Jacob." , A.shape)
print("Taille Résid." , V.shape)

print(" BLs Vraies ")
for cpl in itertools.combinations(PXPtrue_arr,2):
    print(np.linalg.norm(cpl[0] - cpl[1])) 

if with_BL:   
    print("BLs Observées ")
    for bl in ObsBL:
        print(bl)

print(" BL des coords a prioris('Mod') ")
for cpl in itertools.combinations(PXPapri0_arr,2):
    print(np.linalg.norm(cpl[0] - cpl[1]))

print(" BLs News à l'issue de l'inversion")
for cpl in itertools.combinations(PXPnew_arr,2):
    print(np.linalg.norm(cpl[0] - cpl[1]))



#F.close()