# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 13:04:35 2015

@author: psakicki
"""

from megalib import *

path = '/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox/working/test6'
exp = 'test6'
bigdico = acls.give_me_the_path(path,exp)

# ==============================================
# INITIALISATION DES PARAMETRES
# ==============================================

# PARAMETRES ANNEXES
kriter = 10**7
kriterold = 10**10
iiter = 1
iitermax = 2
h = 0 
nbproc = 8
PXPold  = np.array([9999999,9999999,9999999])
PXPnew_stk = []



# PARAMETRES IMPORTANTS
ipxp   = 1 
M      = bigdico['M'][ipxp]['d']
Z      = bigdico['Z']['d']
C      = bigdico['C']['d']
Zold = Z
Cold = C

if 0:
    SSP = np.loadtxt('/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox/exemple/SSP/SSP_NOAA_dep5781_20050827000000')
    Zin = SSP[:,0]
    Cin = SSP[:,1]
    Z = Zin
    C = Cin


Com    = bigdico['M'][ipxp]['c']
ObsASM = M[:,3]
Xbato  = list(M[:,:3])
PXP    = acls.pxp_string_2_array(Com['pxp_coords'])


kmeta_pxp_apri = 1
k_pxp_apri = ipxp + kmeta_pxp_apri
sigma_pxp_apri = 10
R_pxp_apri = np.random.RandomState(k_pxp_apri)
PXPapri = PXP + R_pxp_apri.randn(3) * sigma_pxp_apri

# ===============================
# BOUCLE D'INVERSION
# ===============================

#while np.linalg.norm(PXPapri - PXPold) > 10**-4:
while np.linalg.norm(kriter - kriterold) > 10**-4 and iiter < iitermax:
    print("iter No" , iiter)
    ObsASM , ModASM  = acls.vectorialize_ASM_multi(PXPapri,ObsASM,Xbato,Z,C,nbprocs=nbproc)
    JacobASM = acls.jacob_ASM(PXPapri,ObsASM,Xbato,Z,C,h,nbprocs=nbproc)

    A = JacobASM
    B = ObsASM - ModASM

    # ==== A partir de là on bosse ac des matrix
    A = np.matrix(A)
    B = np.matrix(B).T

    N = A.T * A
    Ninv = scipy.linalg.inv(N) 

    dX = Ninv * A.T * B
    c,resid,rank,sigma = scipy.linalg.lstsq(A,B)
    
    # ==== fin de la section des matrix    
    dX = np.array(dX).squeeze()

    PXPold  = PXPapri    
    PXPnew  = PXPapri + dX
    PXPapri = PXPnew
    print('new :',PXPnew)
    PXPnew_stk.append(PXPnew)
    
    kriterold = kriter
    kriter = np.linalg.norm(PXPnew - PXPold)
    print("kriter , kriterold", kriter , kriterold)
    print("ecart à la postion vraie : " , \
    np.linalg.norm(PXPnew - PXP)) 
    iiter=iiter+1