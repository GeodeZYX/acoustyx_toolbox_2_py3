# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 15:09:22 2015

@author: psakicki
"""

import acouclass as acls
from megalib import *



gradmod = 1

path_gene = "/home/psakicki/THESE/DATA/1506_GEODESEA/WORKING/OPERA_DATA"
exp  = 'gsea_1'
path_gene = '/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox/working'
exp  = 'test6'

path_exp  = os.path.join(path_gene,exp)
bigdico = acls.give_me_the_path(path_exp,exp)#,[2,3,4])

zrec = 4000
Zin = bigdico['Z']['d']
Cin = bigdico['C']['d'] 

Zold = Zin
Cold = Cin

Zreal = Zin
Creal = Cin

if 1:
    SSP = np.loadtxt('/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox/exemple/SSP/SSP_NOAA_dep5781_20050827000000')
    Zin = SSP[:,0]
    Cin = SSP[:,1]
    
    Zreal = Zin
    Creal = Cin

# Avoir les gradients
# zb == 0 => recherche du zb optimal
zb = 800
max_loop = 7
#acls.SSP_2_bilin_grads(Zin,Cin,zrec,zb=zb)
# Front end
Zout , Cout ,Grad = acls.SSPreal_2_SSPbilin(Zin,Cin,zrec,zb=zb,max_loop=max_loop,
                                            outgrads=1)
                                                                             
                                        
Xsrc = np.array([0,0,0])
Xrec = np.array([0,500,-2500])

F = rt.raytrace_seek_bilin_input

vars_of_rt = (Xsrc,Xrec,Grad[0],Grad[1],Cin[0],zb,0,88,0,0)
kwargs_of_rt = {'Xrec' : Xrec,
                'Xsrc' : Xsrc,
                'g1'   : Grad[0],
                'g2'   : Grad[1],
                'cs'   :  Cin[0],
                 'zb'  :  zb   ,
                 "thtmin" : 0,
                 "thtmax" : 88,
                 "verbose":0,
                 "fulloutput":0
                }

rt.raytrace_seek_bilin_input(**kwargs_of_rt)

geok.partial_derive(F,'g1',-1,kwargs_of_rt)
geok.jacobian(F,[2],-1,[kwargs_of_rt],h=10**-6)

F(Xsrc,Xrec,Grad[0],Grad[1],Cin[0],zb,0,88,0,0)

# retour aux données de simu

# PARAMETRES ANNEXES
kriter = 10**7
kriterold = 10**10
iiter = 1
iitermax = 10
h = 0 
nbproc = 8
PXPold  = np.array([9999999,9999999,9999999])
PXPnew_stk = []

# PARAMETRES IMPORTANTS
ipxp   = 1 
M      = bigdico['M'][ipxp]['d']
#Z      = bigdico['Z']['d']
#C      = bigdico['C']['d']
Com    = bigdico['M'][ipxp]['c']
ObsASM = M[:,3]
Xbato  = list(M[:,:3])
PXP    = acls.pxp_string_2_array(Com['pxp_coords'])


kmeta_pxp_apri = 1
k_pxp_apri = ipxp + kmeta_pxp_apri
sigma_pxp_apri = 10
R_pxp_apri = np.random.RandomState(k_pxp_apri)
PXPapri = PXP + R_pxp_apri.randn(3) * sigma_pxp_apri



if gradmod:
    APRI = np.concatenate((PXPapri , Grad))
else:
    APRI = PXPapri

APRInew_stk = []

# ===============================
# BOUCLE D'INVERSION
# ===============================

print("début de l'inversion")

#while np.linalg.norm(PXPapri - PXPold) > 10**-4:
while np.linalg.norm(kriter - kriterold) > 10**-5 and iiter < iitermax:
    if gradmod:
        Z,C = acls.bilin_grads_2_SSP(APRI[-2],APRI[-1],Cin[0],zb)
    else:
        Z,C = Zreal , Creal

    print("iter No" , iiter)
    ObsASM , ModASM  = acls.vectorialize_ASM_multi(APRI[0:3],ObsASM,Xbato,Z,C,nbprocs=nbproc)
    JacobASM = acls.jacob_ASM(APRI[0:3],ObsASM,Xbato,Z,C,h,nbprocs=nbproc)
    
    if gradmod:
        # JACOBGRAD
        kwdic_gene = {'Xrec' : APRI[0:3],
                    'Xsrc' : Xbato,
                    'g1'   : APRI[-2],
                    'g2'   : APRI[-1],
                    'cs'   :  Cin[0],
                    'zb'  :  zb   ,
                    "thtmin" : 0,
                    "thtmax" : 88,
                    "verbose":0,
                    "fulloutput":0}          
        kwdic_var = {'Xsrc' : Xbato}
        KWARGSLIS = geok.kwargs_for_jacobian(kwdic_gene,kwdic_var)
        JacobGRAD = geok.jacobian(F,[2,3],-1,KWARGSLIS,h=0)

    if gradmod:  
        A = np.hstack((JacobASM , JacobGRAD)) # 
    else:
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

#    PXPold  = PXPapri    
#    PXPnew  = PXPapri + dX
#    PXPapri = PXPnew
#    print 'new :',PXPnew
#    PXPnew_stk.append(PXPnew)

    APRIold = APRI
    APRInew = APRI + dX
    APRI    = APRInew
    print('new :' , APRI)
    APRInew_stk.append(APRInew)
    
    kriterold = kriter
    kriter = np.linalg.norm(APRInew[0:3] - APRIold[0:3])
    print("kriter , kriterold", kriter , kriterold)
    print("ecart à la postion vraie : " , \
    np.linalg.norm(APRInew[0:3] - PXP)) 
    iiter=iiter+1



#hhh = 1.05132300233e-09
#
#aq = F(Xsrc,Xrec,Grad[0]-hhh,Grad[1],Cin[0],zb,0,88,0,0)
#zs = F(Xsrc,Xrec,Grad[0]+hhh,Grad[1],Cin[0],zb,0,88,0,0)
#
#(aq[-1] - zs[-1]) / (2*hhh)
#

