# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 15:32:22 2015

@author: psakicki
"""

from megalib import *
path = '/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox/working/data_2_4_PXP'

Mlis = []
for f in glob.glob(os.path.join(path,'Mout*')):
    M = np.loadtxt(f)
    Mlis.append(M)

Z = np.loadtxt(glob.glob(os.path.join(path,'Z*'))[0])
C = np.loadtxt(glob.glob(os.path.join(path,'C*'))[0])

ObsASM_lis = []
PXPapri_lis = [] 

for M in Mlis:
    ObsASM_lis.append(M[:,3])
    PXPapri_lis.append(M[0,10:13])
    
    Xbato = list(M[:,:3])
    
PXPnew_stk = []

PXPapri = np.array([])
for pxp in PXPapri_lis:
    PXPapri = np.append(PXPapri,pxp)   

PXPold  = 9999999 * np.ones(3 * len(PXPapri_lis) )

while np.linalg.norm(PXPapri - PXPold) > 1 * 10**-4:
    
    print('new loop ... critere:', np.linalg.norm(PXPapri - PXPold))
    
    PXPapri = np.array([])
    for pxp in PXPapri_lis:
        PXPapri = np.append(PXPapri,pxp)   
        
    ObsASM , ModASM = acls.linearize_ASM_multi(PXPapri_lis,ObsASM_lis,Xbato,Z,C,nbprocs=8)
    JacobASM = acls.jacob_ASM(PXPapri_lis,ObsASM_lis,Xbato,Z,C,nbproc=8)
    
    A = JacobASM
    B = ObsASM - ModASM
    
    At = A.T
    N = np.dot(At,A)
    Ninv = scipy.linalg.inv(N) 
    
    dX = np.dot(np.dot(Ninv,At),B)
    c,resid,rank,sigma = scipy.linalg.lstsq(A,B)
    PXPold  = PXPapri    
    PXPnew  = PXPapri + dX
    PXPapri = PXPnew
    PXPnew_stk.append(PXPnew)
    
    PXPapri_lis = [ np.array(e) for e in genefun.chunkIt(list(PXPapri),len(PXPapri_lis)) ]
