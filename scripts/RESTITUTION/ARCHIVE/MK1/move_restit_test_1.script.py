# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 15:32:22 2015

@author: psakicki
"""


from megalib import *
path = "/home/pierre/Documents/CODES/acoustyx_toolbox/working/data_2/"
path = '/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox/working/data_20/'


M = np.loadtxt(path + 'Mout2000.txt')
ObsASM = M[:,3]
Xbato = list(M[:,:3])

Z = np.loadtxt(path + 'Zout2000.txt')
C = np.loadtxt(path + 'Cout2000.txt')
PXPapri0 = np.array([-2412,-2892,3010])
PXPapri0 = np.array([-51,14,3878])
PXPapri0 = np.array([-2,1,4001])
PXPapri0 = np.array([-2480,-2515,4015])

PXPapri = PXPapri0
PXPold  = np.array([9999999,9999999,9999999])
kriter = 0
kriterold = 999999


PXPnew_stk = []

while np.linalg.norm(PXPapri - PXPold) > 10**-4:
    kriter = np.linalg.norm(PXPapri - PXPold)
    print('kriter norm :' , np.linalg.norm(PXPapri - PXPold))
    ObsASM , ModASM  = acls.linearize_ASM_multi(PXPapri,ObsASM,Xbato,Z,C,nbprocs=8)
    JacobASM = acls.jacob_ASM(PXPapri,ObsASM,Xbato,Z,C,nbproc=8)

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
    print('new :',PXPnew)
    PXPnew_stk.append(PXPnew)