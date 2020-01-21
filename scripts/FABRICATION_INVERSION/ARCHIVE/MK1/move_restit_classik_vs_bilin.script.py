# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 15:32:22 2015

Ce script est un merge d'une restituton classque avec un SSP 'figé' et d'un SSP
bilinéaire (figé également, mais qui gère la variation)

@author: psakicki
"""


from megalib import *
path = '/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox/working/data_1/'

M = np.loadtxt(path + 'Mout2000.txt')
ObsASM = M[:,3]
Xbato = list(M[:,:3])

Z = np.loadtxt(path + 'Zout2000.txt')
C = np.loadtxt(path + 'Cout2000.txt')
PXPapri0 = np.array([-2412,-2892,3010])
PXPapri0 = np.array([-2480,-2515,4015])


#   _____ _               _ _    
#  / ____| |             (_) |   
# | |    | | __ _ ___ ___ _| | __
# | |    | |/ _` / __/ __| | |/ /
# | |____| | (_| \__ \__ \ |   < 
#  \_____|_|\__,_|___/___/_|_|\_\

print("CLASSIK PART")
                               
                            
PXPapri = PXPapri0
PXPold  = np.array([9999999,9999999,9999999])

PXPnew_stk = []
PXPnew_stk.append(PXPapri)


while np.linalg.norm(PXPapri - PXPold) > 10**-4:
    ObsASM , ModASM  = acls.linearize_ASM_multi(PXPapri,ObsASM,Xbato,Z,C)
    JacobASM = acls.jacob_ASM(PXPapri,ObsASM,Xbato,Z,C)

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
    
    ObsASM_resid = []
    for xbato in Xbato:
        theta , r , t = rt.raytrace_seek(xbato,PXPnew,Z,C,1,89,False,False)
        ObsASM_resid.append(t)
    ObsASM_resid = np.array(ObsASM_resid)
    Delta = ObsASM_resid - ObsASM 
    print('sum(v**2) for classik restit', np.sum(Delta**2))    
    
    PXPnew_stk.append(PXPnew)

PXP_stk_classik = list(PXPnew_stk)
    
#  ____  _ _ _       
# |  _ \(_) (_)      
# | |_) |_| |_ _ __  
# |  _ <| | | | '_ \ 
# | |_) | | | | | | |
# |____/|_|_|_|_| |_|
#
    
print("BILINEAR PART")

Xbato = list(M[:,:3])
PXP = np.array([-2500,-2500,4000])

g1_apri  = -0.05276368770770182 #-1.73109904e-02
g2_apri  =  0.01268617271481548 # 1.73109904e-02
cs = C[0]
zb = Z[np.argmin(C)]

X = np.array([g1_apri,g2_apri])

zmin = 0.
zmax = 9999. #np.max(Z)

Xold = np.array([9999999,9999999])

#Plot
#plt.clf()
#plt.plot(C,-Z)

iloop = 0

fct = acls.raytrace_seek_bi_linear

while np.linalg.norm(X - Xold) > 5* 10**-4:
    iloop = iloop + 1
    print('une nouvelle boucle ...', np.linalg.norm(X - Xold))
    ObsASM_line = []
    g1diff_stk = []
    g2diff_stk = []
    cbdiff_stk = []
    zbdiff_stk = []
    
    Zline = np.array([zmin,zb,zmax+100])
    Cline = acls.fct_bi_linear(Zline,X[0],X[1],cs,zb)[0]

#    plt.plot(Cline,-Zline,label=str(iloop))        
    
    for xbato in Xbato:
        # derivation
        kwargs_bilin = {'Xsrc': xbato,
                      'Xrec':PXP,
                      'g1':X[0],
                      'g2':X[1],
                      'cs':cs,
                      'zb':zb,
                      'zmin':zmin,
                      'zmax':zmax,
                      'thetaminin':1,
                      'thetamaxin':89,
                      'verbose':False,
                      'fulloutput':False}
                  
        theta , r , t = fct(**kwargs_bilin)
        ObsASM_line.append(t)
        
        g1diff = geok.partial_derive(fct,'g1',-1,kwargs_bilin,())
        g2diff = geok.partial_derive(fct,'g2',-1,kwargs_bilin,())
        zbdiff = geok.partial_derive(fct,'zb',-1,kwargs_bilin,())

        g1diff_stk.append(g1diff)
        g2diff_stk.append(g2diff)
    
    aa1 = np.array(g1diff_stk)
    aa2 = np.array(g2diff_stk)
    
    ObsASM_line = np.array(ObsASM_line)
    B = ObsASM - ObsASM_line
    
    A = np.vstack((aa1,aa2)).T
    
    At = A.T
    N = np.dot(At,A)
    Ninv = scipy.linalg.inv(N) 
    
    dX = np.dot(np.dot(Ninv,At),B)
    dX2 = scipy.linalg.lstsq(A,B)[0]
    
    Xold = X
    Xnew = X + dX
    X = Xnew
    
    print('Xold :' ,Xold)
    print('Xnew :' ,Xnew)
    
    ObsASM_resid = []
    for xbato in Xbato:
        theta , r , t = fct(xbato,PXP,X[0],X[1],cs,zb,zmin,zmax,False,False)
        ObsASM_resid.append(t)
        
    ObsASM_resid = np.array(ObsASM_resid)
    
    Delta = ObsASM_resid - ObsASM 
    print('sum(v**2)', np.sum(Delta**2))     

    
PXPapri = PXPapri0
PXPold  = np.array([9999999,9999999,9999999])

PXPnew_stk = []
PXPnew_stk.append(PXPapri)

PXP

Z = np.array([zmin,zmax+100])
C = acls.fct_bi_linear(Z,X[0],X[1],cs,zb)[0]


while np.linalg.norm(PXPapri - PXPold) > 10**-4:
    ObsASM , ModASM  = acls.linearize_ASM_multi(PXPapri,ObsASM,Xbato,Z,C)
    JacobASM = acls.jacob_ASM(PXPapri,ObsASM,Xbato,Z,C)

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

    ObsASM_resid = []
    for xbato in Xbato:
        theta , r , t = rt.raytrace_seek(xbato,PXPnew,Z,C,1,89,False,False)
        ObsASM_resid.append(t)
    ObsASM_resid = np.array(ObsASM_resid)
    Delta = ObsASM_resid - ObsASM 
    print('sum(v**2) for bilin restit', np.sum(Delta**2))    
    
    PXPnew_stk.append(PXPnew)


PXP_stk_bilin = list(PXPnew_stk)

print('the end')
