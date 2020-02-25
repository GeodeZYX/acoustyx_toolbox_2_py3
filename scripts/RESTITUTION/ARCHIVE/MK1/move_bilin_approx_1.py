# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 13:48:43 2015
@author: psakicki
"""

from megalib import *

#  _____  _                     _   _                      _ 
# |  __ \(_)                   | | (_)                    | |
# | |  | |_ ___  ___ ___  _ __ | |_ _ _ __  _   _  ___  __| |
# | |  | | / __|/ __/ _ \| '_ \| __| | '_ \| | | |/ _ \/ _` |
# | |__| | \__ \ (_| (_) | | | | |_| | | | | |_| |  __/ (_| |
# |_____/|_|___/\___\___/|_| |_|\__|_|_| |_|\__,_|\___|\__,_|
# 
#13 juin
#go MK2_bilin_approx_demo.script.py                                                            

# le coin des fonctions

# ============================================================================

#  _      _                       
# | |    (_)                      
# | |     _ _ __   ___  __ _ _ __ 
# | |    | | '_ \ / _ \/ _` | '__|
# | |____| | | | |  __/ (_| | |   
# |______|_|_| |_|\___|\__,_|_|   
#                                 
   
if False:
    path = '/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox/working/data_1/'
    
    M = np.loadtxt(path + 'Mout2000.txt')
    ObsASM = M[:,3]
    Xbato = list(M[:,:3])
    Z = np.loadtxt(path + 'Zout2000.txt')
    C = np.loadtxt(path + 'Cout2000.txt')
    PXP = np.array([-2500,-2500,4000])
                              
    zmin = 0.
    zmax = 9999. #np.max(Z)
    
    Xold = np.array([9999999,9999999])
    g_apri  = 0.01
    cs_apri = 1500

    X = np.array([  1.73109904e-02 ,  1.47089998e+03 ])
    X = np.array([g_apri,cs_apri])
    
    zmin = 0.
    zmax = 9999. #np.max(Z)
    
    while np.linalg.norm(X - Xold) > 10**-3:
        print('une nouvelle boucle ...', np.linalg.norm(X - Xold))
        ObsASM_line = []
        gdiff_stk = []
        csdiff_stk = []
        
        for xbato in Xbato:
            theta , r , t = raytrace_seek_linear(xbato,PXP,X[0],X[1],zmin,zmax,
                                                 1,89,False,False)
            ObsASM_line.append(t)
            # derivation
    #        args = [xbato,PXP,X[0],X[1],zmin,zmax,1,89,False,False]
            kwargs = {'Xsrc': xbato,
                      'Xrec':PXP,
                      'g':X[0],
                      'cs':X[1],
                      'zmin':zmin,
                      'zmax':zmax,
                      'thetaminin':1,
                      'thetamaxin':89,
                      'verbose':False,
                      'fulloutput':False}
    
            gdiff  = geok.partial_derive(raytrace_seek_linear,'g',-1,kwargs,(),h=0)
            csdiff = geok.partial_derive(raytrace_seek_linear,'cs',-1,kwargs,(),h=0)
                        
            gdiff_stk.append(gdiff)
            csdiff_stk.append(csdiff)
        
        aa1 = np.array(gdiff_stk)
        aa2 = np.array(csdiff_stk)
        
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
            theta , r , t = raytrace_seek_linear(xbato,PXP,X[0],X[1],zmin,zmax,False,False)
            ObsASM_resid.append(t)
        ObsASM_resid = np.array(ObsASM_resid)
        Delta = ObsASM_resid - ObsASM 
        print('sum(v**2)', np.sum(Delta**2))



#  ____  _   _      _                       
# |  _ \(_) | |    (_)                      
# | |_) |_  | |     _ _ __   ___  __ _ _ __ 
# |  _ <| | | |    | | '_ \ / _ \/ _` | '__|
# | |_) | | | |____| | | | |  __/ (_| | |   
# |____/|_| |______|_|_| |_|\___|\__,_|_|   
#  


fct = acls.raytrace_seek_bi_linear
Delta_final_stk = []
zb_list = list(range(200,300,300))

outpath = '/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox/working/data_1/Deltas/'


if True:
    for zb in zb_list:
        path = '/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox/working/data_1/'
        
        M = np.loadtxt(path + 'Mout2000.txt')
        ObsASM = M[:,3]
        Xbato = list(M[:,:3])
        Z = np.loadtxt(path + 'Zout2000.txt')
        C = np.loadtxt(path + 'Cout2000.txt')
        PXP = np.array([-2500,-2500,4000])
        
        g1_apri  = -0.05276368770770182 #-1.73109904e-02
        g2_apri  =  0.012686172714815484 #1.73109904e-02
        
        g1_apri  = -0.05 #-1.73109904e-02
        g2_apri  =  0.01 # 1.73109904e-02

        cs = C[0]
    #    zb = Z[np.argmin(C)]
    #    zb = 400
    
        X = np.array([g1_apri,g2_apri])
        
        zmin = 0.
        zmax = 9999. #np.max(Z)
        
        Xold = np.array([9999999,9999999])
        
        #Plot
    #    plt.clf()
#        plt.plot(C,-Z)
        
        iloop = 0
        
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
    
#            plt.plot(Cline,-Zline)        
            
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
                theta , r , t = fct(xbato,PXP,X[0],X[1],cs,zb,zmin,zmax,verbose=False,fulloutput=False)
                ObsASM_resid.append(t)
                
            ObsASM_resid = np.array(ObsASM_resid)
            
            Delta = ObsASM_resid - ObsASM 
            print('sum(v**2)', np.sum(Delta**2))
    
        np.savetxt(outpath + 'zb' + str(zb),Delta)
        Delta_final_stk.append(Delta)

print('fin')
