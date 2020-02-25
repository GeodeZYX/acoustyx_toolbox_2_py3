# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 13:48:43 2015
@author: psakicki
"""

import numpy as np
import scipy
from scipy import interpolate
import scipy.optimize as optimize
import time
import datetime as dt
import dateutil
import os
import inspect
import itertools
import multiprocessing as mp
import glob

import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
#import mayavi.mlab as mlab

#lib persos geodezyx
import genefun
import geodetik as geok

# libs persos acoustic
import acouclass as acls
import raytrace as rt
import SSP as ssp


#  _____  _                     _   _                      _ 
# |  __ \(_)                   | | (_)                    | |
# | |  | |_ ___  ___ ___  _ __ | |_ _ _ __  _   _  ___  __| |
# | |  | | / __|/ __/ _ \| '_ \| __| | '_ \| | | |/ _ \/ _` |
# | |__| | \__ \ (_| (_) | | | | |_| | | | | |_| |  __/ (_| |
# |_____/|_|___/\___\___/|_| |_|\__|_|_| |_|\__,_|\___|\__,_|
# 
#13 juin
#go MK2_bilin_approx_demo.script.py                                                            





#  ____  _   _      _                       
# |  _ \(_) | |    (_)                      
# | |_) |_  | |     _ _ __   ___  __ _ _ __ 
# |  _ <| | | |    | | '_ \ / _ \/ _` | '__|
# | |_) | | | |____| | | | |  __/ (_| | |   
# |____/|_| |______|_|_| |_|\___|\__,_|_|   
#  



# Trouver le zb optimal
# V1
#for zb in zb_list:
#    G , V = find_bilin_grads(Z,C,Xbato,ObsASM,PXP,zb=zb,verb=0)
#    
#    Zline = np.array([min(Z),zb,max(Z)+100])
#    Cline = acls.fct_bi_linear(Zline,G[0],G[1],C[0],zb)[0]
#
#    Zline_stk.append(Zline)
#    Cline_stk.append(Cline)
#
#
#    print G
# V2

#plt.figure()
#ax =plt.gca()
#ax.set_yscale('log')
#plt.plot(zb_list[4:-2],sum_stk[4:-2])

#plt.figure()
#plt.plot(C,-Z)
#for z,c in zip(Zline_stk,Cline_stk):
#    plt.plot(c,-z)

#
#if True:
#    for zb in zb_list:
#        M = np.loadtxt(path + 'Mout2000.txt')
#        ObsASM = M[:,3]
#        Xbato = list(M[:,:3])
#        Z = np.loadtxt(path + 'Zout2000.txt')
#        C = np.loadtxt(path + 'Cout2000.txt')
#
#        PXP = np.array([-2500,-2500,4000])
#        PXP = np.array([-00,-00,4000])
#        
#        g1_apri  = -0.05276368770770182 #-1.73109904e-02
#        g2_apri  =  0.012686172714815484 #1.73109904e-02
#        
##        g1_apri  = -0.05 #-1.73109904e-02
##        g2_apri  =  0.01 # 1.73109904e-02
##        g1_apri = -0.05149047 
##        g2_apri = 0.01300313
#
##        g1_apri = 1
##        g2_apri = 0.1
##        
#        g1_apri  = -0.1 #-1.73109904e-02
#        g2_apri  =  0.01 # 1.73109904e-02
#
#        cs = C[0]
#        
#        X = np.array([g1_apri,g2_apri])
#        
#        zmin = 0.
#        zmax = 9999. #np.max(Z)
#        
#        Xold = np.array([9999999,9999999])
#        
#        iloop = 0
#        
#        while np.linalg.norm(X - Xold) > 1 * 10**-5 and iloop < 20:
#            iloop = iloop + 1
#            print 'une nouvelle boucle ...', np.linalg.norm(X - Xold) , 'no' , iloop
#            ObsASM_line = []
#            g1diff_stk = []
#            g2diff_stk = []
#            cbdiff_stk = []
#            zbdiff_stk = []
#
#            Zline = np.array([zmin,zb,zmax+100])
#            Cline = acls.fct_bi_linear(Zline,X[0],X[1],cs,zb)[0]
#            
#            kwargs_bilin_list = []
#            strt = time.time()
#
#            for xbato in Xbato:
#                # derivation
#                kwargs_bilin = {'Xsrc': xbato,
#                              'Xrec':PXP,
#                              'g1':X[0],
#                              'g2':X[1],
#                              'cs':cs,
#                              'zb':zb,
#                              'zmin':zmin,
#                              'zmax':zmax,
#                              'thetaminin':1,
#                              'thetamaxin':89,
#                              'verbose':False,
#                              'fulloutput':False}
#                              
#                kwargs_bilin_list.append(kwargs_bilin)
#
#                          
#                theta , r , t = fct(**kwargs_bilin)
#                ObsASM_line.append(t)
#                
##                g1diff = geok.partial_derive(fct,'g1',-1,kwargs_bilin,())
##                g2diff = geok.partial_derive(fct,'g2',-1,kwargs_bilin,())
##    
##                g1diff_stk.append(g1diff)
##                g2diff_stk.append(g2diff)
##            
##            aa1 = np.array(g1diff_stk)
##            aa2 = np.array(g2diff_stk)
##            
#            ObsASM_line = np.array(ObsASM_line)
#            B = ObsASM - ObsASM_line
#            
#            strt = time.time()
#            A = geok.jacobian(fct,['g1','g2'],-1,kwargs_f_list=kwargs_bilin_list,nproc=6,h=0)
#            print time.time() - strt
#            
##            A = np.vstack((aa1,aa2)).T
#            
#            At = A.T
#            N = np.dot(At,A)
#            Ninv = scipy.linalg.inv(N) 
#            
#            dX = np.dot(np.dot(Ninv,At),B)
#            dX2 = scipy.linalg.lstsq(A,B)[0]
#            
#            Xold = X
#            Xnew = X + dX
#            X = Xnew
#            
#            print 'Xold :' ,Xold
#            print 'Xnew :' ,Xnew
#
#            Delta = ObsASM_resid - ObsASM 
#            print 'sum(v**2)', np.sum(Delta**2)
#    
#        np.savetxt(outpath + 'zb' + str(zb),Delta)
#        Delta_final_stk.append(np.abs(Delta))
#        print 'sum(v) FINAL', np.sum(Delta) , zb
#        Delta_sum_final_stk.append((zb,np.sum(np.abs(Delta))))
#        Delta_max_final_stk.append((zb,np.max(np.abs(Delta))))
#        
#        g1g2stk.append((X[0],X[1]))
#        

##    
##plt.clf()
##for Deltaa in Delta_final_stk:
##    for x,delta in zip( Xbato,Deltaa):
##        #plt.plot(x[0],x[1],'+')
##        d=np.linalg.norm(x)
##        plt.plot(d,delta,'+')
#        
#
#
#plt.figure()
#plt.yscale('log', nonposy='clip')
#for z,delta in  Delta_max_final_stk:
#    plt.plot(z,delta,'*')    
#
#
##plt.clf()
##plt.plot(C,-Z)
#

print('fin')
