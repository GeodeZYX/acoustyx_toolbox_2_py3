# -*- coding: utf-8 -*-
"""
Created on Sun Mar  1 18:28:47 2015

@author: psakicki
"""

import acouclass as acls
import raytrace as rt
import numpy as np
import scipy.optimize as optimize


# ==================== zone de test ====================

Zbazik = np.arange(0,5000,1)
Cbazik1 = 1500 + 20 * np.cos(Zbazik  * (0.001))
Cbazik2 = 1500 + 20 * np.sin(Zbazik  * (0.001))

SSPbazik1 = acls.SSP(Zbazik,Cbazik1,e=0)
SSPbazik2 = acls.SSP(Zbazik,Cbazik2,e=10)
zm,cm = rt.munk(6000)

SSPmunk=acls.SSP(zm,cm,e=0)

#SSPbazik1.plot()
#SSPbazik2.plot()

#ssf1 = acls.make_SSF3D_from_SSP(SSPbazik1,xmin,xmax,ymin,ymax,xstep,ystep)
#ssf2 = acls.make_SSF3D_from_SSP(SSPbazik2,xmin,xmax,ymin,ymax,xstep,ystep)
ssf_munk = acls.make_SSF3D_from_SSP(SSPmunk)

ssf_munk.Z


#fig = mlab.figure()
#ssf1.plot()
#time.sleep(2)
##mlab.figure()
#mlab.clf()
#ssf2.plot()

#ssf4d1 = acls.SSF4D([ssf2])
#ssf4d1.add_ssf_in_ssf3d_list(ssf1)
#ssf4d1.ssf3d_list

#mlab.figure()
#ssf_out.plot()
#mlab.figure()
#ssf2.plot()

Zextract , Cextract = ssf_munk.make_a_SSP(obj=False)

z0=0
r0=0
theta0 = 80

theta_4_integ =  theta0
theta_4_clasik = - theta0
smax = 3456
smax = 4456
h=100

h1000 = 1000
iposi = [0,0,0]
iangle = [-theta0,90]
restylis = ['rkck']# , 'rk2'  ,  , 'rkck' ]
restylis = ["rk2",'rk4','rkf45','rkck','rkdp','euler']# , 'rk2'  ,  , 'rkck' ]

resultlis = []
result1000lis = []

#resty='euler'
#YYY,s,c,t,HH = acls.raytrace_ODE_2or3d(ssf_munk,iposi,iangle,smax,h,resotype=resty,return_mode='full')

for resty in restylis:
    print('=====' , resty , '=====')
    restuptmp = acls.raytrace_ODE_2or3d(ssf_munk,iposi,iangle,smax,h,resotype=resty,return_mode='full')
#    restuptmp1000 = acls.raytrace_ODE_2or3d(ssf_munk,iposi,iangle,smax,h1000,resotype=resty,return_mode='full')
    resultlis.append(restuptmp)
#    result1000lis.append(restuptmp1000)
    
resultlis[1][1]
       
t =  2.28687167063
t = 6
#BIGTUP   = rt.raytrace_SnellDesc2(Zextract,Cextract,zstart=0,zmax=max(Zextract),tmax=t,theta=theta_4_clasik)
#X_clasik , Z_clasik , T_clasik  = rt.raytrace_SnellDesc_frontend(Zextract,Cextract,0,max(Zextract),tmax=t,theta=theta_4_clasik)
print('=====','classic','=====')
BIGTUP = rt.raytrace_SD3_frontend(Zextract,Cextract,0,max(Zextract),tmax=t,theta=-theta_4_clasik,cumsum = 1)
X_clasik , Z_clasik , T_clasik = BIGTUP

np.savetxt('/home/psakicki/Zmunk1',Zextract)
np.savetxt('/home/psakicki/Cmunk1',Cextract)

Z_ml = np.loadtxt('/home/psakicki/Zmunk_matlab_1')
X_ml = np.loadtxt('/home/psakicki/Xmunk_matlab_1')

X_py = np.cumsum(BIGTUP[0])
Z_py = np.cumsum(BIGTUP[1])

#plt.clf()
#plt.plot(X_ml,-Z_ml,'-g+')
#plt.plot(X_py,-Z_py,'rx-')
#plt.plot(YYY[:,1],-YYY[:,2],'*')
#plt.plot(np.cumsum(BIGTUP[2]),-np.cumsum(BIGTUP[0]),'b+')
#for res,restype in zip(resultlis,restylis):
#    plt.plot(res[0][:,1],-res[0][:,2],'o',label=restype)
#plt.legend()
#plt.axis('equal')

#fig = plt.figure()
#plt.clf()
#ax = fig.add_subplot(111, projection='3d')
#ax.scatter(Y23d[:,0],Y23d[:,1],Y23d[:,2])
#plt.show()

#%%

Xe = [0,0,0]
Xr = [2500,2500,4000]
acls.orthogonal_proj(Xr,[0,0,0],[1000,0,5000])    

xapri = acls.canonical_shooting_angle(Xe,Xr)

Xref =  Xr

import time

# Solution iterative 1
strt = time.time()
sol = optimize.root(acls.raytrace_ODE_wrapper,xapri,args=(Xe,Xr,ssf_munk,3000,'euler'),method='hybr',tol=10)
sol.x
sol2 = optimize.root(acls.raytrace_ODE_wrapper,sol.x,args=(Xe,Xr,ssf_munk,3000,'rkdp'),method='hybr',tol=10**-8)
sol2.x
sol3 = optimize.root(acls.raytrace_ODE_wrapper,sol2.x,args=(Xe,Xr,ssf_munk,3000,'rkdp'),method='hybr',tol=10**-7)
sol3.x
t1 = time.time() - strt
print(t1)

rt.raytrace_seek(np.array(Xe),np.array(Xr),Zextract,Cextract,5,85)

# Solution iterative 2 => donne direct à manger les aprioris à un raytracer opérationnel
strt = time.time()
#sol = optimize.root(acls.raytrace_ODE_wrapper,xapri,args=(Xe,Xr,ssf_munk,3000,'euler'),method='hybr',tol=10)
sol4 = optimize.root(acls.raytrace_ODE_wrapper,xapri,args=(Xe,Xr,ssf_munk,1,'rkf45'),method='hybr',tol=10**-7)
sol4.x
t2 =  time.time() - strt
print(t2)

acls.raytrace_ODE_wrapper(sol4.x,Xe,Xr,ssf_munk,50,'rkf45')

#sollist = []
#methlist = []
#tlist = []
#for meth in ('hybr','lm','broyden1','broyden2','anderson','linearmixing','diagbroyden','excitingmixing'):#,'krylov'):
#    try:
#        t = time.time()
#        print meth
#        sol = optimize.root(raytrace_ODE_wrapper,x,args=(Xe,Xr,ssf_munk,100,'euler'),method=meth,,tol=0.0001)
#        sollist.append(sol)
#        methlist.append(meth)
#        tlist.append(time.time() - t)
#    except:
#        continue
#print 'cest fini'