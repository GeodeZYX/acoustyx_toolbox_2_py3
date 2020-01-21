# -*- coding: utf-8 -*-
"""
Created on Thu Mar  3 16:05:06 2016

@author: psakicki
"""

from megalib import *


STEP = np.arange(0,3,.2)[:-1]
ZCi  = [rt.munk(3001,10**(-i)) for i in STEP]

fig , axlis = plt.subplots(2,2)

switchtup = (0,1)
kkkktup   = (-1,0)

for igraph , (switch , kkk) in enumerate(itertools.product(switchtup,kkkktup)):
    
    ax = axlis.T.flat[igraph]

    if switch:
        ANG  = np.arange(10,80,10)
        ZZZ  = [2000]
        NUM_COLORS = len(ANG)
    else:
        ANG  = [45]
        ZZZ  = np.arange(500,3001,500)
        NUM_COLORS = len(ZZZ)
    
    cm = plt.get_cmap('viridis_r')
    cool = [cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)]
    
    Bstk = []
    icool = -1

    for a in ANG:
        for z in ZZZ:
            icool = icool+1
            stklis = []
            for ik , i in enumerate(STEP):
                Z,C = ZCi[ik]
                #Z ,C  = rt.munk(3001,10**(-i))
                stklis.append(rt.raytrace_SD1_frontend(Z,C,a,z))
            A = np.vstack(stklis)
            B  = A[:,kkk] - A[-1,kkk]
            Bb = A[:,kkk] 
            medi = gf.median_improved(Bb)
            mean = np.mean(Bb)
            B  = Bb - mean
            Bstk.append(B)
            if switch:
                lab = str(a) + 'deg.'
            else:
                lab = str(z) + 'm.'
            ax.plot(np.power(10,-STEP),B,'o-',label=lab,c=cool[icool])
#            ax.scatter(np.power(10,-STEP)[B == 0][0],B[B == 0][0],marker='*',
#                       c=cool[icool],s=500)
                       
            print(B == 0)

            ax.set_xscale('log')

    ax.set_xlabel('z-sampling of the SSP (m)')
    if kkk == -1:
        ax.set_ylabel('time propagation difference to the mean value (s)')
    else:
        ax.set_ylabel('radial propagation difference to the mean value (m)')

    ax.legend(loc=0)

    ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

axlis[0,1].set_title("propagation parameters for different SSP sampling and \n different shooting angles at constant depth (2000m)")
axlis[0,0].set_title("propagation parameters for different SSP sampling and \n different depths at constant shooting angle (45deg)")

fig.suptitle('Comparison of Snell-Descartes raytracings in function of the Sound Speed Profile sampling to the mean value at constant emission angles and final depths')

Bstk2 = np.vstack(Bstk)
np.std(Bstk2,0)
np.argwhere(np.std(Bstk2,0) == np.min(np.std(Bstk2,0)))

for i in range(100):
    rt.raytrace_SD1_frontend(Z,C,a,z)
   
fig.set_size_inches(11.69,8.27)
pltpath = '/home/psakicki/THESE/RENDU/1511_raytracing/graphs/SSP/optimal_SD.pdf'
plt.savefig(pltpath)

