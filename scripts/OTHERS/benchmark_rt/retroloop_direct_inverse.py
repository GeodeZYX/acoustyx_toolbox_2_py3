#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 12:09:57 2017

@author: psakicki
"""

from megalib import *

plt.ion()

# ========================= PARTIE 0 (easy) =========================
# /home/psakicki/Téléchargements/delta_t_post_gradient.py
from megalib import *
import raytrace as rt

Z,C  = rt.munk(6000)

delta = 4 *  10**-5

outtup1 = rt.raytrace_ultimate(0,0,45,2.8,Z,C)
outtup2 = rt.raytrace_ultimate(0,0,45,2.8 + delta,Z,C)

xxf1 , zzf1 , ttf1 , ddf1 , ttf1 = outtup1
xxf2 , zzf2 , ttf2 , ddf2 , ttf2 = outtup2

xxf1[-1] - xxf2[-1]

zzf1[-1] - zzf2[-1]


geok.dist((xxf1[-1],zzf1[-1]) , (xxf2[-1],zzf2[-1]))


# ========================= PARTIE 1 (easy) =========================

import raytrace as rt

Z,C = rt.munk(6000)

TUP0 = rt.raytrace_ultimate(0,0,45,5,Z,C)
X0   = TUP0[0][-1]
Z0   = TUP0[1][-1]

fig , ax = plt.subplots(1,2)

D_T = np.logspace(-5,-4,6)
D_T = [0,2e-5,4e-5,6e-5,8e-5,1e-4,1.2e-4]

for d_t in D_T:
    
    TUP2 = rt.raytrace_ultimate(0,0,45,5 + d_t,Z,C)

    X2   = TUP2[0][-1]
    Z2   = TUP2[1][-1]
    
    ax[0].scatter(d_t , X2 - X0)
    ax[1].scatter(d_t , Z2 - Z0)

# ========================= PARTIE 2 (less easy) =========================
# provient de /home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/scripts/benchmark_rt/exploit_scripts/inverse/exploit_rt_inverse_gradient_realist.py

if gf.get_computer_name() == 'calipso':
    dir_prefix = '/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/'
else:
    dir_prefix = '/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/'


# FIRST VERSION
path = "/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/scripts/benchmark_rt/results/20160316_052748benmk_inverse_grd.pik"
# REBOOT 1703
protopath =  dir_prefix + '/scripts/benchmark_rt/results_REBOOT1703/inverse/*pik'

Lpath = glob.glob(protopath)

ZZZ , CCC = rt.munk(6000)

for path in Lpath:
    L = gf.pickle_loader(path)
    
    plotdir = "/home/psakicki/THESE/RENDU/1603_graphs_benchmk_rt/indirect_gradient"
    plotdir = "/home/psakicki/THESE/RENDU/1703_graphs_benchmk_rt_REBOOT/indirect_gradient"
    plotdir = dir_prefix + "/scripts/benchmark_rt/results_REBOOT1703/PLOT/"
    
    plotdir = gf.create_dir(plotdir)
    
    linlog = "linear"
    linlog = "log"
    
    if linlog == 'linear':
        logsuffix = ''
    else:
        logsuffix = ' (log. scale)'
        
    
    logsuffix = ''
    
    params_lis        = []
    eiko_grad_lis     = []
    eiko_proto_lis    = []
    eiko_x_grad_lis   = []
    eiko_x_proto_lis  = []
    diff_lis          = []
    diff_t_lis        = []
    tup_lis           = []
    
    t_raw_eiko_lis    = []
    t_raw_proto_lis   = []
    
    
    for i_exp,exp in enumerate(L):
        params          = exp[0]
        try:
            eiko        = exp[2]
            eiko_grad   = eiko[0:3]
            eiko_proto  = eiko[3:6]
        except:
            continue
        
        params_lis.append(params)
        
        eiko_grad_lis.append(eiko_grad)
        eiko_proto_lis.append(eiko_proto)
        eiko_x_grad_lis.append(eiko_grad[0].x)
        eiko_x_proto_lis.append(eiko_proto[0].x)   
        diff_lis.append(eiko_grad[0].x - eiko_proto[0].x)
        diff_t_lis.append(eiko_grad[1][-1] - eiko_proto[1][-1])

        t_raw_proto_lis.append(eiko_proto[1][-1])
        t_raw_eiko_lis.append(eiko_grad[1][-1])

        
        tup_lis.append(params[:-1] + tuple(diff_lis[-1]) + (diff_t_lis[-1],) + \
                      (t_raw_proto_lis[-1],) + (t_raw_eiko_lis[-1],) + \
                      (eiko_grad[0].x[0],)   + (eiko_proto[0].x[0],))
    
    DF = pd.DataFrame(tup_lis)
    
    rnamedic = gf.renamedic_fast_4_pandas("h","restype","adap",
                                        "xr","yr","zr","zsmooth",
                                        "smoothtype","xgrad","ygrad",
                                        "diff_ang1","diff_ang2",
                                        "diff_s","diff_t",
                                        "t_raw_proto","t_raw_eiko",
                                        "ang_proto", "ang_eiko")
    
    DF = DF.rename(columns = rnamedic)
    
    ####################################################

#fig_disc  , ax_disc  = plt.subplots(1,1)
fig_discx , ax_discx = plt.subplots(1,1)
fig_discz , ax_discz = plt.subplots(1,2)

fig_discxz , [ax_discx , ax_discz] = plt.subplots(1,2)

#ax_discx2 = ax_discx.twiny()
#ax_discz2 = ax_discz.twiny()

pdw       = copy.deepcopy(DF)
symzr_lis = list(reversed(gf.symbols_list(np.unique(pdw['zr']))))
colzr_lis = list(reversed(gf.color_list(np.unique(pdw['zr']))))


#col_lis = gf.color_list(pdw['t_raw_proto'])
sym_lis = gf.symbols_list(pdw['diff_t'])

import matplotlib.pyplot as plt
cm = plt.get_cmap('viridis')




for zr , colzr in zip(np.unique(pdw['zr']),colzr_lis):

    astk   = []
    dtstk  = []

    xxxstk = []
    zzzstk = []
    
    xxx_as_abs_stk = []
        
    
    pdww = pdw[pdw['zr'] == zr]

    dt_plot_lis = []
    #### PARTIE POUR LES GRAPH DE LA DISCUSSION 
    for ii , (a , tr , dt , xxx_as_abs) in enumerate(zip(np.abs(pdww['ang_eiko'])  ,
                                                   pdww['t_raw_eiko']     ,
                                                   pdww['diff_t']         ,
                                                   pdww['xr'])):
        
        print('angle : ' , a)
        #dt = 4 *  10**-5
        XXX0 , ZZZ0 , TTT0 , DDD0 , _ = rt.raytrace_ultimate(0,0, a , tr + 0 ,ZZZ,CCC)
        XXX1 , ZZZ1 , TTT1 , DDD1 , _ = rt.raytrace_ultimate(0,0, a , tr + dt,ZZZ,CCC)
        
        xxx0 = XXX0[-1]
        xxx1 = XXX1[-1]
        zzz0 = ZZZ0[-1]
        zzz1 = ZZZ1[-1]

#        new_tick_locations = sorted(pdww['diff_t'])
#        
#        def tick_function(X):
#            return ["%.3f" % z for z in X]
#        
#        ax_discx2.set_xlim(ax_discx.get_xlim())
#        ax_discx2.set_xticks(new_tick_locations)
#        ax_discx2.set_xticklabels(tick_function(new_tick_locations))
#        
        
        astk.append(a) 
        dtstk.append(dt) 

        xxxstk.append((xxx1 - xxx0))
        zzzstk.append((zzz1 - zzz0))
        
        xxx_as_abs_stk.append(xxx_as_abs)
        
    ax_discx.plot( xxx_as_abs_stk , xxxstk ,'o-',  c=colzr , label = str(zr) + ' m')
    ax_discz.plot( xxx_as_abs_stk , zzzstk ,'o-',  c=colzr , label = str(zr) + ' m')
    
ax_discx.set_ylabel("coordinate difference along the X component (m)")
ax_discz.set_ylabel("coordinate difference along the Z component (m)")

ax_discx.set_xlabel("reciever depth (m)")
ax_discz.set_xlabel("reciever depth (m)")

ax_discx.legend(loc=0,title=r'$\bf{Depth}$')
ax_discx.legend().get_frame().set_alpha(0.5)

ax_discz.legend(loc=0,title=r'$\bf{Depth}$')
ax_discz.legend().get_frame().set_alpha(0.5)

plotdir = '/home/psakicki/aaa_FOURBI'


fig_discxz.set_size_inches(fig_discxz.get_size_inches() * np.array([2,1]))

fig_discxz.suptitle('Difference on radial and vertical components for a standard Snell Descartes raytracing, \n with or without time propagation delay induced by a sound speed gradient')

plt.figure(2)
pltname = gf.join_improved('_', 'retroloop','X')
pltpath = os.path.join( plotdir , pltname )
plt.savefig(pltpath + '.png')
plt.savefig(pltpath + '.pdf')

plt.figure(3)
pltname = gf.join_improved('_', 'retroloop','Z')
pltpath = os.path.join( plotdir , pltname )
plt.savefig(pltpath + '.png')
plt.savefig(pltpath + '.pdf')

plt.figure(4)
pltname = gf.join_improved('_', 'retroloop','XZ')
pltpath = os.path.join( plotdir , pltname )
plt.savefig(pltpath + '.png')
plt.savefig(pltpath + '.pdf')
 
    #ax_discx2.plot( dtstk , xxxstk ,'o-',  c=colzr)
    #ax_discz2.plot( dtstk , zzzstk ,'o-',  c=colzr)    

#En X : une double échelle ANGLE et TEMPS (parce qu'il sont liées)
#Et des coleurs pour la profondeur mais sans passer par une colorbar
# En Y : bah le delta X et Y

#plt.close('all')