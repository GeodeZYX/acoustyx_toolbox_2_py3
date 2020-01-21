# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 14:05:37 2016

@author: psakicki

Permet d'exploiter les statistiques 
produites par la 2nde partie de 
/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/scripts/analyz_SSP_MOVE_NOAA/move_plot_new_style.py
"""

from megalib import *
reload(gf)

dicopath = '/home/psakicki/THESE/bd_CTD/NEW_DATA_2015/MOVE/MOVE_reform_1507/SIMPLE/statistics_dico.pik'
D = gf.pickle_loader(dicopath)
        
Dall = D['all']

Tall = geok.ymdhms_vectors2dt(Dall['y'],Dall['m'],Dall['d'],
                  Dall['mi'],Dall['h'],Dall['s'])
Tall = np.unique(Tall)

Tlis = []
for k , ite in D.items():
    T = geok.ymdhms_vectors2dt(ite['y'],ite['m'],ite['d'],
                  ite['mi'],ite['h'],ite['s'])
    Tlis.append(T)

Tcommon = gf.find_common_elts(*Tlis)
argsrtstk = []

maxstdstk = []
stdtrustk = []

for k , ite in D.items():
    if k == 'all':
        continue
    T = geok.ymdhms_vectors2dt(ite['y'],ite['m'],ite['d'],
                  ite['mi'],ite['h'],ite['s'])
    print('sensor', k , 'std max :' , np.max(ite['std']) , \
    'amplitude max :' , np.max(ite['amplitude']))
    maxstdstk.append(np.max(ite['std']))
    boolz = np.array([ t in Tcommon for t in T ])
    Ttru   = np.array(T[boolz])
    stdtru = np.array(ite['std'][boolz])
    stdtrustk.append(stdtru)
    
    argsrt = (np.argsort(stdtru))
    argsrtstk.append(argsrt)
    
### 2 moyens de trouver les min/max/med
    
if 0: # je sais pas ce qu'il fait discontinué
    S = np.sum(np.column_stack(argsrtstk) ,1)
    
    zeind_best  = np.argwhere(S == np.max(S)).squeeze()
    print("best")
    Tbest   = (Tcommon)[zeind_best]
    print(Tbest , S[zeind_best])
    zeind_worst = np.argwhere(S == np.min(S)).squeeze()
    print("worst")
    Tworst  = (Tcommon)[zeind_worst]
    print(Tworst , S[zeind_worst])
    zeind_median = np.argwhere(S == np.median(S)).squeeze()
    print("median")
    Tmedian = (Tcommon)[zeind_median]
    print(Tmedian , S[zeind_median])

elif 1:
    S = np.mean(np.column_stack(stdtrustk),1)
    
    zeind_best  = np.argwhere(S == np.min(S)).squeeze()
    print("best")
    Tbest   = (Tcommon)[zeind_best]
    print(Tbest , S[zeind_best])
    zeind_worst = np.argwhere(S == np.max(S)).squeeze()
    print("worst")
    Tworst  = (Tcommon)[zeind_worst]
    print(Tworst , S[zeind_worst])
    zeind_median = np.argwhere(S == np.median(S)).squeeze()
    print("median")
    Tmedian = (Tcommon)[zeind_median]
    print(Tmedian , S[zeind_median])


else: # trouve le même pire jour et meilleur jour que la methode précédente mais trouve un autre jour médian
    S = np.sum(np.stack(stdtrustk[:5]),0)
    
    zeind_best  = np.argwhere(S == np.min(S)).squeeze()
    print("best")
    Tbest   = (Tcommon)[zeind_best]
    print(Tbest , S[zeind_best])
    zeind_worst = np.argwhere(S == np.max(S)).squeeze()
    print("worst")
    Tworst  = (Tcommon)[zeind_worst]
    print(Tworst , S[zeind_worst])
    zeind_median = np.argwhere(S == np.median(S)).squeeze()
    print("median")
    Tmedian = (Tcommon)[zeind_median]
    print(Tmedian , S[zeind_median])







Twork  = Tbest
Twork  = Tworst

Dwork  = D['all']
Dwork2 = Dwork[(Dwork['y'] == Twork.year ) & (Dwork['m'] == Twork.month) & (Dwork['d'] == Twork.day)]

year   = Twork.year
month  = Twork.month
day    = Twork.day

Dwork2 = Dwork[(Dwork['y'] == year) & (Dwork['m'] == month) & (Dwork['d'] == day)]