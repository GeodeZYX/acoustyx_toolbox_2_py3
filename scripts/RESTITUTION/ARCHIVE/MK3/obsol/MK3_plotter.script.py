# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 16:34:16 2015

@author: psakicki
"""

from megalib import *
import matplotlib.pyplot as plt

expdic_in = "/media/pierre/D34F-24E3/CODES/acoustyx_toolbox_2/working/BATC2/batc2_3x100_x2000_y2000_nois1-1e-06_/*.exp"
expdic_in = "/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/working/BATC2/batc2_3x100_x2000_y2000_nois1-1e-06_/*.exp"

#def aaa(expdics_in,variable_lis = ['with_BL', 'with_monoZ', 'with_alternat_SSP','with_barycenter']):
    
variable_lis = ['with_alternat_SSP', 'with_monoZ', 'with_barycenter', 'with_BL']
    
if type(expdic_in) is str:
    filis = glob.glob(expdic_in)
else:
    filis = expdics_in
    
if len(filis) == 0:
    raise Exception('ERR : no files in the list ...')
    
diclis = []
for f in filis:
    diclis.append(genefun.pickle_loader(f))

PtCalc_lis = []
PtVrai_lis = []
Ptdiff_lis = []
lgnd_lis   = []


calcdiffkey = ('nouvelles coords.',
'ecart a la postion vraie en coordonnees')
calcdiffkey = ("Barycentre 'brut' : Sum Xpxp / Npxp" ,
               "ecart au bary brut/vrai en coords.")
               
vars_str = ', '.join(variable_lis)

for d in diclis:
    ilastiter = max(d.keys())
    if (not calcdiffkey[0] in list(d[ilastiter].keys())) or \
       (not calcdiffkey[1] in list(d[ilastiter].keys())):
        continue
        
    PtCalc = d[ilastiter][calcdiffkey[0]]
    Ptdiff = d[ilastiter][calcdiffkey[1]]
        
    PtVrai = PtCalc - Ptdiff
    
    PtCalc_lis.append(PtCalc)
    Ptdiff_lis.append(Ptdiff)
    PtVrai_lis.append(PtVrai)
    
    boolstr=''.join([str(int(d[0][vari])) for vari \
    in variable_lis if vari in list(d[0].keys())])
        
    lgnd_lis.append(boolstr)

Ptdiff_arr = np.array(Ptdiff_lis)
    
X = Ptdiff_arr
Y = Ptdiff_arr

Ptdiff_arr.shape


if np.ndim(Ptdiff_arr) == 3:
    multipts = True
else:
    multipts = False

npts = len(d[1]['nouvelles coords.'])

if multipts:
    fig , axraw = plt.subplots(npts/2 , npts/2)
    figv , axvraw = plt.subplots(1,npts)
    axtup = axraw.flatten()
    axvtup = axvraw.flatten()

else:
    fig , axraw = plt.subplots() 
    figv , axvraw = plt.subplots()
    axtup = [axraw]
    axvtup = [axvraw]

# PLANI
for ipt,ax in enumerate(axtup): 
    if multipts:
        X = Ptdiff_arr[:,ipt,0]
        Y = Ptdiff_arr[:,ipt,1]
    else:
        X = Ptdiff_arr[:,0]
        Y = Ptdiff_arr[:,1]
    NUM_COLORS = len(X)
    cm = plt.get_cmap('gist_rainbow')
    colors = [cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)]
    scat_lis = []
    ax.scatter(0,0,c='r',marker='*',s=150,alpha=1)
    ax.set_ylim([ - 1000 , 1000 ])
    ax.set_xlim([ - 1000 , 1000 ])
    ax.set_ylabel('Y diff. to the true point (m.)')
    ax.set_xlabel('X diff. to the true point (m.)')
    ax.set_yscale('symlog',linscaley=3,linthreshy=0.01,subsy=np.arange(0,10**-2,10**-3))
    ax.set_xscale('symlog',linscalex=3,linthreshx=0.01,subsx=np.arange(0,10**-2,10**-3))
    ax.grid(True)
    if multipts:
        ax.set_title('PXP no ' + str(ipt+1))
    else:
        ax.set_title('barycenter of ' + str(npts) + ' points')
        
    for i,(x,y) in enumerate((list(zip(X,Y)))):
        scat = ax.scatter(x,y,c=colors[i],s=150,alpha=.5)
        scat_lis.append(scat)
        ax.annotate(lgnd_lis[i], (x,y))

fig.legend(scat_lis,lgnd_lis,'upper right',scatterpoints = 1)
fig.suptitle(calcdiffkey[0] + '(PLANI)' + '\n' +d[ilastiter]['nom'] + '\n' + vars_str)
fig.set_size_inches(11.69,8.27)
plt.show()

#ALTI
for ipt,ax in enumerate(axvtup): 
    if multipts:
        Z = Ptdiff_arr[:,ipt,2]
        X = [0] * len(Z)
    else:
        Z = Ptdiff_arr[:,2]
        X = [0] * len(Z)
    NUM_COLORS = len(Z)
    cm = plt.get_cmap('gist_rainbow')
    colors = [cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)]
    scat_lis = []
    ax.scatter(0,0,c='r',marker='*',s=150,alpha=1)
    ax.set_ylim([ - 100 , 100 ])
    ax.set_xlim([ - 100 , 100 ])
    ax.set_ylabel('Z diff. to the true point (m.)')
    ax.set_yscale('symlog',linscaley=1,linthreshy=0.01,subsy=np.arange(0,10**-2,10**-3))
    ax.set_xscale('symlog',linscalex=1,linthreshx=0.01,subsx=np.arange(0,10**-2,10**-3))
    if multipts:
        ax.set_title('PXP no ' + str(ipt+1))
    else:
        ax.set_title('barycenter of ' + str(npts) + ' points')
    

    ax.grid(True)
    
    for i,(x,z) in enumerate((list(zip(X,Z)))):
        scat = ax.scatter(x,z,c=colors[i],s=150,alpha=.5)
        scat_lis.append(scat)
        ax.annotate(lgnd_lis[i], (x,z))
        
figv.legend(scat_lis,lgnd_lis,'upper right',scatterpoints = 1)
figv.suptitle(calcdiffkey[0] + '(comp. Z)' + '\n'+ d[ilastiter]['nom'] + '\n' + vars_str)
figv.set_size_inches(11.69,8.27)

plt.show()