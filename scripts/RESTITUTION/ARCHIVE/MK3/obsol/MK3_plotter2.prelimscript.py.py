# -*- coding: utf-8 -*-
"""
Created on Tue Sep  8 14:27:14 2015

@author: psakicki
"""

from megalib import *
import matplotlib.pyplot as plt


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



dickey2D = 'ecart 3D a la postion vraie en distance'


dickey2D = 'ecart 2D bary brut/vrai en distance'
dickey3D = 'ecart 3D bary brut/vrai en distance'


               
vars_str = ', '.join(variable_lis)

diff2D_lis = []
diff3D_lis = []
lgnd_lis   = []

for d in diclis:
    ilastiter = max(d.keys())
    if not dickey2D in list(d[ilastiter].keys()):
        continue
    if not dickey3D in list(d[ilastiter].keys()):
        continue
    
    diff2D_lis.append(d[ilastiter][dickey2D])
    diff3D_lis.append(d[ilastiter][dickey3D])
    
    boolstr=''.join([str(int(d[0][vari])) for vari \
    in variable_lis if vari in list(d[0].keys())])
    
    lgnd_lis.append(boolstr)


fig,(ax2d,ax3d) = plt.subplots(2,1)

NUM_COLORS = len(diff2D_lis)
cm = plt.get_cmap('gist_rainbow')
colors = [cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)]
ax2d.set_xscale('symlog',linthreshx=0.0001,linscalex=1,subsx=np.arange(0,10,1))
ax3d.set_xscale('symlog',linthreshx=0.0001,linscalex=1,subsx=np.arange(0,10,1))
ax2d.set_xlim([ 0 , 1000 ])
ax3d.set_xlim([ 0 , 1000 ])
ax2d.set_title('2D distance from the true point')
ax3d.set_title('3D distance from the true point')
ax2d.set_xlabel('2D distance (m)')
ax3d.set_xlabel('3D distance (m)')
ax2d.set_ylabel('PXP ID')
ax3d.set_ylabel('PXP ID')

scat2d_lis = []
scat3d_lis = []

for i,d2d in enumerate(diff2D_lis):
    if not genefun.is_iterable(d2d):
        d2d = [d2d]
    
    for j , dd in enumerate(d2d):
        scat = ax2d.scatter(dd,j+1,c=colors[i],s=150,alpha=.5)
        scat2d_lis.append(scat)
        ax2d.annotate(lgnd_lis[i], (dd,j+1))
        ax2d.scatter(0,j+1,c='r',marker='*',s=350,alpha=1)


for i,d3d in enumerate(diff3D_lis):
    if not genefun.is_iterable(d3d):
        d3d = [d3d]
    
    for j , dd in enumerate(d3d):       
        scat = ax3d.scatter(dd,j+1,c=colors[i],s=150,alpha=.5)
        scat3d_lis.append(scat)
        ax3d.annotate(lgnd_lis[i], (dd,j+1))
        ax3d.scatter(0,j+1,c='r',marker='*',s=350,alpha=1)



fig.legend(scat2d_lis,lgnd_lis,'upper right',scatterpoints = 1)
fig.suptitle(dickey2D + '/' + dickey3D + '\n'+ d[ilastiter]['nom'] + '\n' + vars_str)
fig.set_size_inches(11.69,8.27)

plt.show()