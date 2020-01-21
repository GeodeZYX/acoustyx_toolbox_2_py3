# -*- coding: utf-8 -*-
"""
Created on Sat Jun 13 18:27:22 2015

@author: pierre
"""

from megalib import *

path_gene = '/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox/working'
exp  = 'test2'
exp  ='compar_SSPorig_SSPbilin_nonoise'
exp  ='compar_SSPorig_SSPbilin_noised'
exp  = 'test6'

path_exp  = os.path.join(path_gene,exp)
bigdico = acls.give_me_the_path(path_exp,exp)#,[2,3,4])

zrec = 4000
Zin = bigdico['Z']['d']
Cin = bigdico['C']['d'] 


# Avoir les gradients
# zb == 0 => recherche du zb optimal
zb = 700
zb = 0
max_loop = 7
#acls.SSP_2_bilin_grads(Zin,Cin,zrec,zb=zb)
# Front end
Zout , Cout = acls.SSPreal_2_SSPbilin(Zin,Cin,zrec,zb=zb,max_loop=max_loop)

plt.plot(Cout,-Zout)
plt.plot(Cin , -Zin)

Zout_path = '/'.join(( path_exp , exp + '.2lin_Z.dat'))
Cout_path = '/'.join(( path_exp , exp + '.2lin_C.dat'))
header = bigdico['Z']['c']['raw']
np.savetxt(Zout_path,Zout,header=header)
header = bigdico['C']['c']['raw']
np.savetxt(Cout_path,Cout,header=header)