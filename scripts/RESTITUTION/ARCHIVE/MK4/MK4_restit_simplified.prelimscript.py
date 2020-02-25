# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 09:30:12 2015

@author: pierre
"""

from megalib import *
import raytrace as rt
import sympy

M = np.loadtxt("/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/exemple/SSP/SSP_NOAA_dep5781_20050827000000")

filin = '/home/psakicki/THESE/DATA/1506_GEODESEA/CTD_GEODESEA/CTD/PROCESSED/CTD4/CTD_geodesea_200615_16h15.asc'
Zssp,Cssp = ssp.read_CTD_gedesea(filin)


Vstk = []
Dstk = []
Tstk = []
Estk = np.arange(-60,60,1)
Estkstk = []
Zstk = np.arange(4000 , 5400, 100)

with_delta_extern = 1 # a laisser vrai !!!
with_tempo_evol   = 1
with_time_window  = 1

for z in Zstk:
    for a in Estk:
        X,Z,T = rt.raytrace_SD1_frontend(M[:,0],M[:,1],a,z)
        D = np.linalg.norm([X,Z])
        V = D / T
        Estkstk.append(a)
        Vstk.append(V)
        Dstk.append(D)
        Tstk.append(T)
        
def SSP_mean(Z,C,a,zmax):
    Z = np.array(Z)
    C = np.array(C)

    X,Z,T = rt.raytrace_SD1_frontend(Z,C,a,zmax)
    D = np.linalg.norm([X,Z])
    cout = D / T
    return cout

plt.clf()
plt.plot(Estkstk,Vstk,'+')

sym_var_tup = sympy.symbols('xbato ybato zbato xpxp ypxp zpxp kc c0 tau t t0')
xbato , ybato , zbato , xpxp , ypxp , zpxp , kc , c0 , tau , t , t0 = sym_var_tup
sym_var_lis = list(sym_var_tup)

expr_2  = sympy.sqrt((xbato - xpxp) ** 2 + (ybato - ypxp) ** 2 + (zbato - zpxp) ** 2)  / (kc * (t-t0) + c0)
expr_1a = sympy.sqrt((xbato - xpxp) ** 2 + (ybato - ypxp) ** 2 + (zbato - zpxp) ** 2)  / (c0)


if not with_tempo_evol:
    if with_delta_extern:
        fctobs = expr_1a
    else:
        fctobs = expr_1
else:
    fctobs = expr_2

fctobs_lambda = sympy.lambdify(sym_var_lis,fctobs, "numpy")  

fctobs_diff_list = []
fctobs_diff_lambda_list = []

for sym_var in sym_var_lis:
    fctobs_diff = sympy.diff(fctobs,sym_var)
    fctobs_diff_lambda = sympy.lambdify(sym_var_lis,fctobs_diff, "numpy") 
    fctobs_diff_list.append(fctobs_diff)  
    fctobs_diff_lambda_list.append(fctobs_diff_lambda)

gene_path = "/home/psakicki/THESE/DATA/1506_GEODESEA/WORKING/OPERA_DATA"
exp  = 'IUEM_LAPTOP-3Beacontracking'

exp_path = os.path.join(gene_path,exp)
bigdico  = acls.give_me_the_path(exp_path,exp)#,[1,3,5])

ObsASM_lis  = []
Xbato_lis   = [] 
Tbato_lis   = []  
TTT_lis     = []
PXP_lis     = []
PXPapri_lis = []

kmeta_pxp_apri = 1

with_monoZ = 0

for ipxp,Mpxp in bigdico['N'].items():
    Mdata  = Mpxp['d']
    Mcomm  = Mpxp['c']
    
    print('Mdata' , Mdata.shape)
    
    Xbato  = list(Mdata[:,1:4])
    Tbato  = list(Mdata[:,0])

    Xbato_lis.append(Xbato)
    Tbato_lis.append(Tbato)
    
    ObsASM_load = Mdata[:,-1]

    ObsASM_lis.append(ObsASM_load)
        
    PXP = acls.pxp_string_2_array(Mcomm['pxp_coords'])
    PXP_lis.append(PXP)
    nPXP = len(PXP_lis)

    k_pxp_apri = ipxp + kmeta_pxp_apri
    R_pxp_apri = np.random.RandomState(k_pxp_apri)
    if with_monoZ:
        PXPapri_mono = PXP + np.concatenate((R_pxp_apri.randn(3)[0:2] * sigma_pxp_apri , err_z))
    else:
        PXPapri_mono = PXP
    PXPapri_lis.append(PXPapri_mono)

c0_lis0 = [SSP_mean(Zssp,Cssp,10,PXPapri_lis[-1][-1])] * 3
c0_lis  = list(c0_lis0)
kc_lis  = [0] * 3 


if with_time_window:
    ObsASM_lis_orig = list(ObsASM_lis)
    Xbato_lis_orig  = list(Xbato_lis)
    ObsASM_lis  = []
    Xbato_lis   = []   
    
    ssss = dt.datetime(2015,6,22,3)
    eeee = dt.datetime(2015,6,22,7)

    ssss = dt.datetime(2015,6,21,22)
    eeee = dt.datetime(2015,6,22,3)


    
    ssssp = geok.dt2posix(ssss)
    eeeep = geok.dt2posix(eeee)
    
    for ttt,xbato in zip(Tbato_lis,Xbato_lis_orig):
        _ , outxbato = geok.time_win_basic(ssssp,eeeep,ttt,xbato)
        Xbato_lis.append(outxbato)
    for ttt,obsasm in zip(Tbato_lis,ObsASM_lis_orig):
        _ , outobsasm = geok.time_win_basic(ssssp,eeeep,ttt,obsasm)
        ObsASM_lis.append(outobsasm)

for iiter in range(5):
    Jacobstk  = []
    Bstk      = []
    for ipxp , (Xbato , ObsASM , pxpapri , Tbato) in enumerate(zip(Xbato_lis ,
                                                         ObsASM_lis ,
                                                         PXPapri_lis,
                                                         Tbato_lis)):
        lines_stk = []
    
        for xbato , oasm , t in zip(Xbato , ObsASM , Tbato):
            # xbato xpxp kc c0 tau t
            if not with_tempo_evol:
                argstup = (xbato[0],xbato[1],xbato[2],pxpapri[0],pxpapri[1],
                           pxpapri[2],0,c0_lis[ipxp],oasm,0,0)
            else:
                argstup = (xbato[0],xbato[1],xbato[2],pxpapri[0],pxpapri[1],
                           pxpapri[2],kc_lis[ipxp],c0_lis[ipxp],oasm,t,Tbato[0])
    
            if not with_delta_extern:
                Bstk.append(fctobs_lambda(*argstup))
            else:
                Bstk.append(oasm - fctobs_lambda(*argstup))
    
            xpxpdiff = fctobs_diff_lambda_list[3](*argstup)
            ypxpdiff = fctobs_diff_lambda_list[4](*argstup)
            zpxpdiff = fctobs_diff_lambda_list[5](*argstup)  

            kcdiff   = fctobs_diff_lambda_list[6](*argstup) 
            c0diff   = fctobs_diff_lambda_list[7](*argstup) 
            
            if with_tempo_evol:
                line = [ xpxpdiff , ypxpdiff , zpxpdiff , c0diff , kcdiff ]
            else:
                line = [ xpxpdiff , ypxpdiff , zpxpdiff , c0diff ]
            
            lines_stk.append(line)
        
        Jacobstk.append(np.vstack(lines_stk))
            
    A    = scipy.linalg.block_diag(*Jacobstk)
    B    = np.array(Bstk)
    
    
    Ninv = scipy.linalg.inv((A.T).dot(A))
    dX   = Ninv.dot(A.T).dot(B)


    

    
    print(dX)
    if not with_tempo_evol:    
        dXreshape = dX.reshape((3,4))
    else:
        dXreshape = dX.reshape((3,5))

    PXPapri_lis = list(np.array(np.array(PXPapri_lis) + dXreshape[:,0:3]))
    c0_lis      = np.array(c0_lis) + dXreshape[:,3]
    if with_tempo_evol:
        kc_lis = np.array(kc_lis) + dXreshape[:,4]

V = B - A.dot(dX)

print("final PXP",PXPapri_lis)

fuv = geok.fuv_calc(V,A,1/((10 **-6) **2),1)
Sigma = np.sqrt((np.diag(Ninv) * fuv))
