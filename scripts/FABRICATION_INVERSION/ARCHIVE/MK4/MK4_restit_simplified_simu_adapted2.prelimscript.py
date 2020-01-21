# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 09:30:12 2015

@author: pierre
"""

from megalib import *
import raytrace as rt
import sympy

##Cette config permet des petits dX
## (pas touche et modifie les paramètres opera ci apres)
#with_delta_extern = 1 # a laisser vrai !!!
#with_tempo_evol   = 0
#with_time_window  = 0
#with_monoc0       = 1
#with_alternat_SSP = 0
#with_c0_estim     = 0

with_delta_extern = 1 # a laisser vrai !!!
with_tempo_evol   = 0
with_time_window  = 0
with_monoc0       = 1 # il est vivement recommandé et logique de n'avoir qu'un c0
                      # si on est pas en mode monoc0 on 
                      # on a des results batards
                      # diff true bary 3D (norm) 2154.14831413
                      # diff true bary 2D (norm) 29.7034726047
with_alternat_SSP = 0
with_c0_estim     = 1

gene_path = '/media/psakicki/D32E-8E7E/CODES/acoustyx_toolbox_2/working'
gene_path = '/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/working'

exp = 'batc2c_deriv_decalZ_nopingnoise_5000_R50_noisFalse-0__'

exp = 'batc2c_deriv_decalZ_nopingnoise_trajnoise5decim_5000_R50_noisFalse-0__'
exp = "batc2c_deriv_decalZ_5000_R50_noisTrue-1e-06__"  

exp = 'batc2c_deriv_decalZ_notrajnoise_stdpingnoise_5000_R50_noisTrue-1e-06__'

exp = "batc2c_deriv_decalZ_noping_notraj_noise_5000_R50_noisFalse-0__"
exp = 'batc4_nadir_5000_R50_noisFalse-0__'



exp_path = os.path.join(gene_path,exp)
bigdico  = acls.give_me_the_path(exp_path,exp)#,[1,3,5])

Zssp  = bigdico['Z']['d']
Cssp  = bigdico['C']['d']

if with_alternat_SSP:
    alternat_SSP_path='/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/exemple/SSP/SSP_NOAA_dep5781_20050827000000'
    alternat_SSP = np.loadtxt(alternat_SSP_path)
    Zssp = alternat_SSP[:,0]
    Cssp = alternat_SSP[:,1]




ObsASM_lis  = []
Xbato_lis   = [] 
Tbato_lis   = []  
TTT_lis     = []
PXP_lis     = []
PXPapri_lis = []

kmeta_pxp_apri = 1

with_monoZ = 0

for ipxp,Mpxp in bigdico['M'].items():
    Mdata  =  Mpxp['d']
    Mcomm  =  Mpxp['c']
    
    print('Mdata' , Mdata.shape)
    
    Xbato  = list(Mdata[:,:3])
#    Tbato  = list(Mdata[:,0])

    Xbato_lis.append(Xbato)
#    Tbato_lis.append(Tbato)
    
    ObsASM_load = Mdata[:,3]

    ObsASM_lis.append(ObsASM_load)
        
    PXP = acls.pxp_string_2_array(Mcomm['pxp_coords'])
    PXP_lis.append(PXP)
    Npxp = len(PXP_lis)

    k_pxp_apri = ipxp + kmeta_pxp_apri
    R_pxp_apri = np.random.RandomState(k_pxp_apri)
    
    sigma_pxp_apri = 10.0
    
    if with_monoZ:
        PXPapri_mono = PXP + np.concatenate((R_pxp_apri.randn(Npxp)[0:2] * sigma_pxp_apri , err_z))
    else:
        PXPapri_mono = PXP
        PXPapri_mono = PXP + R_pxp_apri.randn(3) * sigma_pxp_apri

    PXPapri_lis.append(PXPapri_mono)

if not with_monoc0:
    c0_lis0 = [ssp.SSP_mean(Zssp,Cssp,10,zpxp) for zpxp in [e[-1] for e  in PXPapri_lis]]
    c0_lis  = list(c0_lis0)
    kc_lis  = [0] * Npxp
else:
    c0_lis0 = [ ssp.SSP_mean(Zssp,Cssp,0,PXPapri_lis[-1][-1]) ]    
    c0_lis  = list(c0_lis0) 
    
#if with_time_window:
#    ObsASM_lis_orig = list(ObsASM_lis)
#    Xbato_lis_orig  = list(Xbato_lis)
#    ObsASM_lis  = []
#    Xbato_lis   = []   
#    
#    ssss = dt.datetime(2015,6,22,3)
#    eeee = dt.datetime(2015,6,22,7)
#
#    ssss = dt.datetime(2015,6,21,22)
#    eeee = dt.datetime(2015,6,22,3)
#
#    ssssp = geok.dt2posix(ssss)
#    eeeep = geok.dt2posix(eeee)
#    
#    for ttt,xbato in zip(Tbato_lis,Xbato_lis_orig):
#        _ , outxbato = geok.time_win_basic(ssssp,eeeep,ttt,xbato)
#        Xbato_lis.append(outxbato)
#    for ttt,obsasm in zip(Tbato_lis,ObsASM_lis_orig):
#        _ , outobsasm = geok.time_win_basic(ssssp,eeeep,ttt,obsasm)
#        ObsASM_lis.append(outobsasm)

##################  ESTIMATION DU SSP EQUIV ##################
##############################################################

# estimation du dc en fct de l'angle
# ========================================================

zref = 4000
zref = PXPapri_lis[-1][-1]

def func(x,a,b,c):
    return a * np.exp(b * x) + c

cccstk = []
aaastk = []
#for a in range(30):
for a in np.arange(-50,50,1):
    ccc = ssp.SSP_mean(Zssp,Cssp,a,zref)
    cccstk.append(ccc)
    aaastk.append(a)

np.savetxt("/home/psakicki/aaa_FOURBI/cccstk.txt" ,cccstk)
np.savetxt("/home/psakicki/aaa_FOURBI/aaastk.txt" ,aaastk)

popt, pcov = scipy.optimize.curve_fit(func, aaastk, cccstk)


polydeg = 12
niter = 3
Poly = np.polyfit(aaastk, cccstk, polydeg )

fited = []
for a in aaastk:
    #    fited.append(func(a,*popt))
    fited.append(np.polyval(Poly,a))
    
#plt.clf()
#plt.plot(aaastk, cccstk,"*")    
#plt.plot(aaastk,fited,"-+")

# estimation du dc en fct de la prof
# ========================================================

cccstk = []
zzzstk = []
for z in np.arange(-200,200,10):
    ccc = ssp.SSP_mean(Zssp,Cssp,0,zref + z)
    cccstk.append(ccc)
    zzzstk.append(z)
    
k_dz , _ = geok.linear_regression(zzzstk,cccstk)

# ========================================================

sym_var_tup = sympy.symbols('xbato ybato zbato xpxp ypxp zpxp kc c0 tau t t0')
xbato , ybato , zbato , xpxp , ypxp , zpxp , kc , c0 , tau , t , t0 = sym_var_tup
sym_var_lis = list(sym_var_tup)

expr_D    = sympy.sqrt((xbato - xpxp) ** 2 + (ybato - ypxp) ** 2 + (zbato - zpxp) ** 2)  
expr_D2d  = sympy.sqrt((xbato - xpxp) ** 2 + (ybato - ypxp) ** 2)  

expr_angl = sympy.asin((xpxp - zbato) / expr_D) * (180./np.pi)
expr_angl = sympy.acos((zpxp - zbato) / expr_D) * (180./np.pi)

#expr_C = c0 + popt[0] * sympy.exp(expr_angl * popt[1]) - k_dz * (zpxp - zref)

#PTDMODIF B
expr_C   = c0

expr_C   = 1500
expr_C   = np.polyval(Poly,41.47)
Polymod  = np.array(Poly)
Polymod[-1] = 0
expr_C   = np.polyval(Poly,expr_angl) + k_dz * (zpxp - zref)
expr_C   = c0 + np.polyval(Polymod,41.47) + k_dz * (zpxp - zref)
expr_C   = c0 + np.polyval(Polymod,expr_angl) + k_dz * (zpxp - zref)

#PTDMODIF B

#PTDMODIF D
expr_2   = expr_D / (kc * (t-t0) + c0)
expr_1a  = expr_D / (c0)
expr_1b  = expr_D / expr_C
#PTDMODIF D

#expr_3  = expr_D_sym / ((c0 + dc) 

#PTDMODIF C
if not with_tempo_evol:
    if with_delta_extern:
        fctobs = expr_1a
else:
    fctobs = expr_2
    
fctobs = expr_1b
#PTDMODIF C

expr_angl_lambda = sympy.lambdify(sym_var_lis,expr_angl, "numpy")  
fctobs_lambda    = sympy.lambdify(sym_var_lis,fctobs, "numpy")  

fctobs_diff_list = []
fctobs_diff_lambda_list = []

for sym_var in sym_var_lis:
    fctobs_diff = sympy.diff(fctobs,sym_var)
    fctobs_diff_lambda = sympy.lambdify(sym_var_lis,fctobs_diff, "numpy") 
    fctobs_diff_list.append(fctobs_diff)  
    fctobs_diff_lambda_list.append(fctobs_diff_lambda)


################## FIN ESTIMATION DU SSP EQUIV ##################
#################################################################

for iiter in range(niter):
    print("===== iter no " , iiter , "=====")
    Jacobstk  = []
    Bstk      = []
    c0column  = []
    
    tau_klassic_stk = []
    tau_line_stk    = []
    tau_obs_stk     = []
    
#    for ipxp , (Xbato , ObsASM , pxpapri , Tbato) in enumerate(zip(Xbato_lis ,
#                                                         ObsASM_lis ,
#                                                         PXPapri_lis,
#                                                         Tbato_lis)):
    for ipxp , (Xbato , ObsASM , pxpapri ) in enumerate(zip(Xbato_lis ,
                                                            ObsASM_lis ,
                                                            PXPapri_lis)):

        lines_stk = []
        c_v_stk   = []
    
#        for xbato , oasm , t in zip(Xbato , ObsASM , Tbato):
        for xbato_v , oasm  in zip(Xbato , ObsASM):

            # xbato xpxp kc c0 tau t
            if not with_monoc0:
                if not with_tempo_evol:
                    argstup = (xbato_v[0],xbato_v[1],xbato_v[2],
                               pxpapri[0],pxpapri[1],pxpapri[2],
                               0,c0_lis0[ipxp],oasm,0,0)
                else:
                    argstup = (xbato_v[0],xbato_v[1],xbato_v[2],
                               pxpapri[0],pxpapri[1],pxpapri[2],
                               kc_lis[ipxp],c0_lis[ipxp],oasm,t,Tbato[0])
            else:
                ang = np.arccos((pxpapri[2] - xbato_v[2]) / \
                      np.linalg.norm(xbato_v - pxpapri)) * (180./np.pi)
                c_precalc = ssp.SSP_mean(Zssp,Cssp,ang,pxpapri[2])
                
                # PTDMODIF A  
                c_v = c_precalc
                c_v = c0_lis[0]
                # PTDMODIF A
                c_v_stk.append(c_v)
                
                argstup = (xbato_v[0],xbato_v[1],xbato_v[2],
                           pxpapri[0],pxpapri[1],pxpapri[2],
                           0,c_v,oasm,0,0) 
                           
            #a titre informatif 
            if 1:
                tau_line        = fctobs_lambda(xbato_v[0],xbato_v[1],xbato_v[2],
                               pxpapri[0],pxpapri[1],pxpapri[2],
                               0,c_v,oasm,0,0)
                try:      
                    _,_,tau_klassic = rt.raytrace_seek(np.array([xbato_v[0],xbato_v[1],xbato_v[2]]),
                                   np.array([pxpapri[0],pxpapri[1],pxpapri[2]]),
                                   Zssp,Cssp,fulloutput=0,verbose=0)
                except:
                    tau_klassic = np.nan
                               
                tau_line_stk.append(tau_line)
                tau_klassic_stk.append(tau_klassic)
                tau_obs_stk.append(oasm)
                               
            if 0:              
                xbatotst = np.array([0,0,0])
                xpxptst  = np.array([0,0,4000])
                
                aaa = fctobs_lambda(xbatotst[0],xbatotst[1],xbatotst[2],
                                                 xpxptst[0],xpxptst[1],xpxptst[2],
                               0,c_v,oasm,0,0)
                _,_,bbb = rt.raytrace_seek(xbatotst,xpxptst,
                               Zssp,Cssp,fulloutput=0,verbose=0)
                               
                ccc = np.linalg.norm(xbatotst - xpxptst) / \
                (ssp.SSP_mean(Zssp,Cssp,xbatotst[-1],xpxptst[-1]) )
                               
                dtauu = aaa - bbb
                           
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
            if with_monoc0:
                line = [ xpxpdiff , ypxpdiff , zpxpdiff ]
                c0column.append( c0diff )
            
            lines_stk.append(line)
        
        Jacobstk.append(np.vstack(lines_stk))
        
    A    = scipy.linalg.block_diag(*Jacobstk)
    

    if with_c0_estim:
        if with_monoc0:
            A = np.column_stack((A,np.array(c0column)))
            

    # AJOUR DE LA CONTRAINTE SUR c0
    if with_c0_estim:
        nobstmp  = A.shape[0]
        nincotmp = A.shape[1]
        Ac0 = np.zeros(nincotmp)
        Ac0[-1] = 1
        A = np.vstack((A,Ac0))
        Bstk.append(c0_lis0[0] - c_v)
        # et du coup faut une matrice de poids
        _,_,P = geok.weight_mat([10**-6,10**-4],[nobstmp,1],sparsediag=1)
    else:
        nobstmp  = A.shape[0]
        _,_,P = geok.weight_mat([1*10**-6],[nobstmp],sparsediag=1)
        
    B     = np.array(Bstk)
    
    Bsprs = sprs.csc_matrix(B)
    Asprs = sprs.csc_matrix(A)
    Psprs = sprs.csc_matrix(P)
        
    Ninv = sprs.linalg.inv(sprs.csc_matrix((Asprs.T).dot(Psprs).dot(Asprs)))
    
    
    dX   = Ninv.dot(Asprs.T).dot(Psprs).dot(Bsprs.T)
    dX   = np.squeeze(dX.toarray())
    Ninv = Ninv.toarray()

    print('interim dX' , dX)
    print('interim sum dX' , np.sum(dX))

    if not with_monoc0:
        if not with_tempo_evol:    
            dXreshape = dX.reshape((Npxp,4))
        else:
            dXreshape = dX.reshape((Npxp,5))
        if with_c0_estim:
            c0_lis    = np.array(c0_lis) + dXreshape[:,3]
    else:
        if not with_c0_estim:
            dXreshape = dX.reshape((Npxp,3))
        else:
            dXreshape = dX[:-1].reshape((Npxp,3))
            c0_lis    = np.array(c0_lis) + dX[-1]

    PXPapri_lis = list(np.array(np.array(PXPapri_lis) + dXreshape[:,0:3]))
    if with_tempo_evol:
        kc_lis = np.array(kc_lis) + dXreshape[:,4]

    print("interm PXP" , PXPapri_lis)
    print("interm barycenter" , geok.barycenter(PXPapri_lis))    
    print("diff true barycenter    " , geok.barycenter(PXPapri_lis) - geok.barycenter(PXP_lis))
    print("diff true bary 3D (norm)" , np.linalg.norm(geok.barycenter(PXPapri_lis) - geok.barycenter(PXP_lis)))
    print("diff true bary 2D (norm)" , np.linalg.norm(geok.barycenter(PXPapri_lis)[0:2] - geok.barycenter(PXP_lis)[0:2]))
    


V = B - A.dot(dX)

print("final PXP               " , PXPapri_lis)
print("final barycenter        " , geok.barycenter(PXPapri_lis))
print("diff true barycenter    " , geok.barycenter(PXPapri_lis) - geok.barycenter(PXP_lis))
print("diff true bary 3D (norm)" , np.linalg.norm(geok.barycenter(PXPapri_lis) - geok.barycenter(PXP_lis)))
print("diff true bary 2D (norm)" , np.linalg.norm(geok.barycenter(PXPapri_lis)[0:2] - geok.barycenter(PXP_lis)[0:2]))
if with_c0_estim:
    print("c0" , dX[-1] , c0_lis[0] , c0_lis0[0])


aaaaa = (V.T).dot(P.toarray()).dot(V) / (A.shape[0] - A.shape[1])

fuv = geok.fuv_calc(V,A,P,1)
Sigma = np.sqrt((np.diag(Ninv) * fuv))

print("fuv    : " , fuv)
print("sigma  : " , Sigma)

plt.figure()
plt.hist(V[:-1],100)
plt.figure()
plt.plot(c_v_stk)

tau_klassic_stk = np.array(tau_klassic_stk)
tau_line_stk    = np.array(tau_line_stk)
tau_obs_stk     = np.array(tau_obs_stk)

np.nanmax(tau_klassic_stk - tau_obs_stk)
np.nanmax(tau_line_stk - tau_obs_stk)


dtau = tau_klassic_stk - tau_line_stk

plt.figure()
plt.plot(dtau)
plt.title(str(polydeg))
