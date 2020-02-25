#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 16:44:41 2017

@author: psakicki
"""

from megalib import *
import pandas as pds
from scipy.interpolate import interp1d

det = "/home/psakicki/THESE/1705_CANOPUS/DATA_ACOUSTIC_DECODED/donnees 14 mai/Detection14mai.txt"
nav = "/home/psakicki/THESE/1705_CANOPUS/DATA_ACOUSTIC_DECODED/donnees 14 mai/Navigation14mai.txt"
pos = "/home/psakicki/THESE/1705_CANOPUS/DATA_ACOUSTIC_DECODED/donnees 14 mai/Position14mai.txt"
nav = "/home/psakicki/THESE/1705_CANOPUS/DATA_ACOUSTIC_DECODED_2/transfer.ixblue.com_download_20170622-1710/20170514_INSAcoustic50Hz.txt"

det = "/home/psakicki/THESE/1705_CANOPUS/DATA_ACOUSTIC_DECODED/donnees 16 mai/Detection16mai.txt"
nav = "/home/psakicki/THESE/1705_CANOPUS/DATA_ACOUSTIC_DECODED/donnees 16 mai/Navigation16mai.txt"
pos = "/home/psakicki/THESE/1705_CANOPUS/DATA_ACOUSTIC_DECODED/donnees 16 mai/Position16mai.txt"
nav = "/home/psakicki/THESE/1705_CANOPUS/DATA_ACOUSTIC_DECODED_2/transfer.ixblue.com_download_20170622-1710/2017051516_INSacoustic50Hz.txt"

det = "/home/psakicki/THESE/1705_CANOPUS/DATA_ACOUSTIC_DECODED/donnees 14 mai/Detection14mai.txt"
nav = "/home/psakicki/THESE/1705_CANOPUS/DATA_ACOUSTIC_DECODED/donnees 14 mai/Navigation14mai.txt"
pos = "/home/psakicki/THESE/1705_CANOPUS/DATA_ACOUSTIC_DECODED/donnees 14 mai/Position14mai.txt"
nav = "/home/psakicki/THESE/1705_CANOPUS/DATA_ACOUSTIC_DECODED_2/transfer.ixblue.com_download_20170622-1710/20170514_INSAcoustic50Hz.txt"

det = "/home/psakicki/THESE/1705_CANOPUS/DATA_ACOUSTIC_DECODED/donnees 15 mai/Detection15mai.txt"
nav = "/home/psakicki/THESE/1705_CANOPUS/DATA_ACOUSTIC_DECODED/donnees 15 mai/Navigations15mai.txt"
pos = "/home/psakicki/THESE/1705_CANOPUS/DATA_ACOUSTIC_DECODED/donnees 15 mai/Positions15mai.txt"
nav = "/home/psakicki/THESE/1705_CANOPUS/DATA_ACOUSTIC_DECODED_2/transfer.ixblue.com_download_20170622-1710/2017051516_INSacoustic50Hz.txt"

path = det

###### LOADING DATA
if 1:
    # READ AS TABLE
    P = pds.read_table(pos,header=15)
    N = pds.read_table(nav,header=26)
    D = pds.read_table(det,header=24)

if 1:
    # SAVE TABLE AS PICKLE
    gf.pickle_saver(P,full_path=pos + '.pick')
    gf.pickle_saver(N,full_path=nav + '.pick')
    gf.pickle_saver(D,full_path=det + '.pick')
    
if 0:
    # RELOAD THE PICKLES (FASTER)
    P = gf.pickle_loader(pos + '.pick')
    N = gf.pickle_loader(nav + '.pick')
    D = gf.pickle_loader(det + '.pick')    


# CONVERSION DES EPOCHS DE LA NAV EN DATETIME & EN POSIX
Tnav = [geok.string_date2dt(str(d) + ' ' + str(t)) for d,t in zip(N['#date'] , N['time'])]
Tnav_posix = geok.dt2posix(Tnav)

# CONVERSION DES COORDS GEO EN XYZ
X,Y,Z = geok.GEO2XYZ(N['latitude'] , N['longitude'] , N['altitude'] , 'deg')

# DES PLOTS DE LA TRAJECTOIRE POUR LA ROUTE 
fig,ax = plt.subplots()
ax.plot(N['longitude']  , N['latitude'] )

ax.plot(6.864637854495319, 43.199540947296811,'^')
ax.plot(6.8743234521733232,43.206186361307864, '^')
ax.plot(6.883538876197024, 43.199717946703331,'^')

ax.plot(6.8745318,43.1930143,'v')
ax.plot(6.883237,43.1995814,'v')
ax.plot(6.8746002,43.2062435,'v')


P = np.column_stack((X,Y,Z,Tnav))

#ts = gcls.TimeSeriePoint()

# RECHERCHE DU POINT MEDIAN DU CHANTIER
x0,y0,z0 = np.nanmedian(X),np.nanmedian(Y),np.nanmedian(Z)
x0,y0,z0 = 4623265.68527  , 557450.9997   , 4343743.20306

# PASSAGE DES COORDS EN COORDONNÉES LOCALES CENTRÉES SUR LE POINT MEDIAN
E,N,U    = geok.XYZ2ENU_2(X,Y,Z,x0,y0,z0)

# CREATION D'UN INTERPOLATEUR POUR TROUVER LA POSITION EXACTE AU MOMENT DE L'EMISSION/RECEPTION
print('INFO : DEBUT generation interpolateurs')
EfT , NfT , UfT = gcls.interpolator_light(Tnav_posix,E,N,U)
print('INFO : FIN generation interpolateurs')

#for p in P:
#    pt = gcls.Point(*p,initype='XYZ')
#    ts.add_point(pt)
#    
#print 'INFO : fin du chargement de la TS'
#
#ts.ENUcalc_from_mean_posi()
#
#ts.interp_set()

pingdic = dict()
pingstddic = dict()
temidic = dict()
trecdic = dict()


t0 = geok.string_date2dt(D['#date'][0] + ' ' + D['time'][0])

tat = 300 * 10**-3

# DEFINTION  DU TAT POUR CHACUNE DES BEACON DANS UN DICO
tatdic = dict()
tatdic[11] = 150 * 10**-3
tatdic[12] = 200 * 10**-3
tatdic[13] = 250 * 10**-3
tatdic[14] = 300 * 10**-3

# definition des coords APRIORI
xyzapri       = dict()
xyzapri[11] = [ -802.69331721 , -131.81250275 , 1250.25481608]
xyzapri[13] = [  -56.38172483 , 560.12604956  , 1283.72197572]
xyzapri[14] = [  686.75297138 , -130.18414091 , 1204.48841112]

xyzapri[11] =  [ -863.32338053 ,  -77.73029269 , 1361.60493971]
xyzapri[13] =  [  -75.9346007  ,  660.65228102 , 1320.61130829]
xyzapri[14] =  [  673.19675219 ,  -58.0821559  , 1211.93784559]

# coordonnées issues de iXblue mais avec un passage dans le LSQ
xyzapri[11] = [  -78.01707448 , -854.88852201 , 1354.75153031]
xyzapri[13] = [  657.46200569 ,  -75.80416804 , 1312.0076598 ]
xyzapri[14] = [  -55.18166279 ,  666.13423397 , 1207.30815759]

# POSITION iXBLUE
PXPosiflh= dict()
PXPosiflh[11]=[43.1930143,6.8745318,-1263.785]
PXPosiflh[13]=[43.1995814,6.883237,-1282.492]
PXPosiflh[14]=[43.2062435,6.8746002,-1214.483]

#activate iXblue coordinates
if 0:
    for ibeacon , val in list(PXPosiflh.items()):
        xyzapri[ibeacon]     = np.hstack(geok.XYZ2ENU_2(*geok.GEO2XYZ(*val,angle='deg'),x0=x0,y0=y0,z0=z0))
        xyzapri[ibeacon][-1] = - xyzapri[ibeacon][-1]  


# POUR CHANQUE BALISE ON TROUVE LES POSITIONS QUI VONT BIEN AUX INSTANTS D'EMISSION RECEPTION
for idbeacon in np.unique(D['beaconId']):
    count_diff_4 = 0
    if idbeacon < 10:
        continue
    Dwork = D[idbeacon == D['beaconId']]
        
    meanstk = []
    stdstk  = []
    temistk = []
    trecstk = []
    
    pingdic[idbeacon]    = meanstk
    pingstddic[idbeacon] = stdstk
    temidic[idbeacon]    = temistk
    trecdic[idbeacon]    = trecstk
    
    print('INFO :  boucle interne pour la balise', idbeacon)
    ##BOUCLE INTERNE POUR TESTER LA VALIDITÉ DES PINGS
    for t1 , t2  ,t3 , t4 , interdat in zip(Dwork['timeH1'],Dwork['timeH2'], 
                                            Dwork['timeH3'],Dwork['timeH4'], 
                                            Dwork['interrogationDate']):
        
        ping = np.array([t1 , t2  ,t3 , t4])
        
        temi = t0 + dt.timedelta(seconds=interdat)
        
        
        ping = ping[ping != 0]
        
        if 0 : #len(ping) != 0
            count_diff_4 += 1
            continue
        
        # ON FAIT LA MOYENNE DE 4 HYDROS
        meanstk.append(np.mean(ping))
        stdstk.append(np.std(ping))
        
        temistk.append(temi)
        trecstk.append(temi + dt.timedelta(seconds=np.mean(ping)))
    
    meanstk = np.array(meanstk)
    
    T_emi_posix = geok.dt2posix(temistk,1)
    N_emi       = NfT(T_emi_posix)
    E_emi       = EfT(T_emi_posix)
    U_emi       = UfT(T_emi_posix)
    
    T_rec_human = np.column_stack([(e.year , e.month , e.day , e.hour , e.minute , e.second , e.microsecond) for e in trecstk]).T
    T_emi_human = np.column_stack([(e.year , e.month , e.day , e.hour , e.minute , e.second , e.microsecond) for e in temistk]).T

    T_rec_posix = geok.dt2posix(trecstk,1)
    N_rec = NfT(T_rec_posix)
    E_rec = EfT(T_rec_posix)
    U_rec = UfT(T_rec_posix)
    
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    U_rec = U_rec
    U_emi = U_emi
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #U_rec = U_emi
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #U_rec = np.zeros(len(U_rec))
    #U_emi = np.zeros(len(U_emi))
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    tatarr =  np.array([tatdic[idbeacon]] * len(meanstk))
    
    tau_m_tat = np.array(meanstk) - tatarr
    
    
    tau_m_tat_clean , bool_mad = geok.outiler_mad(tau_m_tat)
    
    print(np.sum(np.logical_not(bool_mad)) ," outliers removed ")
    
    # On GENERE LA "MATRICE"  i.e. le tableau des données pour le O file
    tup = ( T_emi_posix[bool_mad] ,
            E_emi[bool_mad] ,
            N_emi[bool_mad] ,
            - U_emi[bool_mad],
            T_emi_human[bool_mad] ,
            T_rec_posix[bool_mad] ,
            E_rec[bool_mad] ,
            N_rec[bool_mad] ,
            - U_rec[bool_mad] , 
            T_rec_human[bool_mad] ,
            meanstk[bool_mad] ,
            tatarr[bool_mad] ,
            tau_m_tat[bool_mad] )
    
    M = np.column_stack(tup)
    
    prefix = os.path.basename(path).split('.')[0]
    
    # ON ECRIT LE O FILE
    savedir =  '/home/psakicki/THESE/1705_CANOPUS/O_files/' + prefix
    gf.create_dir(savedir)
    
    print('INFO : ' , count_diff_4 , 'pings removed bc not 4 valid ping') 
    print("iNFO : will be saved in " + savedir)
    print(prefix)
    
    name = genefun.join_improved('.',prefix,'PXP'+ str(int(idbeacon)),'O','dat')

    if idbeacon < 10:
        idbeacon2 = idbeacon + 10.
    else:
        idbeacon2 = idbeacon
        
    com = "pxp_coords : " + str(np.array(xyzapri[idbeacon2]))
    np.savetxt(savedir + '/' + name , M , header = com)
    
    plt.figure()
    plt.plot(trecstk,U_rec,'+')
    plt.plot(temistk,U_emi,'x')
    
pathctd = '/home/psakicki/THESE/1705_CANOPUS/CTD/0514_0800Z/SBE19plus_01906689_2017_05_14_0006.cnv.txt'
ZC = np.loadtxt(pathctd)

#nameZ = genefun.join_improved('.',prefix,'PXP'+ str(int(idbeacon)),'Z','dat')
#np.savetxt(savedir + '/' + nameZ , ZC[:,0])
#nameC = genefun.join_improved('.',prefix,'PXP'+ str(int(idbeacon)),'C','dat')
#np.savetxt(savedir + '/' + nameC , ZC[:,1])

plt.figure()

#UUU = ts.UfT(ts.to_list(specific_output=3))

matplotlib.rcParams.update({'font.size': 20})

plt.clf()
plt.figure()
plt.plot(Tnav,U,'.')
plt.plot(trecstk,UfT(geok.dt2posix(trecstk)),'x')
plt.plot(temistk,UfT(geok.dt2posix(temistk)),'+')
plt.gca().set_ylabel('Up component (m)')

plt.figure()
plt.plot(Tnav,np.array(E),'.')
plt.plot(trecstk,EfT(geok.dt2posix(trecstk)),'x')
plt.plot(temistk,EfT(geok.dt2posix(temistk)),'+')
plt.gca().set_ylabel('East component (m)')

plt.figure()
plt.plot(Tnav,np.array(N),'.')
plt.plot(trecstk,NfT(geok.dt2posix(trecstk)),'x')
plt.plot(temistk,NfT(geok.dt2posix(temistk)),'+')
plt.gca().set_ylabel('North component (m)')

#plt.figure()
#plt.plot(geok.posix2dt(ts.to_list(specific_output=3)),UUU,'+')


