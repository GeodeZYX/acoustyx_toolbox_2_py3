# -*- coding: utf-8 -*-
"""
Created on Fri May  6 18:15:39 2016

@author: psakicki
"""

import fabrik_8_fct_mk7 as fab
from   megalib import *


################ TOOLBOX PATH
toolbox_path = "/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/"
toolbox_path = "/dsk/ggsp_pf/PLAYGROUND/psakicki/scripts_PS/zyx_TOOLBOXs/acoustyx_toolbox_2_py3"
toolbox_path = "/home/psakicki/CODES/acoustyx_toolbox_2_py3/"


################ EXEMPLE RUN
if 0:
    exp_prefix   = 'EXEMPLE_RUN'

    # STEP 1 : here we define the changeable parameters
    #          values for each parameters are stored in lists
    Ldesign = ['grand-carree','triangle']
    Ltraj   = ['station','derive']
    Lnbobs  = [5000,10000]
    Ldate = [dt.datetime(2009, 10, 0o2, 00, 00)]

    # here is made the cartesian product of the parameters lists
    # we obtain a list of combinations
    ITER = list(itertools.product(Ltraj,Lnbobs,Ldate,Ldesign))

    # We run a loop for each combination
    print(len(ITER),' experiences en batch')
    for i,itup in enumerate(ITER):
        print(i , '/' , len(ITER) ,' experiences en batch')

        # STEP 2 : here the parameters are recovered in temporary variables
        traj  = itup[0]
        nbobs = itup[1]
        dat   = itup[2]
        design= itup[3]

        dico = dict()

        # STEP 3 : temporary variables are assigned to the official variables
        #          name as in the fabrik_fct function/script,
        #          trought a dictionary
        dico['prm_design']   =  design

        dico['prm_zones_bound_input']  = [1000]
        dico['prm_sigma_zones_input']  = [10**-4,10**-5]

        dico['prm_with_temporal']   = 1
        dico['prm_tempor_start']    = dat
        dico['tempor_ctd_time']     = dat
        dico['prm_tempor_ctd_time'] = dat

        dico['prm_null_water_column_noise'] = 1
        dico['prm_coef_sigma_zones']        = 0

        dico['prm_with_ssf']           = 0
        dico['prm_with_pseudo_ssf']    = 0

        dico['prm_sigma_x'] = 0.0
        dico['prm_sigma_y'] = 0.0
        dico['prm_sigma_z'] = 0.0

        if traj == 'droite':
            dico['prm_drt_nb_obs']    =  nbobs/3
        if traj == 'station':
            dico['prm_sta_nb_obs']    =  nbobs
        else:
            dico['prm_drv_nb_obs']    =  nbobs

        dico["prm_add_offset"]    =  False
        dico["imporved_noising"]  =  1
        dico["prm_traject"]       =  traj
        dico['prm_drv_R']         =  10
        dico['prm_vit']           =  1

        dico['procs']             = 7

        # the midfix is a customisable name
        midfix = 'date' + str(dat.strftime('%Y%m%d')) + '_' + design

        # the experience with custom parameters is launched
        fab.fabrik_fct(exp_prefix,midfix,dico)



################ OPERATIONAL RUNs

if 0:
    exp_prefix  = 'TRAJECT_NOISE_mk1b'
    L = [0,0.01,0.05,0.1,0.5]
    ITER = list(itertools.product(L,L,L))


    ## LA BOUCLE
    print(len(ITER),' experiences en batch')
    for i,itup in enumerate(ITER):
        print(i , '/' , len(ITER) ,' experiences en batch')

        dico = dict()
        dico['prm_sigma_x'] = itup[0]
        dico['prm_sigma_y'] = itup[1]
        dico['prm_sigma_z'] = itup[2]
        dico['prm_nb_obs']  = 100

        midfix = gf.join_improved('_','traject_noise',*itup)
        fab.fabrik_fct(exp_prefix,midfix,dico)

# ================================================================

if 0:
    exp_prefix  = 'OFFSETS_mk1b'
    prm_add_offset = True
    L = [0,0.01,0.05,0.1,0.5]
    ITER = list(itertools.product(L,L,L))

    ## LA BOUCLE
    print(len(ITER),' experiences en batch')
    for i,itup in enumerate(ITER):
        print(i , '/' , len(ITER) ,' experiences en batch')

        dico = dict()
        dico['prm_offset'] = np.array(itup)
        dico['prm_nb_obs']          = 100
        dico["prm_add_offset"]  =  prm_add_offset

        midfix = gf.join_improved('_','offsets',*itup)
        fab.fabrik_fct(exp_prefix,midfix,dico)

# ================================================================

if 0:
    exp_prefix  = 'BUOYstyle_mk3'

    nbdays                = 365
    L_prm_tempor_start    = [dt.datetime(2012, 1, 1, 12, 00) + dt.timedelta(days=i) for i in range(nbdays)]

    intern_random_seed_x = 789
    intern_random_seed_y = 456
    sigma_posi_x  = 50
    sigma_posi_y  = 50

    Rx = np.random.RandomState(intern_random_seed_x)
    Ry = np.random.RandomState(intern_random_seed_y)

    L_prm_x0 = 10000 + Rx.randn(365) * sigma_posi_x
    L_prm_y0 = 10000 + Ry.randn(365) * sigma_posi_y

    prm_tempor_ctd_time = L_prm_tempor_start[0]

    ITER = list(zip(L_prm_tempor_start,L_prm_x0,L_prm_y0))

    error_counter = 0

    # LA BOUCLE
    print(len(ITER),' experiences en batch')
    for i,itup in enumerate(ITER):
        print(i , '/' , len(ITER) ,' experiences en batch')

        dico = dict()
        dico['prm_with_temporal']   = True
        dico['prm_tempor_start']    = itup[0]
        dico['prm_x0']              = itup[1]
        dico['prm_y0']              = itup[2]
        dico['prm_tempor_ctd_time'] = prm_tempor_ctd_time
        dico['prm_tempor_len']      = 3600*4
        dico['prm_nb_obs']          = 100
        dico['prm_R']               = 500
        dico['prm_K_derive']        = (dico['prm_nb_obs'] + i) * 1
        dico['prm_K_xyz']           = (dico['prm_nb_obs'] + i) * 42


        dico['procs']               = 7

        midfix = gf.join_improved('_','BUOYstyle',itup[0].date(),
                                  *[np.round(e) for e in itup[1:]])

        try:
            fab.fabrik_fct(exp_prefix,midfix,dico)
        except Exception as e:
            print("ERR : for " ,  itup)
            error_counter += 1
            print("error_counter" ,  error_counter)
            continue

# ================================================================

if 0:
    exp_prefix  = 'BUOYstyle_mk4'

    nbdays                = 365
    L_prm_tempor_start    = [dt.datetime(2012, 1, 1, 12, 00) + dt.timedelta(days=i) for i in range(nbdays)]

    intern_random_seed_x = 789
    intern_random_seed_y = 456
    sigma_posi_x  = 50
    sigma_posi_y  = 50

    Rx = np.random.RandomState(intern_random_seed_x)
    Ry = np.random.RandomState(intern_random_seed_y)

    L_prm_x0 = 10000 + Rx.randn(365) * sigma_posi_x
    L_prm_y0 = 10000 + Ry.randn(365) * sigma_posi_y

    prm_tempor_ctd_time = L_prm_tempor_start[0]

    ITER = list(zip(L_prm_tempor_start,L_prm_x0,L_prm_y0))

    error_counter = 0

    # LA BOUCLE
    print(len(ITER),' experiences en batch')
    for i,itup in enumerate(ITER):
        print(i , '/' , len(ITER) ,' experiences en batch')

        dico = dict()
        dico['prm_with_temporal']   = True
        dico['prm_tempor_start']    = itup[0]
        dico['prm_x0']              = itup[1]
        dico['prm_y0']              = itup[2]
        dico['prm_tempor_ctd_time'] = itup[0] #prm_tempor_ctd_time
        dico['prm_tempor_len']      = 3600*4
        dico['prm_nb_obs']          = 100
        dico['prm_R']               = 500
        dico['prm_K_derive']        = (dico['prm_nb_obs'] + i) * 1
        dico['prm_K_t_zone']        = (dico['prm_nb_obs'] + i) * 42


        dico['procs']               = 7

        midfix = gf.join_improved('_','BUOYstyle',itup[0].date(),
                                  *[np.round(e) for e in itup[1:]])

        try:
            fab.fabrik_fct(exp_prefix,midfix,dico)
        except Exception as e:
            print("ERR : for " ,  itup)
            error_counter += 1
            print("error_counter" ,  error_counter)
            continue

# ================================================================

if 0:
    exp_prefix   = 'ULTRAKLASSIK_mk9'
    midfix_str   = 'obs'
    L = [100,500,1000,2000,5000,10000,20000]

    ITER = list(L)

    ## LA BOUCLE
    print(len(ITER),' experiences en batch')
    for i,itup in enumerate(ITER):
        print(i , '/' , len(ITER) ,' experiences en batch')

        dico = dict()
        dico['prm_with_temporal'] =  False
        dico['prm_nb_obs']        =  itup
        dico["prm_add_offset"]    =  False
        dico["imporved_noising"]  =  True
        dico["prm_traject"]       =  'derive'
        dico['prm_R'] = 10

        dico['procs']             = 7

        midfix = gf.join_improved('_',midfix_str,itup)
        fab.fabrik_fct(exp_prefix,midfix,dico)

# ================================================================

if 0:
    exp_prefix   = 'STRAIGHT_TRAJ_OFFSET_mk1'
    midfix_str   = 'angle'

    L = np.arange(0,91,10)
    L = (0,90)
    L = np.arange(10,81,10)

    ITER = list(L)

    ## LA BOUCLE
    print(len(ITER),' experiences en batch')
    for i,itup in enumerate(ITER):
        print(i , '/' , len(ITER) ,' experiences en batch')

        dico = dict()
        dico['prm_with_temporal'] = False
        dico['prm_nb_obs']        = 100
        dico["prm_add_offset"]    =  True
        dico["prm_offset"]        = [0.1,0,0]
        dico["imporved_noising"]  =  True
        dico["prm_null_water_column_noise"]  =  True
        dico["prm_traject"]       =  'droite'
        dico["prm_angle"]   =  itup
        dico["prm_sigma_x"] =  0.0
        dico["prm_sigma_y"] =  0.0
        dico["prm_sigma_z"] =  0.0

        dico['procs']               = 7

        midfix = gf.join_improved('_',midfix_str,dico["prm_angle"])
        fab.fabrik_fct(exp_prefix,midfix,dico)

# ================================================================

if 0:
    exp_prefix  = 'KLASSIK_DECAL_CENTER_mk1'
    midfix_str   = 'obs'
    L    = np.array([10000,10500,11000,12000])
    ITER = list(itertools.product(L,L))

    print(ITER)

    ## LA BOUCLE
    print(len(list(ITER)),' experiences en batch')
    for i,itup in enumerate(ITER):
        print(i , '/' , len(ITER) ,' experiences en batch')

        dico = dict()
        dico['prm_with_temporal'] = False
        dico['prm_nb_obs']        = 1000
        dico['prm_x0']        = itup[0]
        dico['prm_y0']        = itup[1]
        dico['prm_xcenter']        = itup[0]
        dico['prm_ycenter']        = itup[1]
        dico["prm_add_offset"]    =  False
        dico["imporved_noising"]  =  True
        dico["prm_traject"]       =  'derive'

        dico['procs']               = 7

        midfix = gf.join_improved('_',midfix_str,itup)
        fab.fabrik_fct(exp_prefix,midfix,dico)

# ================================================================

if 0:
    exp_prefix   = 'BATCH_RC_mkADEL'
    midfix_str   = 'obs'
    Ltraj = ['derive']
    LR    = [10,50,100,1000]


    ITER = list(itertools.product(Ltraj,LR))

    ## LA BOUCLE
    print(len(ITER),' experiences en batch')
    for i,itup in enumerate(ITER):
        print(i , '/' , len(ITER) ,' experiences en batch')

        traj = itup[0]
        R    = itup[1]

        if traj == 'droite' and R != 10:
            continue

        dico = dict()
        dico['prm_with_temporal'] =  False
        if traj == 'droite':
            dico['prm_nb_obs']    =  1000/3
        else:
            dico['prm_nb_obs']    =  1000

        dico["prm_add_offset"]    =  False
        dico["imporved_noising"]  =  True
        dico["prm_traject"]       =  traj
        dico['prm_R']             = R

        dico['procs']             = 7

        midfix = '_'
        fab.fabrik_fct(exp_prefix,midfix,dico)


# ================================================================

if 0:
    exp_prefix   = 'BATCH_RC_mk5b_authentik_dekal_zonal_noise_wo_10koef'
    midfix_str   = 'obs'
    Ltraj = ['derive']
    LR     = [1000,100,10]
    LR     = [50,500]
    Lnbobs = [1000,5000]

    ITER = list(itertools.product(Ltraj,LR,Lnbobs))

    ## LA BOUCLE
    print(len(ITER),' experiences en batch')
    for i,itup in enumerate(ITER):
        print(i , '/' , len(ITER) ,' experiences en batch')

        traj = itup[0]
        R    = itup[1]

        if traj == 'droite' and R != 10:
            continue

        dico = dict()
        dico['prm_with_temporal'] =  0
        if traj == 'droite':
            dico['prm_nb_obs']    =  itup[2]/3
        else:
            dico['prm_nb_obs']    =  itup[2]

        dico["prm_add_offset"]    =  False
        dico["imporved_noising"]  =  1
        dico["prm_traject"]       =  traj
        dico['prm_R']             =  R

        dico['procs']             = 6

        midfix = '_'
        fab.fabrik_fct(exp_prefix,midfix,dico)

# ================================================================

if 0:
    exp_prefix   = 'BATCH_RC_mk2_OPERASCENAR_tempor_AVEC_pseudograd_10p6'
    Ltraj = ['derive']
    Lnbobs = [100,500,1000,5000,10000]
    Ldate = [dt.datetime(2007, 4, 24, 00, 00),
             dt.datetime(2010, 2, 25, 00, 00),
             dt.datetime(2009, 6, 23, 00, 00)]
    ITER = list(itertools.product(Ltraj,Lnbobs,Ldate))

    ## LA BOUCLE
    print(len(ITER),' experiences en batch')
    for i,itup in enumerate(ITER):
        print(i , '/' , len(ITER) ,' experiences en batch')

        traj  = itup[0]
        nbobs = itup[1]
        dat   = itup[2]

        dico = dict()
        dico['prm_with_temporal']   =  1
        dico['prm_tempor_start']    = dat
        dico['tempor_ctd_time']     = dat
        dico['prm_tempor_ctd_time'] = dat

        if traj == 'droite':
            dico['prm_drt_nb_obs']    =  nbobs/3
        else:
            dico['prm_drv_nb_obs']    =  nbobs

        dico["prm_add_offset"]    =  False
        dico["imporved_noising"]  =  1
        dico["prm_traject"]       =  traj
        dico['prm_drv_R']             =  50

        dico['procs']             = 3

        midfix = 'date' + str(dat.strftime('%Y%m%d'))
        fab.fabrik_fct(exp_prefix,midfix,dico)

# ================================================================

if 0:
    exp_prefix   = 'FINALb_VariTempor_basic'
    exp_prefix   = 'FINALb_VariTemporCORIGEDb_basic'
    exp_prefix   = 'FINALb_NO-VariTemporCORIGEDb_basic_REBOOT1803'

    Ldesign = ['grand-carree',"petit-carree","triangle","grand-trigrangle"]
    Ltraj   = ['derive']
    Ltraj   = ["droite"]
    Ltraj   = ['derive',"droite"]
    Ldesign = ['grand-carree',"petit-carree"]
    Ldesign = ['grand-carree']

    Lnbobs  = [10000]
    Lnbobs  = [100,500,1000,5000,10000]

    Ldate = [dt.datetime(2007, 4, 24, 00, 00),
             dt.datetime(2010, 2, 25, 00, 00),
             dt.datetime(2009, 10, 0o2, 00, 00)]

    Ldate = [dt.datetime(2009, 10, 0o2, 00, 00)]

    ITER = list(itertools.product(Ltraj,Lnbobs,Ldate,Ldesign))

    ## LA BOUCLE
    print(len(ITER),' experiences en batch')
    for i,itup in enumerate(ITER):
        print(i , '/' , len(ITER) ,' experiences en batch')

        traj  = itup[0]
        nbobs = itup[1]
        dat   = itup[2]
        design= itup[3]

        dico = dict()

        dico['prm_design']   =  design

        dico['prm_zones_bound_input']  = [1000]
        dico['prm_sigma_zones_input']  = [10**-4,10**-5]

        dico['prm_with_temporal']   =  1
        dico['prm_tempor_start']    = dat
        dico['tempor_ctd_time']     = dat
        dico['prm_tempor_ctd_time'] = dat

        dico['prm_null_water_column_noise'] = 1
        dico['prm_coef_sigma_zones']        = 0

        dico['prm_with_ssf']           = 0
        dico['prm_with_pseudo_ssf']    = 0

        dico['prm_sigma_x'] = 0.03
        dico['prm_sigma_y'] = 0.03
        dico['prm_sigma_z'] = 0.05

        if traj == 'droite':
            dico['prm_drt_nb_obs']    =  int(nbobs/3)
        else:
            dico['prm_drv_nb_obs']    =  nbobs

        dico["prm_add_offset"]    =  False
        dico["imporved_noising"]  =  1
        dico["prm_traject"]       =  traj
        dico['prm_drv_R']         =  50

        dico['procs']             = 36

        midfix = 'date' + str(dat.strftime('%Y%m%d')) + '_' + design

        fab.fabrik_fct(exp_prefix,midfix,dico)

if 0:
    exp_prefix   = 'FINALb_VariTempor_basic_dekal'

    Ldesign = ['grand-carree',"petit-carree"]
    Ltraj   = ['derive']
    Ltraj   = ['derive',"droite"]

    Lnbobs  = [100,1000,10000]
    Ldate = [dt.datetime(2007, 4, 24, 00, 00),
             dt.datetime(2010, 2, 25, 00, 00),
             dt.datetime(2009, 6, 23, 00, 00)]
    ITER = list(itertools.product(Ltraj,Lnbobs,Ldate,Ldesign))

    ## LA BOUCLE
    print(len(ITER),' experiences en batch')
    for i,itup in enumerate(ITER):
        print(i , '/' , len(ITER) ,' experiences en batch')

        traj  = itup[0]
        nbobs = itup[1]
        dat   = itup[2]
        design= itup[3]

        dico = dict()

        dico['prm_design']   =  design

        dico['prm_zones_bound_input']  = [1000]
        dico['prm_sigma_zones_input']  = [10**-4,10**-5]

        dico['prm_with_temporal']   =  1
        dico['prm_tempor_start']    = dat
        dico['tempor_ctd_time']     = dat
        dico['prm_tempor_ctd_time'] = dat

        dico['prm_null_water_column_noise'] = 1
        dico['prm_coef_sigma_zones']        = 0

        dico['prm_with_ssf']           = 0
        dico['prm_with_pseudo_ssf']    = 0

        dico['prm_sigma_x'] = 0.03
        dico['prm_sigma_y'] = 0.03
        dico['prm_sigma_z'] = 0.05

        if traj == 'droite':
            dico['prm_drt_nb_obs']    =  nbobs/3
            dico['prm_drt_xcenter']   = 11000
        else:
            dico['prm_drv_nb_obs']    =  nbobs
            dico['prm_drv_xcenter']   = 11000
            dico['prm_drv_x0']        = 11000

        dico["prm_add_offset"]    =  False
        dico["imporved_noising"]  =  1
        dico["prm_traject"]       =  traj
        dico['prm_drv_R']         =  50

        dico['procs']             = 6

        midfix = 'date' + str(dat.strftime('%Y%m%d')) + '_' + design

        fab.fabrik_fct(exp_prefix,midfix,dico)

if 0:
    exp_prefix   = 'FINALb_VariTempor_ss-bruit-surface'

    Ltraj   = ['derive',"droite"]
    Ldesign = ['grand-carree',"petit-carree",
               "triangle","grand-trigrangle"]

    Lnbobs  = [100,1000,10000]
    Ldate = [dt.datetime(2007, 4, 24, 00, 00),
             dt.datetime(2010, 2, 25, 00, 00),
             dt.datetime(2009, 6, 23, 00, 00)]
    ITER = list(itertools.product(Ltraj,Lnbobs,Ldate,Ldesign))

    ## LA BOUCLE
    print(len(ITER),' experiences en batch')
    for i,itup in enumerate(ITER):
        print(i , '/' , len(ITER) ,' experiences en batch')

        traj  = itup[0]
        nbobs = itup[1]
        dat   = itup[2]
        design= itup[3]

        dico = dict()

        dico['prm_design']   =  design

        dico['prm_zones_bound_input']  = [1000]
        dico['prm_sigma_zones_input']  = [10**-4,10**-5]

        dico['prm_with_temporal']   =  1
        dico['prm_tempor_start']    = dat
        dico['tempor_ctd_time']     = dat
        dico['prm_tempor_ctd_time'] = dat

        dico['prm_null_water_column_noise'] = 1
        dico['prm_coef_sigma_zones']        = 0

        dico['prm_with_ssf']           = 0
        dico['prm_with_pseudo_ssf']    = 0

        dico['prm_sigma_x'] = 0.0
        dico['prm_sigma_y'] = 0.0
        dico['prm_sigma_z'] = 0.0

        if traj == 'droite':
            dico['prm_drt_nb_obs']    =  nbobs/3
        else:
            dico['prm_drv_nb_obs']    =  nbobs

        dico["prm_add_offset"]    =  False
        dico["imporved_noising"]  =  1
        dico["prm_traject"]       =  traj
        dico['prm_drv_R']         =  50

        dico['procs']             = 7

        midfix = 'date' + str(dat.strftime('%Y%m%d')) + '_' + design

        fab.fabrik_fct(exp_prefix,midfix,dico)

if 0:
    exp_prefix   = 'FINALa2b_VariTempor_traj-perturbed'

    Ltraj   = ['derive',"droite"]
    Ldesign = ['grand-carree',"petit-carree","triangle","grand-trigrangle"]

    Lnbobs  = [100,1000,10000]
    Ldate = [dt.datetime(2007, 4, 24, 00, 00),
             dt.datetime(2010, 2, 25, 00, 00),
             dt.datetime(2009, 6, 23, 00, 00)]
    ITER = list(itertools.product(Ltraj,Lnbobs,Ldate,Ldesign))

    ## LA BOUCLE
    print(len(ITER),' experiences en batch')
    for i,itup in enumerate(ITER):
        print(i , '/' , len(ITER) ,' experiences en batch')

        traj  = itup[0]
        nbobs = itup[1]
        dat   = itup[2]
        design= itup[3]

        dico = dict()

        dico['prm_design']   =  design

        dico['prm_zones_bound_input']  = [1000]
        dico['prm_sigma_zones_input']  = [10**-4,10**-5]

        dico['prm_with_temporal']   =  1
        dico['prm_tempor_start']    = dat
        dico['tempor_ctd_time']     = dat
        dico['prm_tempor_ctd_time'] = dat

        dico['prm_null_water_column_noise'] = 1
        dico['prm_coef_sigma_zones']        = 0

        dico['prm_with_ssf']           = 0
        dico['prm_with_pseudo_ssf']    = 0

        dico['prm_sigma_x'] = 0.50
        dico['prm_sigma_y'] = 0.50
        dico['prm_sigma_z'] = 1.

        if traj == 'droite':
            dico['prm_drt_nb_obs']    =  nbobs/3
        else:
            dico['prm_drv_nb_obs']    =  nbobs

        dico["prm_add_offset"]    =  False
        dico["imporved_noising"]  =  1
        dico["prm_traject"]       =  traj
        dico['prm_drv_R']         =  50

        dico['procs']             = 7

        midfix = 'date' + str(dat.strftime('%Y%m%d')) + '_' + design

        fab.fabrik_fct(exp_prefix,midfix,dico)

if 0:
    exp_prefix   = 'FINALa2b_VariTempor_traj-perturbed'
    exp_prefix   = 'RB_EGU18_VariTempor_traj-perturbed'

    Ltraj   = ['derive',"droite"]
    Ldesign = ['grand-carree']

    Lnbobs  = [100,500,1000,5000,10000]
    Ldate = [dt.datetime(2007, 4, 24, 00, 00),
             dt.datetime(2010, 2, 25, 00, 00),
             dt.datetime(2009, 6, 23, 00, 00)]

    Ldate = [dt.datetime(2009, 6, 23, 00, 00)]
    LKdrift = [10,20,30]

    ITER = list(itertools.product(Ltraj,Lnbobs,Ldate,Ldesign,LKdrift))

    ## LA BOUCLE
    print(len(ITER),' experiences en batch')
    for i,itup in enumerate(ITER):
        print(i , '/' , len(ITER) ,' experiences en batch')

        traj  = itup[0]
        nbobs = itup[1]
        dat   = itup[2]
        design= itup[3]
        Kdrift = itup[4]


        dico = dict()

        dico['prm_design']   =  design

        dico['prm_zones_bound_input']  = [1000]
        dico['prm_sigma_zones_input']  = [10**-4,10**-5]

        dico['prm_with_temporal']   =  1
        dico['prm_tempor_start']    = dat
        dico['tempor_ctd_time']     = dat
        dico['prm_tempor_ctd_time'] = dat

        dico['prm_null_water_column_noise'] = 1
        dico['prm_coef_sigma_zones']        = 0

        dico['prm_with_ssf']           = 0
        dico['prm_with_pseudo_ssf']    = 0

        dico['prm_sigma_x'] = 0.50
        dico['prm_sigma_y'] = 0.50
        dico['prm_sigma_z'] = 1.

        if traj == 'droite':
            dico['prm_drt_nb_obs']    =  iter(nbobs/3)
        else:
            dico['prm_drv_nb_obs']    =  nbobs
            dico["prm_drv_K_derive"]  = Kdrift

        dico["prm_add_offset"]    =  False
        dico["imporved_noising"]  =  1
        dico["prm_traject"]       =  traj
        dico['prm_drv_R']         =  50

        dico['procs']             = 7

        midfix = 'date' + str(dat.strftime('%Y%m%d')) + '_' + design + '_' + "Kdrift" + str(Kdrift)

        fab.fabrik_fct(exp_prefix,midfix,dico)


if 0:
    exp_prefix   = 'RB_EGU18_TEMPO-D_DriftSmall_Same-dZ'
    exp_prefix   = 'RB_EGU18_TEMPO-E_DriftVerySmall_Same-dZ'

    Ltraj   = ['derive',"droite"]
    Ltraj   = ['station']
    Ltraj   = ['derive']
    Ltraj   = ['droite']

    Ldesign = ['grand-carree']
    Ldesign = ['grand-carree-no_dZ']

    Lnbobs  = [100]
    Lnbobs  = [100,500,1000,5000,10000]
    Lnbobs  = [10000]

    Ldate = [dt.datetime(2007, 4, 24, 00, 00),
             dt.datetime(2010, 2, 25, 00, 00),
             dt.datetime(2009, 6, 23, 00, 00)]

    Ldate = [dt.datetime(2009, 6, 23, 00, 00)]
    Ldate = [dt.datetime(2010, 2, 25, 00, 00)]

    LKdrift = [1111]
    LKdrift = [10,20,30]
    LKdrift = [10]

    ITER = list(itertools.product(Ltraj,Lnbobs,Ldate,Ldesign,LKdrift))

    ## LA BOUCLE
    print(len(ITER),' experiences en batch')
    for i,itup in enumerate(ITER):
        print(i , '/' , len(ITER) ,' experiences en batch')

        traj   = itup[0]
        nbobs  = itup[1]
        dat    = itup[2]
        design = itup[3]
        Kdrift = itup[4]

        dico = dict()

        dico['prm_design']             =  design

        dico['prm_zones_bound_input']  = [1000]
        dico['prm_sigma_zones_input']  = [10**-4,10**-5]

        dico['prm_with_temporal']   =  1
        dico['prm_tempor_start']    = dat
        dico['tempor_ctd_time']     = dat
        dico['prm_tempor_ctd_time'] = dat

        dico['prm_null_water_column_noise'] = 1
        dico['prm_coef_sigma_zones']        = 0

        dico['prm_with_ssf']           = 0
        dico['prm_with_pseudo_ssf']    = 0

        dico['prm_sigma_x'] = 0.
        dico['prm_sigma_y'] = 0.
        dico['prm_sigma_z'] = 0.

        if traj == 'droite':
            dico['prm_drt_nb_obs']    =  int(nbobs/3)
        elif traj == "station":
            dico['prm_sta_nb_obs']    =  nbobs
        else:
            dico['prm_drv_nb_obs']    =  nbobs
            dico["prm_drv_K_derive"]  = Kdrift

        dico["prm_add_offset"]    =  False
        dico["imporved_noising"]  =  1
        dico["prm_traject"]       =  traj
        #For a small circle
        dico['prm_drv_R']         =  .0001
        dico["prm_drv_step_size"] =  .00001

        dico['procs']             = 36

        midfix = 'date' + str(dat.strftime('%Y%m%d')) + '_' + design + '_' + "Kdrift" + str(Kdrift)

        fab.fabrik_fct(exp_prefix,midfix,dico)


if 0:
    exp_prefix   = 'RB_EGU18_KEDALnoiszA2-tempo'

    Ltraj   = ['derive']
    Ldesign = ['grand-carree']

    Lnbobs  = [100,500,1000,5000,10000]
    Lnbobs  = [100]

    Ldate = [dt.datetime(2007, 4, 24, 00, 00),
             dt.datetime(2010, 2, 25, 00, 00),
             dt.datetime(2009, 6, 23, 00, 00)]

    Ldate = [dt.datetime(2009, 6, 23, 00, 00)]
    LKdrift = [10,20,30]
    LKdrift = [1111]

    ITER = list(itertools.product(Ltraj,Lnbobs,Ldate,Ldesign,LKdrift))

    ## LA BOUCLE
    print(len(ITER),' experiences en batch')
    for i,itup in enumerate(ITER):
        print(i , '/' , len(ITER) ,' experiences en batch')

        traj  = itup[0]
        nbobs = itup[1]
        dat   = itup[2]
        design= itup[3]
        Kdrift = itup[4]


        dico = dict()

        dico['prm_design']   =  design

        dico['prm_zones_bound_input']  = [1000]
        dico['prm_sigma_zones_input']  = [10**-4,10**-5]

        dico['prm_with_temporal']   =  1
        dico['prm_tempor_start']    = dat
        dico['tempor_ctd_time']     = dat
        dico['prm_tempor_ctd_time'] = dat

        dico['prm_null_water_column_noise'] = 0
        dico['prm_coef_sigma_zones']        = 0

        dico['prm_with_ssf']           = 0
        dico['prm_with_pseudo_ssf']    = 0

        dico['prm_sigma_x'] = 0.
        dico['prm_sigma_y'] = 0.
        dico['prm_sigma_z'] = 0.

        if traj == 'droite':
            dico['prm_drt_nb_obs']    =  int(nbobs/3)
        else:
            dico['prm_drv_nb_obs']    =  nbobs
            dico["prm_drv_K_derive"]  = Kdrift

        dico["prm_add_offset"]    =  False
        dico["imporved_noising"]  =  1
        dico["prm_traject"]       =  traj
        dico['prm_drv_R']         =  50

        dico['procs']             = 7

        midfix = 'date' + str(dat.strftime('%Y%m%d')) + '_' + design + '_' + "Kdrift" + str(Kdrift)

        fab.fabrik_fct(exp_prefix,midfix,dico)




if 0:
    exp_prefix   = 'FINALa_VariTempor_zonal1056'
    exp_prefix   = 'FINALb_VariTempor_zonal1056_ss-bruit-surface'
    exp_prefix   = 'FINALb_VariTempor_zonal1045_ss-bruit-surface'

    Ltraj   = ['derive',"droite"]
    Ldesign = ['grand-carree',"petit-carree","triangle","grand-trigrangle"]

    Lnbobs  = [100,1000,10000]
    Ldate = [dt.datetime(2007, 4, 24, 00, 00),
             dt.datetime(2010, 2, 25, 00, 00),
             dt.datetime(2009, 6, 23, 00, 00)]
    ITER = list(itertools.product(Ltraj,Lnbobs,Ldate,Ldesign))

    ## LA BOUCLE
    print(len(ITER),' experiences en batch')
    for i,itup in enumerate(ITER):
        print(i , '/' , len(ITER) ,' experiences en batch')

        traj  = itup[0]
        nbobs = itup[1]
        dat   = itup[2]
        design= itup[3]

        dico = dict()

        dico['prm_design']   =  design

        dico['prm_zones_bound_input']  = [1000]
        dico['prm_sigma_zones_input']  = [10**-4,10**-5]

        dico['prm_with_temporal']   =  1
        dico['prm_tempor_start']    = dat
        dico['tempor_ctd_time']     = dat
        dico['prm_tempor_ctd_time'] = dat

        dico['prm_null_water_column_noise'] = 0
        dico['prm_coef_sigma_zones']        = 1

        dico['prm_with_ssf']           = 0
        dico['prm_with_pseudo_ssf']    = 0

        dico['prm_sigma_x'] = 0.0
        dico['prm_sigma_y'] = 0.0
        dico['prm_sigma_z'] = 0.0

        if traj == 'droite':
            dico['prm_drt_nb_obs']    =  nbobs/3
        else:
            dico['prm_drv_nb_obs']    =  nbobs

        dico["prm_add_offset"]    =  False
        dico["imporved_noising"]  =  1
        dico["prm_traject"]       =  traj
        dico['prm_drv_R']         =  50

        dico['procs']             = 7

        midfix = 'date' + str(dat.strftime('%Y%m%d')) + '_' + design

        fab.fabrik_fct(exp_prefix,midfix,dico)

if 0:
    #FOR GENERATION OF Gdic
    # Le Udic est en input

    exp_prefix         = 'GdicGENERATOR__triangleTST'
    exp_prefix_proto   = 'GdicGENERATOR_2_'

    exp_prefix_proto   = 'GdicGENERATOR_2_REBOOT0317_FOR_UNITARY_BENCHMARKING_'

    Lkouran = (True , False)
    Ldesign = ('triangle','petit-carree','grand-carree','grand-trigrangle')
    Lship   = ('atal','meteor')


    Ldesign = ('grand-carree',)


    for boolkouran,design,ship in itertools.product(Lkouran,Ldesign,Lship):

        dico = dict()

        exp_prefix = exp_prefix_proto + ship

        dico['prm_design'] = design

        dico['prm_Udic_path'] = toolbox_path + "/exemple/Udic/UdicFINAL_" + ship + '.pik'

        dico['prm_sigma_x'] = 0
        dico['prm_sigma_y'] = 0
        dico['prm_sigma_z'] = 0
        dico['prm_sigma_t'] = 0

        dico['prm_add_offset'] = 0
        dico['prm_null_water_column_noise'] = 1
        dico['prm_coef_sigma_zones']        = 0

        dico['prm_with_ssf']           = True
        dico['prm_bypass_ssf_rec']     = 1 #
        dico['prm_ssf_realistik_grad'] = 0
        dico['prm_with_pseudo_ssf']    = 0

        dico['prm_SSF_based_on_interpolation']                 = True
        dico['prm_SSF_based_on_interpolation_cental_SSP_ONLY'] = boolkouran

        dico['prm_eiko_h']          = 1000
        dico['prm_eiko_resotype']   = 'rkck'
        dico['prm_x_grad']          = 0
        dico['prm_y_grad']          = 0
        dico['prm_z_mixed']         = 3500 #500

        dico['prm_traject']     = 'droite'

        dico['prm_drt_xsize']   = 4100
        dico['prm_drt_ysize']   = 4100
        dico['prm_drt_xcenter'] = 10000
        dico['prm_drt_ycenter'] = 10000
        dico['prm_drt_zcenter'] = 0
        dico['prm_drt_nb_pass'] = 10
        dico['prm_drt_nb_obs']  = 10
        dico['prm_drt_angle']   = 90

        dico['procs']             = 6

        midfix = 'design_' + design + '_kouran' + str(int(boolkouran))

        fab.fabrik_fct(exp_prefix,midfix,dico)

if 0:
    exp_prefix   = 'FINALb_GradientEffect'
    Lnbobs  = [100,1000,10000]
    Lship   = ['ATAL','METEOR']
    Ltraj   = ['derive',"droite"]
    Ldesign = ['grand-carree',"petit-carree","triangle",'grand-trigrangle']
    Ldesign = ['grand-carree',"petit-carree"]
    Lnbobs  = [50000]

    ITER = list(itertools.product(Ltraj,Lnbobs,Ldesign,Lship))

    ## LA BOUCLE
    print(len(ITER),' experiences en batch')
    for i,itup in enumerate(ITER):
        print(i , '/' , len(ITER) ,' experiences en batch')

        traj  = itup[0]
        nbobs = itup[1]
        design= itup[2]
        ship  = itup[3]

        dico = dict()

        if design == 'grand-carree':
            if ship == 'ATAL':
                Gdicpath = toolbox_path + "/exemple/Gdic/Gdic_FINALa_20_ATAL_grand-carree.pik"
            if ship == 'METEOR':
                Gdicpath = toolbox_path + "/exemple/Gdic/Gdic_FINALa_20_METEOR_grand-carree.pik"

        elif design == "petit-carree":
            if ship == 'ATAL':
                Gdicpath = toolbox_path + "/exemple/Gdic/Gdic_FINALa_20_ATAL_petit-carree.pik"
            if ship == 'METEOR':
                Gdicpath = toolbox_path + "/exemple/Gdic/Gdic_FINALa_20_METEOR_petit-carree.pik"

        elif design == "triangle":
            if ship == 'ATAL':
                Gdicpath = toolbox_path + "/exemple/Gdic/Gdic_FINALa_20_ATAL_triangle.pik"
            if ship == 'METEOR':
                Gdicpath = toolbox_path + "/exemple/Gdic/Gdic_FINALa_20_METEOR_triangle.pik"

        elif design == "grand-trigrangle":
            if ship == 'ATAL':
                Gdicpath = toolbox_path + "/exemple/Gdic/Gdic_FINALa_20_ATAL_triangle.pik"
            if ship == 'METEOR':
                Gdicpath = toolbox_path + "/exemple/Gdic/Gdic_FINALa_20_METEOR_triangle.pik"

        Gdicpath = toolbox_path + "/exemple/Gdic/Gdic_FINALb_10_" + ship.lower() + '_' + design + '.pik'

        dico['prm_design']   =  design
        dico['prm_Gdic_path']   =  Gdicpath


        dico['prm_sigma_x'] = 0
        dico['prm_sigma_y'] = 0
        dico['prm_sigma_z'] = 0
        dico['prm_sigma_t'] = 0

        dico['prm_add_offset'] = 0
        dico['prm_null_water_column_noise'] = True
        dico['prm_coef_sigma_zones']        = 0

        dico['prm_with_ssf']           = False
        dico['prm_with_pseudo_ssf']    = True

        dico['prm_with_temporal']   =  False

        dico['prm_drt_angle']   = 90

        if traj == 'droite':
            dico['prm_drt_nb_obs']    =  nbobs/3
        else:
            dico['prm_drv_nb_obs']    =  nbobs

        dico["prm_add_offset"]    =  0
        dico["imporved_noising"]  =  1
        dico["prm_traject"]       =  traj
        dico['prm_drv_R']         =  50

        dico['procs']             = 7

        midfix =  ship + '_' + design

        fab.fabrik_fct(exp_prefix,midfix,dico)

if 0:
    exp_prefix   = 'FINALb_GradientEffect_rotat'
    Ltraj   = ["droite"]
    Ldesign = ['grand-carree',"petit-carree","triangle",'grand-trigrangle']
    Lnbobs  = [100,1000,10000]
    Lship  = ['ATAL','METEOR']

    ITER = list(itertools.product(Ltraj,Lnbobs,Ldesign,Lship))

    ## LA BOUCLE
    print(len(ITER),' experiences en batch')
    for i,itup in enumerate(ITER):
        print(i , '/' , len(ITER) ,' experiences en batch')

        traj  = itup[0]
        nbobs = itup[1]
        design= itup[2]
        ship  = itup[3]

        dico = dict()

        if design == 'grand-carree':
            if ship == 'ATAL':
                Gdicpath = toolbox_path + "/exemple/Gdic/Gdic_FINALa_20_ATAL_grand-carree.pik"
            if ship == 'METEOR':
                Gdicpath = toolbox_path + "/exemple/Gdic/Gdic_FINALa_20_METEOR_grand-carree.pik"

        elif design == "petit-carree":
            if ship == 'ATAL':
                Gdicpath = toolbox_path + "/exemple/Gdic/Gdic_FINALa_20_ATAL_petit-carree.pik"
            if ship == 'METEOR':
                Gdicpath = toolbox_path + "/exemple/Gdic/Gdic_FINALa_20_METEOR_petit-carree.pik"

        elif design == "triangle":
            if ship == 'ATAL':
                Gdicpath = toolbox_path + "/exemple/Gdic/Gdic_FINALa_20_ATAL_triangle.pik"
            if ship == 'METEOR':
                Gdicpath = toolbox_path + "/exemple/Gdic/Gdic_FINALa_20_METEOR_triangle.pik"

        elif design == "grand-trigrangle":
            if ship == 'ATAL':
                Gdicpath = toolbox_path + "/exemple/Gdic/Gdic_FINALa_20_ATAL_triangle.pik"
            if ship == 'METEOR':
                Gdicpath = toolbox_path + "/exemple/Gdic/Gdic_FINALa_20_METEOR_triangle.pik"

        Gdicpath = toolbox_path + 'acoustyx_toolbox_2/exemple/Gdic/Gdic_FINALb_10_' + ship.lower() + '_' + design + '.pik'

        dico['prm_design']   =  design
        dico['prm_Gdic_path']   =  Gdicpath

        dico['prm_sigma_x'] = 0
        dico['prm_sigma_y'] = 0
        dico['prm_sigma_z'] = 0
        dico['prm_sigma_t'] = 0

        dico['prm_add_offset'] = 0
        dico['prm_null_water_column_noise'] = True
        dico['prm_coef_sigma_zones']        = 0

        dico['prm_with_ssf']           = False
        dico['prm_with_pseudo_ssf']    = True

        dico['prm_with_temporal']   =  False

        dico['prm_drt_angle']   = 90

        if traj == 'droite':
            dico['prm_drt_nb_obs']    =  nbobs/3
        else:
            dico['prm_drv_nb_obs']    =  nbobs

        dico["prm_add_offset"]    =  0
        dico["imporved_noising"]  =  1
        dico["prm_traject"]       =  traj
        dico['prm_drv_R']         =  50

        dico['procs']             = 7

        midfix =  ship + '_' + design

        fab.fabrik_fct(exp_prefix,midfix,dico)

if 0:
    exp_prefix  = 'FINALb_BUOYstyle_mk6d_zonal1056'

    nbdays                = 365
    L_prm_tempor_start    = [dt.datetime(2012, 1, 1, 12, 00) + dt.timedelta(days=i) for i in range(nbdays)]

    intern_random_seed_x = 789
    intern_random_seed_y = 456
    sigma_posi_x  = 50
    sigma_posi_y  = 50

    Rx = np.random.RandomState(intern_random_seed_x)
    Ry = np.random.RandomState(intern_random_seed_y)

    L_prm_x0 = 10000 + Rx.randn(365) * sigma_posi_x
    L_prm_y0 = 10000 + Ry.randn(365) * sigma_posi_y

    prm_tempor_ctd_time = L_prm_tempor_start[0]

    Ldesign = ['grand-carree',"petit-carree","triangle",'grand-trigrangle']

    ITER = list(zip(L_prm_tempor_start,L_prm_x0,L_prm_y0))

    error_counter = 0

    # LA BOUCLE
    print(len(ITER),' experiences en batch')
    for design in Ldesign:
        for i,itup in enumerate(ITER):
            print(i , '/' , len(ITER) ,' experiences en batch')

            dico = dict()
            dico['prm_with_temporal']   = True
            dico['prm_tempor_start']    = itup[0]
            dico['prm_x0']              = itup[1]
            dico['prm_y0']              = itup[2]
            dico['prm_tempor_ctd_time'] = itup[0] #prm_tempor_ctd_time
            dico['prm_tempor_len']      = 3600*4
            dico['prm_drv_nb_obs']      = 100
            dico['prm_drv_R']           = 50
            dico['prm_drv_K_derive']    = (dico['prm_drv_nb_obs'] + i) * 1
            dico['prm_K_t_zone']        = (dico['prm_drv_nb_obs'] + i) * 42

            dico['prm_design']   =  design


            dico["prm_add_offset"]    =  False
            dico["imporved_noising"]  =  1

            dico["prm_traject"]       =  'derive'

            dico['procs']               = 2

            dico['prm_with_ssf']           = 0
            dico['prm_with_pseudo_ssf']    = 0
            dico['prm_null_water_column_noise'] = 1
            dico['prm_coef_sigma_zones']        = 0


            if 0: #Si l'on veut appliquer un gradient
                dico['prm_with_pseudo_ssf']    = 1
                ship = 'ATAL'
                if design == 'grand-carree':
                    if ship == 'ATAL':
                        Gdicpath = toolbox_path + "/exemple/Gdic/Gdic_FINALa_20_ATAL_grand-carree.pik"
                    if ship == 'METEOR':
                        Gdicpath = toolbox_path + "/exemple/Gdic/Gdic_FINALa_20_METEOR_grand-carree.pik"

                elif design == "petit-carree":
                    if ship == 'ATAL':
                        Gdicpath = toolbox_path + "/exemple/Gdic/Gdic_FINALa_20_ATAL_petit-carree.pik"
                    if ship == 'METEOR':
                        Gdicpath = toolbox_path + "/exemple/Gdic/Gdic_FINALa_20_METEOR_petit-carree.pik"

                elif design == "triangle":
                    if ship == 'ATAL':
                        Gdicpath = toolbox_path + "/exemple/Gdic/Gdic_FINALa_20_ATAL_triangle.pik"
                    if ship == 'METEOR':
                        Gdicpath = toolbox_path + "/exemple/Gdic/Gdic_FINALa_20_METEOR_triangle.pik"

                elif design == "grand-trigrangle":
                    if ship == 'ATAL':
                        Gdicpath = toolbox_path + "/exemple/Gdic/Gdic_FINALa_20_ATAL_triangle.pik"
                    if ship == 'METEOR':
                        Gdicpath = toolbox_path + "/exemple/Gdic/Gdic_FINALa_20_METEOR_triangle.pik"

                Gdicpath = toolbox_path + 'acoustyx_toolbox_2/exemple/Gdic/Gdic_FINALb_10_' + ship.lower() + '_' + design + '.pik'

                dico['prm_Gdic_path']   =  Gdicpath
            elif 1:

                dico['prm_null_water_column_noise'] = 0
                dico['prm_coef_sigma_zones']        = 1
                dico['prm_zones_bound_input']  = [1000]
                dico['prm_sigma_zones_input']  = [10**-5,10**-6]


            midfix = gf.join_improved('_','BUOYstyle',itup[0].date().strftime('%Y%m%d'),
                                      design , *[np.round(e) for e in itup[1:3]])

            try:
                fab.fabrik_fct(exp_prefix,midfix,dico)
            except Exception as e:
                print("ERR : for " ,  itup)
                error_counter += 1
                print("error_counter" ,  error_counter)
                continue

if 0:
    exp_prefix  = 'FINALb_BUOYstyle_mk6ef_no-zonal-no-traj-noise'
    exp_prefix  = 'FINALb_BUOYstyle_mk6g_no-zonal-traj-noise-new-twister'
    exp_prefix  = 'FINALb_BUOYstyle_mk8i2_zonal1056'
    exp_prefix  = 'FINALb_BUOYstyle_mk8h2_Gradient'

    exp_prefix  = 'FINALb_BUOYstyle_mk9bordertraj_Gradient'
    exp_prefix  = 'FINALb_BUOYstyle_mk9bordertraj_zonal1056'
    exp_prefix  = 'FINALb_BUOYstyle_mk8e2_no-zonal'


    nbdays                = 365
    L_prm_tempor_start    = [dt.datetime(2012, 0o1, 1, 12, 00) + dt.timedelta(days=i) for i in range(nbdays)]

    intern_random_seed_x = 789
    intern_random_seed_y = 456
    sigma_posi_x  = 500
    sigma_posi_y  = 500

    Rx = np.random.RandomState(intern_random_seed_x)
    Ry = np.random.RandomState(intern_random_seed_y)

    if 0: # bouée centrée
        L_prm_x0 = 10000 + Rx.randn(365) * sigma_posi_x
        L_prm_y0 = 10000 + Ry.randn(365) * sigma_posi_y
    else: # bouée sur les bords
        L_prm_x0 ,L_prm_y0 = geok.points_circle_border(1000,1000,250,1,
                                        np.pi/2,dir_range=1,
                                        seed=intern_random_seed_x)
        L_prm_x0 = L_prm_x0 + 10000
        L_prm_y0 = L_prm_y0 + 10000


    prm_tempor_ctd_time = L_prm_tempor_start[0]

    Ldesign = ['grand-carree']

    ITER = list(zip(L_prm_tempor_start,L_prm_x0,L_prm_y0))

    error_counter = 0

    # LA BOUCLE
    print(len(ITER),' experiences en batch')
    for design in Ldesign:
        for i,itup in enumerate(ITER):
            print(i , '/' , len(ITER) ,' experiences en batch')

            dico = dict()
            dico['prm_with_temporal']   = True
            dico['prm_tempor_start']    = itup[0]
            dico['prm_drv_x0']          = itup[1]
            dico['prm_drv_y0']          = itup[2]
            dico['prm_drv_xcenter']     = itup[1]
            dico['prm_drv_ycenter']     = itup[2] # center of the circle
            dico['prm_tempor_ctd_time'] = itup[0] #prm_tempor_ctd_time
            dico['prm_tempor_len']      = 3600*4
            dico['prm_drv_nb_obs']      = 100
            dico['prm_drv_R']           = 500
            dico['prm_drv_K_derive']    = (dico['prm_drv_nb_obs'] + i) * 1
            dico['prm_K_t_zone']        = (dico['prm_drv_nb_obs'] + i) * 42

            dico['prm_K_xyz']           = i * 17

            dico['prm_design']   =  design


            dico["prm_add_offset"]    =  False
            dico["imporved_noising"]  =  1

            dico["prm_traject"]       =  'derive'

            dico['procs']               = 6

            dico['prm_with_ssf']           = 0
            dico['prm_with_pseudo_ssf']    = 0
            dico['prm_null_water_column_noise'] = 1
            dico['prm_coef_sigma_zones']        = 0


            if 0: #Si l'on veut appliquer un gradient
                dico['prm_with_pseudo_ssf']    = 1
                ship = 'ATAL'
                if design == 'grand-carree':
                    if ship == 'ATAL':
                        Gdicpath = toolbox_path + "/exemple/Gdic/Gdic_FINALa_20_ATAL_grand-carree.pik"
                    if ship == 'METEOR':
                        Gdicpath = toolbox_path + "/exemple/Gdic/Gdic_FINALa_20_METEOR_grand-carree.pik"

                elif design == "petit-carree":
                    if ship == 'ATAL':
                        Gdicpath = toolbox_path + "/exemple/Gdic/Gdic_FINALa_20_ATAL_petit-carree.pik"
                    if ship == 'METEOR':
                        Gdicpath = toolbox_path + "/exemple/Gdic/Gdic_FINALa_20_METEOR_petit-carree.pik"

                elif design == "triangle":
                    if ship == 'ATAL':
                        Gdicpath = toolbox_path + "/exemple/Gdic/Gdic_FINALa_20_ATAL_triangle.pik"
                    if ship == 'METEOR':
                        Gdicpath = toolbox_path + "/exemple/Gdic/Gdic_FINALa_20_METEOR_triangle.pik"

                elif design == "grand-trigrangle":
                    if ship == 'ATAL':
                        Gdicpath = toolbox_path + "/exemple/Gdic/Gdic_FINALa_20_ATAL_triangle.pik"
                    if ship == 'METEOR':
                        Gdicpath = toolbox_path + "/exemple/Gdic/Gdic_FINALa_20_METEOR_triangle.pik"

                Gdicpath = toolbox_path + 'acoustyx_toolbox_2/exemple/Gdic/Gdic_FINALb_10_' + ship.lower() + '_' + design + '.pik'

                dico['prm_Gdic_path']   =  Gdicpath
            elif 0:
                dico['prm_null_water_column_noise'] = 0
                dico['prm_coef_sigma_zones']        = 1
                dico['prm_zones_bound_input']  = [1000]
                dico['prm_sigma_zones_input']  = [10**-5,10**-6]

            elif 1:
                dico['prm_null_water_column_noise'] = 1
                dico['prm_coef_sigma_zones']        = 0
                dico['prm_zones_bound_input']  = [1000]
                dico['prm_sigma_zones_input']  = [0,0]

            dico['prm_sigma_x'] = 0.03
            dico['prm_sigma_y'] = 0.03
            dico['prm_sigma_z'] = 0.05

            midfix = gf.join_improved('_','BUOYstyle',itup[0].date().strftime('%Y%m%d'),
                                      design , *[np.round(e) for e in itup[1:3]])

            try:
                fab.fabrik_fct(exp_prefix,midfix,dico)
            except Exception as e:
                print("ERR : for " ,  itup)
                error_counter += 1
                print("error_counter" ,  error_counter)
                continue

if 0:
    exp_prefix   = 'FINALc_VariTemporTest'

    Ldesign = ['grand-carree']
    Ltraj   = ['station']
    Lnbobs  = [10000]
    Ldate = [dt.datetime(2009, 10, 0o2, 00, 00)]

    ITER = list(itertools.product(Ltraj,Lnbobs,Ldate,Ldesign))

    ## LA BOUCLE
    print(len(ITER),' experiences en batch')
    for i,itup in enumerate(ITER):
        print(i , '/' , len(ITER) ,' experiences en batch')

        traj  = itup[0]
        nbobs = itup[1]
        dat   = itup[2]
        design= itup[3]

        dico = dict()

        dico['prm_design']   =  design

        dico['prm_zones_bound_input']  = [1000]
        dico['prm_sigma_zones_input']  = [10**-4,10**-5]

        dico['prm_with_temporal']   = 1
        dico['prm_tempor_start']    = dat
        dico['tempor_ctd_time']     = dat
        dico['prm_tempor_ctd_time'] = dat

        dico['prm_null_water_column_noise'] = 1
        dico['prm_coef_sigma_zones']        = 0

        dico['prm_with_ssf']           = 0
        dico['prm_with_pseudo_ssf']    = 0

        dico['prm_sigma_x'] = 0.0
        dico['prm_sigma_y'] = 0.0
        dico['prm_sigma_z'] = 0.0

        if traj == 'droite':
            dico['prm_drt_nb_obs']    =  nbobs/3
        if traj == 'station':
            dico['prm_sta_nb_obs']    =  nbobs
        else:
            dico['prm_drv_nb_obs']    =  nbobs

        dico["prm_add_offset"]    =  False
        dico["imporved_noising"]  =  1
        dico["prm_traject"]       =  traj
        dico['prm_drv_R']         =  10
        dico['prm_vit']           =  1

        dico['procs']             = 7

        midfix = 'date' + str(dat.strftime('%Y%m%d')) + '_' + design

        fab.fabrik_fct(exp_prefix,midfix,dico)


if 0:
    exp_prefix   = 'FINALc_VariTemporTest_REBOOT1706_1m-diam_cental'

    Ldesign = ['grand-carree']
    Ltraj   = ['station']
    Lnbobs  = [1000]
    Ldate   = [dt.datetime(2009, 10, 0o2, 00, 00)]
    Lcentonoff = [True,False]

    ITER = list(itertools.product(Ltraj,Lnbobs,Ldate,Ldesign,Lcentonoff))

    ## LA BOUCLE
    print(len(ITER),' experiences en batch')
    for i,itup in enumerate(ITER):
        print(i , '/' , len(ITER) ,' experiences en batch')

        traj  = itup[0]
        nbobs = itup[1]
        dat   = itup[2]
        design= itup[3]
        centalonoff = itup[4]

        dico = dict()

        dico['prm_design']   =  design

        dico['prm_zones_bound_input']  = [1000]
        dico['prm_sigma_zones_input']  = [10**-4,10**-5]

        dico['prm_with_temporal']   = 0
        dico['prm_tempor_start']    = dat
        dico['tempor_ctd_time']     = dat
        dico['prm_tempor_ctd_time'] = dat

        dico['prm_SSF_based_on_interpolation_cental_SSP_ONLY'] = centalonoff

        dico['prm_null_water_column_noise'] = 1
        dico['prm_coef_sigma_zones']        = 0

        dico['prm_with_ssf']           = 1
        dico['prm_with_pseudo_ssf']    = 0

        dico['prm_sigma_x'] = 0.0
        dico['prm_sigma_y'] = 0.0
        dico['prm_sigma_z'] = 0.0

        if traj == 'droite':
            dico['prm_drt_nb_obs']    =  nbobs/3
        if traj == 'station':
            dico['prm_sta_nb_obs']    =  nbobs
        else:
            dico['prm_drv_nb_obs']    =  nbobs

        dico["prm_add_offset"]    =  False
        dico["imporved_noising"]  =  1
        dico["prm_traject"]       =  traj
        dico['prm_drv_R']         =  1
        dico['prm_vit']           =  1

        dico['procs']             = 4


        midfix = 'date' + str(dat.strftime('%Y%m%d')) + '_' + design

        fab.fabrik_fct(exp_prefix,midfix,dico)

if 0:
    exp_prefix   = '1709A_10m-diam-central-withorwithoutLATERALgrad'

    Ldesign = ['grand-carree']
    Ltraj   = ['derive']
    Lnbobs  = [1000]
    Ldate   = [dt.datetime(2009, 10, 0o2, 00, 00)]


    Lcentonoff = [True,False]

    ITER = list(itertools.product(Ltraj,Lnbobs,Ldate,Ldesign,Lcentonoff))

    ## LA BOUCLE
    print(len(ITER),' experiences en batch')
    for i,itup in enumerate(ITER):
        print(i , '/' , len(ITER) ,' experiences en batch')

        traj  = itup[0]
        nbobs = itup[1]
        dat   = itup[2]
        design= itup[3]
        centalonoff = itup[4]

        dico = dict()


        dico['prm_design']   =  design

        dico['prm_zones_bound_input']  = [1000]
        dico['prm_sigma_zones_input']  = [10**-4,10**-5]

        dico['prm_with_temporal']   = 0
        dico['prm_tempor_start']    = dat
        dico['tempor_ctd_time']     = dat
        dico['prm_tempor_ctd_time'] = dat

        dico['prm_SSF_based_on_interpolation_cental_SSP_ONLY'] = centalonoff

        dico['prm_null_water_column_noise'] = 1
        dico['prm_coef_sigma_zones']        = 0

        dico['prm_with_ssf']           = 1
        dico['prm_with_pseudo_ssf']    = 0

        dico['prm_sigma_x'] = 0.0
        dico['prm_sigma_y'] = 0.0
        dico['prm_sigma_z'] = 0.0

        if traj == 'droite':
            dico['prm_drt_nb_obs']    =  nbobs/3
        if traj == 'station':
            dico['prm_sta_nb_obs']    =  nbobs
        else:
            dico['prm_drv_nb_obs']    =  nbobs

        dico["prm_add_offset"]    =  False
        dico["imporved_noising"]  =  1
        dico["prm_traject"]       =  traj
        dico['prm_drv_R']         =  10
        dico['prm_vit']           =  1

        dico['procs']             = 7


        midfix = 'date' + str(dat.strftime('%Y%m%d')) + '_' + design + '_' + 'LAT_GRAD'  + str(not centalonoff)

        fab.fabrik_fct(exp_prefix,midfix,dico)

if 0:
    # will produce only one exp, with the SSP mode, as a basic and stable exemple
    exp_prefix   = '1709B_10m-diam-central-SSP-as-BASIS'

    Ldesign = ['grand-carree']
    Ltraj   = ['derive']
    Lnbobs  = [1000]
    Ldate   = [dt.datetime(2009, 10, 0o2, 00, 00)]


    Lcentonoff = [False]

    ITER = list(itertools.product(Ltraj,Lnbobs,Ldate,Ldesign,Lcentonoff))

    ## LA BOUCLE
    print(len(ITER),' experiences en batch')
    for i,itup in enumerate(ITER):
        print(i , '/' , len(ITER) ,' experiences en batch')

        traj  = itup[0]
        nbobs = itup[1]
        dat   = itup[2]
        design= itup[3]
        centalonoff = itup[4]

        dico = dict()


        dico['prm_design']   =  design

        dico['prm_zones_bound_input']  = [1000]
        dico['prm_sigma_zones_input']  = [10**-4,10**-5]

        dico['prm_with_temporal']   = 0
        dico['prm_tempor_start']    = dat
        dico['tempor_ctd_time']     = dat
        dico['prm_tempor_ctd_time'] = dat

        dico['prm_SSF_based_on_interpolation_cental_SSP_ONLY'] = centalonoff

        dico['prm_null_water_column_noise'] = 1
        dico['prm_coef_sigma_zones']        = 0

        dico['prm_with_ssf']           = 0
        dico['prm_with_pseudo_ssf']    = 0

        dico['prm_sigma_x'] = 0.0
        dico['prm_sigma_y'] = 0.0
        dico['prm_sigma_z'] = 0.0

        if traj == 'droite':
            dico['prm_drt_nb_obs']    =  nbobs/3
        if traj == 'station':
            dico['prm_sta_nb_obs']    =  nbobs
        else:
            dico['prm_drv_nb_obs']    =  nbobs

        dico["prm_add_offset"]    =  False
        dico["imporved_noising"]  =  1
        dico["prm_traject"]       =  traj
        dico['prm_drv_R']         =  10
        dico['prm_vit']           =  1

        dico['procs']             = 7


        midfix = 'date' + str(dat.strftime('%Y%m%d')) + '_' + design + '_' + 'LAT_GRAD'  + str(not centalonoff)

        fab.fabrik_fct(exp_prefix,midfix,dico)

if 0:
    exp_prefix   = '1709C_10m-diam-central-withorwithoutLATERALgrad-butbackwardBYPASSED'

    Ldesign = ['grand-carree']
    Ltraj   = ['derive']
    Lnbobs  = [1000]
    Ldate   = [dt.datetime(2009, 10, 0o2, 00, 00)]


    Lcentonoff = [True,False]

    ITER = list(itertools.product(Ltraj,Lnbobs,Ldate,Ldesign,Lcentonoff))

    ## LA BOUCLE
    print(len(ITER),' experiences en batch')
    for i,itup in enumerate(ITER):
        print(i , '/' , len(ITER) ,' experiences en batch')

        traj  = itup[0]
        nbobs = itup[1]
        dat   = itup[2]
        design= itup[3]
        centalonoff = itup[4]

        dico = dict()


        dico['prm_design']   =  design

        dico['prm_zones_bound_input']  = [1000]
        dico['prm_sigma_zones_input']  = [10**-4,10**-5]

        dico['prm_with_temporal']   = 0
        dico['prm_tempor_start']    = dat
        dico['tempor_ctd_time']     = dat
        dico['prm_tempor_ctd_time'] = dat

        dico['prm_SSF_based_on_interpolation_cental_SSP_ONLY'] = centalonoff

        dico['prm_null_water_column_noise'] = 1
        dico['prm_coef_sigma_zones']        = 0

        dico['prm_with_ssf']           = 1
        dico['prm_with_pseudo_ssf']    = 0

        dico['prm_bypass_ssf_rec'] = 1

        dico['prm_sigma_x'] = 0.0
        dico['prm_sigma_y'] = 0.0
        dico['prm_sigma_z'] = 0.0

        if traj == 'droite':
            dico['prm_drt_nb_obs']    =  nbobs/3
        if traj == 'station':
            dico['prm_sta_nb_obs']    =  nbobs
        else:
            dico['prm_drv_nb_obs']    =  nbobs

        dico["prm_add_offset"]    =  False
        dico["imporved_noising"]  =  1
        dico["prm_traject"]       =  traj
        dico['prm_drv_R']         =  10
        dico['prm_vit']           =  1

        dico['procs']             = 7


        midfix = 'date' + str(dat.strftime('%Y%m%d')) + '_' + design + '_' + 'LAT_GRAD'  + str(not centalonoff)

        fab.fabrik_fct(exp_prefix,midfix,dico)
if 0:
    exp_prefix   = '1709E_10m-diam-central-withorwithoutLATERALgrad-pseudoSSFavecDic'

    Ldesign = ['grand-carree']
    Ltraj   = ['derive']
    Lnbobs  = [1000]
    Ldate   = [dt.datetime(2009, 10, 0o2, 00, 00)]


    Lcentonoff = [True,False]

    ITER = list(itertools.product(Ltraj,Lnbobs,Ldate,Ldesign,Lcentonoff))

    ## LA BOUCLE
    print(len(ITER),' experiences en batch')
    for i,itup in enumerate(ITER):
        print(i , '/' , len(ITER) ,' experiences en batch')

        traj  = itup[0]
        nbobs = itup[1]
        dat   = itup[2]
        design= itup[3]
        centalonoff = itup[4]

        dico = dict()


        dico['prm_design']   =  design

        dico['prm_zones_bound_input']  = [1000]
        dico['prm_sigma_zones_input']  = [10**-4,10**-5]

        dico['prm_with_temporal']   = 0
        dico['prm_tempor_start']    = dat
        dico['tempor_ctd_time']     = dat
        dico['prm_tempor_ctd_time'] = dat

        #dico['prm_SSF_based_on_interpolation_cental_SSP_ONLY'] = centalonoff

        dico['prm_null_water_column_noise'] = 1
        dico['prm_coef_sigma_zones']        = 0

        dico['prm_with_ssf']           = 0
        dico['prm_with_pseudo_ssf']    = centalonoff

        dico['prm_bypass_ssf_rec'] = 1

        dico['prm_sigma_x'] = 0.0
        dico['prm_sigma_y'] = 0.0
        dico['prm_sigma_z'] = 0.0

        if traj == 'droite':
            dico['prm_drt_nb_obs']    =  nbobs/3
        if traj == 'station':
            dico['prm_sta_nb_obs']    =  nbobs
        else:
            dico['prm_drv_nb_obs']    =  nbobs

        dico["prm_add_offset"]    =  False
        dico["imporved_noising"]  =  1
        dico["prm_traject"]       =  traj
        dico['prm_drv_R']         =  10
        dico['prm_vit']           =  1

        dico['procs']             = 7


        midfix = 'date' + str(dat.strftime('%Y%m%d')) + '_' + design + '_' + 'pseudoLAT_GRAD'  + str(not centalonoff)

        fab.fabrik_fct(exp_prefix,midfix,dico)
