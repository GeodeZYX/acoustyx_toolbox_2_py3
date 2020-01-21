#  -*- coding: utf-8 -*-
"""
Éditeur de Spyder

Ceci est un script temporaire.
"""

from megalib import *

exp_path = '/home/pierre/Documents/CODES/acoustyx_toolbox_2/working/batch_3x50_x2000_y2000_nois1-1e-06_'
exp_path = '/home/pierre/aaa_FOURBI/formatD/'
exp_path = '/home/pierre/aaa_FOURBI/formatE/'
exp_path = '/home/pierre/Documents/CODES/acoustyx_toolbox_2/working/batch_3x1500_x0_y0_nois1-1e-06_'
exp_path = '/home/pierre/Documents/CODES/acoustyx_toolbox_2/working/batch_3x1500_x2000_y2000_nois1-1e-06_'
exp_path = '/home/pierre/Documents/CODES/acoustyx_toolbox_2/working/batch_3x500_x2000_y2000_nois1-1e-06_'
exp_path = "/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/working/batc2_3x1500_x0_y0_nois1-1e-06_"
exp_path = "/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/working/batc2c_deriv_decalZ_1000_R50_noisTrue-1e-06__/Results/compar_av_et_ss_dzmaster"
exp_path = "/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/working/batc2c_deriv_decalZ_1000_R50_noisTrue-1e-06__" 

exp_path_lis = glob.glob('/home/pierre/Documents/CODES/acoustyx_toolbox_2/working/batch*')
exp_path_lis = glob.glob("/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/working/batc2b_deriv_*")
exp_path_lis = glob.glob('/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/working/BEFORE_150915/batc2_3x1500_x2000_y2000_nois1-1e-06_')
exp_path_lis = glob.glob('/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/working/BEFORE_150915/batc2_3x1500_x0_y0_nois1-1e-06_*')
exp_path_lis = glob.glob('/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/working/batc2b*')
exp_path_lis = glob.glob('/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/working/BEFORE_150915/batc2_*')
exp_path_lis = glob.glob('/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/working/batc2c_deriv_decalZ_10000_R50_noisTrue-1e-06__')
exp_path_lis = glob.glob('/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/working/batc2*')
exp_path_lis = glob.glob("/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/working/BEFORE_150915/BATC2_1508/batc2*")
exp_path_lis = glob.glob('/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/working/RESULTS_END_1510/*')
exp_path_lis = glob.glob('/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/working/batc2c_deriv_decalZ_10000_R50_noisTrue-1e-06__')
exp_path_lis = glob.glob("/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/working/batc3*")
exp_path_lis = glob.glob("/home/psakicki/THESE/DATA/1506_GEODESEA/WORKING/OPERA3/*")
exp_path_lis = glob.glob("/home/psakicki/THESE/CODES/CodePython/acoustyx_toolbox_2/working/ADEL9993*")

#exp_path_lis = [exp_path]

if gf.spyder_run_check():
    plt.ioff()

# LES PLOTS
if 1:
    for exp_path in exp_path_lis:
        print(exp_path)
        if 0:
            out_sum_fil = acls.exp_summary(exp_path)
        if 1:
            try:
                #fig , ax = acls.plot_cntr_from_expdics(exp_path)
                fig , ax = acls.plot_cntr_from_expdics(exp_path,centroid_ref_point=1,variable_lis=['with_jackknife','with_V_4_P_reinject','jackknife inverted'])
                acls.plot_dist_from_expdics(exp_path)
            except:
                continue
            
#LE TABLEAU
if 1:
    for exp_input in exp_path_lis:
        F = genefun.Tee_frontend( exp_input , 'tab')   
        for exp_input in exp_path_lis:
            bigtitle = ' '.join(('**************', os.path.basename(exp_input),'**************'))
            print(bigtitle)
            
            outdir=''
            path_input=True
            verbose=False
                
            if path_input and outdir == '':
                outdir = exp_input
                
            if path_input:
                explis = glob.glob(exp_input + '/*exp')
            else:
                explis = exp_input     
            
            keys_lis_lis = []
            keys_lis = []
            
            if explis == []:
                print('ERR : exp_summary : explis is empty ... abort')
            
            final_iter_dic_list = []
            exp_name_lis        = []
            
            for exp in explis:
                dic = genefun.pickle_loader(exp)
                if -1 in list(dic.keys()):
                    final_iter_dic = dic[-1]
                else:
                    final_iter_dic = dic[max(dic.keys())]
                for k in list(final_iter_dic.keys()):
                    if not k in keys_lis:
                        keys_lis.append(k)
                keys_lis_lis.append(list(final_iter_dic.keys()))
                final_iter_dic_list.append(final_iter_dic)
                exp_name_lis.append(os.path.basename(exp))
                
            bary2d_diff_stk     = []
            params_stk          = []
            for fid in final_iter_dic_list:
                print('INFO : ' , list(fid.keys()) , len(list(fid.keys())))
                bary2d_diff_stk.append(fid['ecart 2D bary brut/vrai en distance'])
                params_stk.append(fid["params var. de l'exp."])
            
            bary2d_diff_stk = np.array(bary2d_diff_stk)
            good_bary       = bary2d_diff_stk < .1
            
            all_varparams_keys = []
            for trexp in params_stk:
                all_varparams_keys = all_varparams_keys + list(trexp.keys())
            all_varparams_keys     = list(set(all_varparams_keys))   
            
            all_varparams_keys.append('ecart 2D bary brut/vrai en distance')
            all_varparams_keys.append('ecart 2D bary MC/vrai en distance')
            
            # (np.max([len(e) for e in exp_name_lis])+1)
            for gob in range(2):
                title = ' ' * 25  + ' '.join(all_varparams_keys)
                print(title)
                gob = np.logical_not(gob)
                for i ,(exp , fid) in enumerate(zip(exp_name_lis , final_iter_dic_list)):
                    if gob == good_bary[i]:
                        line = exp[-24:] + ' '
                        for k in all_varparams_keys:
                            protofmt = '{:^' + str(len(k)) + '}'
                            if k in list(fid.keys()):
                                linestk = protofmt.format(fid[k])  
                            elif k in list(fid["params var. de l'exp."].keys()):
                                linestk = protofmt.format(fid["params var. de l'exp."][k])                      
                            else:
                                linestk = protofmt.format('N/A')
                            line = line +  linestk + ' ' 
                        print(line)
        F.stop()
                        