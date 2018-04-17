import numpy as np
import glob
import os
def load_summary(filename):
    dtype=[('minr', 'f8'),
           ('maxr', 'f8'), 
           ('ca_ratio', 'f8'),
           ('ba_ratio', 'f8'),
           ('a', 'f8'),
           ('center', 'f8'),
           ('width', 'f8'),
           ('mu', 'f8')]
    summary = np.loadtxt(filename, dtype=dtype)    
    return summary

def load_experiment(input_path="../data/mstar_selected_summary/vmax_sorted/", n_sat=11, full_data=False):
    files = glob.glob(input_path+"M31_group_*_nsat_{}.dat".format(n_sat))
    group_id = []
    for f in files:
        i = int(f.split("_")[-3])
        if i not in group_id:
            group_id.append(i)
    print(group_id, len(group_id))

    n_groups = len(group_id)
    
    fields = ['width','mu', 'a', 'ba_ratio', 'ca_ratio']
    M31_all = {}
    MW_all = {}
    if full_data:
        for field in fields:
            M31_all[field] = np.empty((0))
            MW_all[field] = np.empty((0))
            M31_all[field+'_random'] = np.empty((0))
            MW_all[field+'_random'] = np.empty((0))
    else:
        for field in fields:
            M31_all[field] = np.ones(n_groups)
            MW_all[field] = np.ones(n_groups)
            M31_all[field+'_sigma'] = np.ones(n_groups)
            MW_all[field+'_sigma'] = np.ones(n_groups)
        
            M31_all[field+'_random'] = np.ones(n_groups)
            MW_all[field+'_random'] = np.ones(n_groups)
            M31_all[field+'_random_sigma'] = np.ones(n_groups)
            MW_all[field+'_random_sigma'] = np.ones(n_groups)

    MW_summary = {}
    M31_summary = {}
    for g in range(n_groups):
        filename_MW = os.path.join(input_path,"MW_group_{}_nsat_{}.dat".format(group_id[g],n_sat))
        filename_M31 = os.path.join(input_path,"M31_group_{}_nsat_{}.dat".format(group_id[g],n_sat))

        MW_summary[g] = load_summary(filename_MW)
        M31_summary[g] = load_summary(filename_M31)
    
    
    for field in fields:
        a = np.empty((0))
        b = np.empty((0))
        a_random = np.empty((0))
        b_random = np.empty((0))
        
        for g in range(n_groups):
            data = M31_summary[g]
            a = data[field][0]
            a_random = data[field][1:101]
        
            data = MW_summary[g]
            b = data[field][0]
            b_random = data[field][1:101]
                #print('a_random {} iter: {} {}'.format(field, i, a_random))
           
            if full_data:
                M31_all[field] = np.append(M31_all[field], a)
                MW_all[field] = np.append(MW_all[field], b)
                M31_all[field+'_random'] = np.append(M31_all[field+'_random'], a_random)
                MW_all[field+'_random'] = np.append(MW_all[field+'_random'], b_random)
        
            else:
                M31_all[field][g] = np.average(a)
                MW_all[field][g] = np.average(b)
                M31_all[field+'_sigma'][g] = np.std(a)
                MW_all[field+'_sigma'][g] = np.std(b)
                M31_all[field+'_random'][g] = np.average(a_random)
                MW_all[field+'_random'][g] = np.average(b_random)
                M31_all[field+'_random_sigma'][g] = np.std(a_random)
                MW_all[field+'_random_sigma'][g] = np.std(b_random)
                
    return M31_all, MW_all


fields = ['width','ca_ratio', 'ba_ratio']
names = {'width':'Plane width (kpc)', 'ca_ratio':'$c/a$ ratio', 'ba_ratio':'$b/a$ ratio'}

for n_sat in range(11,16):
    print("OBSERVATIONS - NSAT = {}".format(n_sat))
    in_path = "../data/obs_summary/"
    M31_obs_stats, MW_obs_stats = load_experiment(input_path=in_path, n_sat=n_sat, full_data=False)

    print("M31(phys) - MW(phys) | M31(rand) | MW (rand) | M31(norm) | MW (norm)|")
    for field in fields:
        print("{} & ${:.2f}$ & ${:.2f}$ & ${:.2f}\pm{:.2f}$ & ${:.2f}\pm{:.2f}$ & ${:.2f}$ & ${:.2f}$\\\\\\hline".format(
            names[field], 
            M31_obs_stats[field][0], MW_obs_stats[field][0],
              M31_obs_stats[field+'_random'][0], M31_obs_stats[field+'_random_sigma'][0],
                MW_obs_stats[field+'_random'][0],MW_obs_stats[field+'_random_sigma'][0],
              (M31_obs_stats[field][0]-M31_obs_stats[field+'_random'][0])/M31_obs_stats[field+'_random_sigma'][0],
            (MW_obs_stats[field][0]-M31_obs_stats[field+'_random'][0])/MW_obs_stats[field+'_random_sigma'][0]))
    print()
    
for n_sat in range(11,16):
    print("Illustris - NSAT = {}".format(n_sat))
    
    in_path = "../data/illustris1_mstar_selected_summary/"
    M31_illu_stats, MW_illu_stats = load_experiment(input_path=in_path, n_sat=n_sat)
    in_path = "../data/illustris1dark_mstar_selected_summary/"
    M31_illudark_stats, MW_illudark_stats = load_experiment(input_path=in_path, n_sat=n_sat)
    in_path = "../data/elvis_mstar_selected_summary/"
    M31_elvis_stats, MW_elvis_stats = load_experiment(input_path=in_path, n_sat=n_sat)
    
    print("M31(phys) - MW(phys) | M31(rand) | MW (rand) | M31(norm) | MW (norm)|")
    for field in fields:
        print("{} & {:.2f} {:.2f} & {:.2f} {:.2f} ".format(
            names[field],
            np.mean(M31_illu_stats[field]), np.std(M31_illu_stats[field]),
            np.mean(MW_illu_stats[field]), np.mean(MW_illu_stats[field])))
    print()
             