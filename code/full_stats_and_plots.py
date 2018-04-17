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

def load_obs(obs_name):
    dtype=[('name','|S20'),
           ('x', 'f8'),
           ('y', 'f8'), 
           ('z', 'f8'),
           ('delta_plus', 'f8'),
           ('delta_minus', 'f8'),
           ('vmag', 'f8'),
           ('delta_vmag', 'f8')]
    obs_data = np.loadtxt(obs_name, dtype=dtype)    
    return obs_data

def main_halos(snap_data, sort_column='mstar', single_reverse=False):
    id_sorted = np.argsort(snap_data[sort_column])
    if single_reverse:
        id_halo_A = id_sorted[0]
        main_halos_data = snap_data[id_halo_A]
    else:
        id_halo_A = id_sorted[-1]
        id_halo_B = id_sorted[-2]
        main_halos_data = snap_data[[id_halo_A, id_halo_B]]
    return main_halos_data

def satellite_halos(halo_data, halo_center, radius=300, 
                    max_points=-1, max_initial_points=-1,
                    sort_column='vmag', reverse=False, randomize=False):
        
    jj = np.argsort(halo_data[sort_column])
    if reverse:
        jj = jj[::-1]
        
    tmp_halo = halo_data[jj]

    for col in list(['x','y','z']):
        tmp_halo[col] = tmp_halo[col] - halo_center[col]
    
    r = np.sqrt(tmp_halo['x']**2 + tmp_halo['y']**2 + tmp_halo['z']**2)
    ii = (r < radius) & (r>1E-6)
    
    tmp_halo = tmp_halo[ii]
    r = np.sqrt(tmp_halo['x']**2 + tmp_halo['y']**2 + tmp_halo['z']**2)
    
    if max_initial_points>0 :
        tmp_halo = tmp_halo[:max_initial_points]
        r = np.sqrt(tmp_halo['x']**2 + tmp_halo['y']**2 + tmp_halo['z']**2)

    if max_points > 0:
        if randomize:
            jj = np.argsort(np.random.random(len(tmp_halo)))
            return tmp_halo[jj[:max_points]], np.min(r[jj[:max_points]]), np.max(r[jj[:max_points]])
        else:
            return tmp_halo[:max_points], np.min(r[:max_points]), np.max(r[:max_points])
    else:
        return tmp_halo, np.min(r), np.max(r)

def gen_random_sphere(n_points):
    """
    Sets of points in the 3D sphere
    """
    r = np.random.random(n_points)**(1.0/3.0)
    phi = np.random.random(n_points) * 2.0 * np.pi
    costheta = 2.0*(np.random.random(n_points) -0.5)
    theta = np.arccos(costheta)
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(th
    
def compile_randomized(group_id=0, iter_id=0, n_sat=11, n_random=1000, 
                         elvis = False,
                         obs_data=False,
                         reverse=False,
                         sort_column='vmax',
                         randomize=False,
                         data_path = "../data/mstar_selected/", 
                         obs_data_path = "../data/obs/",
                         output_path = "../data/mstar_selected_summary/"):
    
    if obs_data:
        MW_data = load_obs(os.path.join(obs_data_path, "MW_satellites.txt"))
        M31_data = load_obs(os.path.join(obs_data_path, "M31_satellites.txt"))
    
        M31 = main_halos(M31_data, sort_column='vmag', single_reverse=True)
        MW = main_halos(MW_data, sort_column='vmag', single_reverse=True)
    else:
        if elvis:
            base_name = 'elvis'
        else:
            base_name = 'Illustris_group'
        M31_data = load_snapshot(os.path.join(data_path, "{}_{}.dat".format(base_name, group_id)), elvis=elvis)
        MW_data = load_snapshot(os.path.join(data_path, "{}_{}.dat".format(base_name, group_id)), elvis=elvis)
        LG_data = load_snapshot(os.path.join(data_path, "{}_{}.dat".format(base_name, group_id)), elvis=elvis)
        M31, MW = main_halos(LG_data, sort_column=sort_column, single_reverse=False)
   
        
    satellite_data_A, min_r_M31, max_r_M31 = satellite_halos(M31_data, M31,sort_column=sort_column)
    satellite_data_B, min_r_MW, max_r_MW = satellite_halos(MW_data, MW, sort_column=sort_column)
    
    N_A = len(satellite_data_A)
    N_B = len(satellite_data_B)
   

    if((N_A < n_sat_min) | (N_B <n_sat_min)):
        print('Failed Check Groupid, N bright:!', group_id, N_A, N_B)
        return
    
    satellite_data_A, min_r_M31, max_r_M31 = satellite_halos(M31_data,M31,sort_column=sort_column, 
                                                             max_points=n_sat, reverse=reverse,   randomize=randomize)
    satellite_data_B, min_r_MW, max_r_MW = satellite_halos(MW_data, MW, sort_column=sort_column,
                                                            max_points=n_sat, reverse=reverse, randomize=randomize)
 
    #print(satellite_data_B)
    #number of bright satellites
    N_A = len(satellite_data_A)
    N_B = len(satellite_data_B)
    #print('MW', N_B, satellite_data_B['name'])
    #print('M31', N_A, satellite_data_A['name'])
    
    output_A = open(os.path.join(output_path, 
                                 "M31_group_{}_nsat_{}_iter_{}.dat".format(group_id, n_sat_max, iter_id)), "w")
    output_B = open(os.path.join(output_path, 
                                 "MW_group_{}_nsat_{}_iter_{}.dat".format(group_id, n_sat_max, iter_id)), "w")
    #print(satellite_data_B)
    #minimum and maximum radius for the satellites
    output_A.write("{:2f} {:2f}\t".format(min_r_M31, max_r_M31))
    output_B.write("{:2f} {:2f}\t".format(min_r_MW, max_r_MW))

    r_AB = write_center_info(output_A, M31, MW)
    r_AB = write_center_info(output_B, M31, MW)
    write_inertia_plane(output_A, satellite_data_A, M31, unit_vector=r_AB, randomize=False)
    write_inertia_plane(output_B, satellite_data_B, MW, unit_vector=r_AB, randomize=False)
    output_A.write("\n")
    output_B.write("\n")
    
    for i in range(n_random):
        output_A.write("{:2f} {:2f}\t".format(min_r_M31, max_r_M31))
        output_B.write("{:2f} {:2f}\t".format(min_r_MW, max_r_MW))
        r_AB = write_center_info(output_A, M31, MW)
        r_AB = write_center_info(output_B, M31, MW)
        write_inertia_plane(output_A, satellite_data_A, M31, unit_vector=r_AB, randomize=True)
        write_inertia_plane(output_B, satellite_data_B, MW, unit_vector=r_AB, randomize=True)

        output_A.write("\n")
        output_B.write("\n")

    output_A.close()
    output_B.close()