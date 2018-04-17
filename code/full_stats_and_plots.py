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
                    max_points=-1, max_initial_points=-1,sort_column='vmag', reverse=False, randomize=False):
        
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