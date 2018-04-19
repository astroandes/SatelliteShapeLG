import numpy as np
import glob
import os
import matplotlib.pyplot as plt
import corner
import asphericity_stats as asp

def plot_randoms_and_sims(simulation, n_sat):
    print('simulation {}'.format(simulation))
    in_path = "../data/{}_mstar_selected_summary/".format(simulation)
    M31_sim_stats, MW_sim_stats = asp.load_experiment(input_path=in_path, n_sat=n_sat, full_data=True)

    in_path = "../data/obs_summary/"
    M31_obs_stats, MW_obs_stats = asp.load_experiment(input_path=in_path, n_sat=n_sat, full_data=False)
    M31_obs = asp.get_data_obs(M31_obs_stats, normed=False)
    MW_obs = asp.get_data_obs(MW_obs_stats, normed=False)

    M31_obs_stats, MW_obs_stats = asp.load_experiment(input_path=in_path, n_sat=n_sat, full_data=True)
    
        
    data_random_obs_M31 = np.array([M31_obs_stats['width_random'],
                                    M31_obs_stats['ca_ratio_random'], 
                                    M31_obs_stats['ba_ratio_random']]).T
    data_random_obs_MW = np.array([MW_obs_stats['width_random'],
                                    MW_obs_stats['ca_ratio_random'], 
                                    MW_obs_stats['ba_ratio_random']]).T
    
    print(np.shape(data_random_obs_MW))
    print(MW_obs)
        
    plt.figure(figsize=(8,5))
    plt.rc('text', usetex=True,)
    plt.rc('font', family='serif', size=25)
    figure = corner.corner(data_random_obs_M31, 
                      quantiles=[0.16, 0.5, 0.84],
                      labels=[r"$w$ M31", r"$c/a$ M31", r"$b/a$ M31"],
                      show_titles=True, title_kwargs={"fontsize": 12}, 
                      truths=M31_obs['data_obs'])
        
        
    min_w = 10; max_w = 90; min_ac = 0.0; max_ac = 1.0 ; min_ab = 0.6; max_ab = 1.0
    ndim = 3
    axes = np.array(figure.axes).reshape(ndim,ndim)
    ax = axes[1,1];ax.set_xlim(min_ac, max_ac)
    ax = axes[2,1];ax.set_xlim(min_ac, max_ac);ax.set_ylim(min_ab, max_ab)
    ax = axes[2,2];ax.set_xlim(min_ab, max_ab)
    ax = axes[0,0];ax.set_xlim(min_w, max_w)
    ax = axes[1,0];ax.set_xlim(min_w, max_w);ax.set_ylim(min_ac, max_ac)
    ax = axes[2,0];ax.set_xlim(min_w, max_w);ax.set_ylim(min_ab, max_ab)
    
    #    for i in range(ndim):
#        ax = axes[i, i]
#        ax.set_xlim(-5,5)
            
#    for yi in range(ndim):
#        for xi in range(yi):
#            ax = axes[yi, xi]
#            ax.set_xlim(-5,5)
#            ax.set_ylim(-5,5)

    filename = "../paper/input_{}_M31_n_{}.pdf".format(simulation, n_sat)
    print('saving figure to {}'.format(filename))
    plt.savefig(filename, bbox_inches='tight')
    plt.clf()
        
    plt.figure(figsize=(8,5))
    plt.rc('text', usetex=True,)
    plt.rc('font', family='serif', size=25)
    figure = corner.corner(data_random_obs_MW, 
                      quantiles=[0.16, 0.5, 0.84],
                      labels=[r"$w$ MW", r"$c/a$ MW", r"$b/a$ MW"],
                      show_titles=True, title_kwargs={"fontsize": 12}, 
                      truths=MW_obs['data_obs'])
    
    ndim = 3
    axes = np.array(figure.axes).reshape(ndim,ndim)
    ax = axes[1,1];ax.set_xlim(min_ac, max_ac)
    ax = axes[2,1];ax.set_xlim(min_ac, max_ac);ax.set_ylim(min_ab, max_ab)
    ax = axes[2,2];ax.set_xlim(min_ab, max_ab)
    ax = axes[0,0];ax.set_xlim(min_w, max_w)
    ax = axes[1,0];ax.set_xlim(min_w, max_w);ax.set_ylim(min_ac, max_ac)
    ax = axes[2,0];ax.set_xlim(min_w, max_w);ax.set_ylim(min_ab, max_ab)
#    axes = np.array(figure.axes).reshape(ndim,ndim)
#    for i in range(ndim):
#        ax = axes[i, i]
#        ax.set_xlim(-5,5)
            
#    for yi in range(ndim):
#        for xi in range(yi):
#            ax = axes[yi, xi]
#            ax.set_xlim(-5,5)
#            ax.set_ylim(-5,5)
    filename = "../paper/input_{}_MW_n_{}.pdf".format(simulation, n_sat)
    print('saving figure to {}'.format(filename))
    plt.savefig(filename, bbox_inches='tight')
    plt.clf()
    
    plt.close('all')


#plot_randoms_and_sims('illustris1', 11)


