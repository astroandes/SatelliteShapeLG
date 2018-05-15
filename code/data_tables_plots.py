import compile_randomized_data as crd

obs = False
illustris1 = False
illustris1dark = False
elvis = False

if obs:        
    print('Compiling stats for the observations')
    output_path = "../data/obs_summary/"
    for j in range(11,16):
        print('\t Nsat = {}'.format(j))
        crd.compile_stats(group_id=0, n_sat = j, n_random=10000, 
                        output_path=output_path, sort_column='vmag', randomize=False, reverse=False, obs_data=True)
    
if illustris1:
    print('Compiling stats for illustris1')
    input_path = "../data/illustris1_mstar_selected/"
    output_path = "../data/illustris1_mstar_selected_summary/"
    for j in range(11,16):
        print('\t Nsat = {}'.format(j))
        for i in range(27):
            crd.compile_stats(group_id=i, n_sat=j, n_random=10000, data_path=input_path,
                            output_path=output_path, sort_column='vmax', reverse=True, randomize=False)

if illustris1dark:
    print('Compiling stats for illustris1dark')
    input_path = "../data/illustris1dark_mstar_selected/"
    output_path = "../data/illustris1dark_mstar_selected_summary/"
    for j in range(11,16):
        print('\t Nsat = {}'.format(j))
        for i in range(27):
            crd.compile_stats(group_id=i, n_sat=j, n_random=10000, data_path=input_path,
                            output_path=output_path, sort_column='vmax', reverse=True, randomize=False)
if elvis:
    print('Compiling stats for ELVIS')
    input_path = "../data/elvis_mstar_selected/"
    output_path = "../data/elvis_mstar_selected_summary/"
    for j in range(11,16):
        print('\t Nsat = {}'.format(j))
        for i in range(12):
            crd.compile_stats(group_id=i, n_sat=j, n_random=10000, data_path=input_path,
                            output_path=output_path, sort_column='vmax', reverse=True, randomize=False, elvis=True)
            
import asphericity_stats as asp

ObsTable = False
SimTable = False
CovPlots = False
ObsAsphPlots = False
PrintModelNumbers = False
PlotModelNumbers = False
PlotShapeRandomObs = False
PlotShapeObsSim = False
PlotShapeObsSimNormed = False
PrintCovarianceTables = True

if ObsTable:
    asp.print_table_obs_shape()
    
if SimTable:
    asp.print_table_sim_shape()
    
if CovPlots:
    for nsat in range(11,12):
        asp.plot_covariance('illustris1dark', nsat)
        asp.plot_covariance('illustris1', nsat)
        asp.plot_covariance('elvis', nsat)
    
if ObsAsphPlots:
    asp.plot_asphericity_obs(0)
    asp.plot_asphericity_obs(1)
    asp.plot_asphericity_obs(2)
    
if PrintModelNumbers:
    asp.print_numbers()

if PlotModelNumbers:
    asp.plot_numbers()

if PlotShapeRandomObs:
    for nsat in range(11,16):
        asp.plot_shape_obs_randoms(nsat)
    
if PlotShapeObsSim:
    for nsat in range(11,12):
        asp.plot_shape_obs_sims('illustris1dark', nsat)
        asp.plot_shape_obs_sims('illustris1', nsat)
        asp.plot_shape_obs_sims('elvis', nsat)

if PlotShapeObsSimNormed:
    for nsat in range(11,12):
        asp.plot_shape_obs_sims_normed('illustris1dark', nsat)
        asp.plot_shape_obs_sims_normed('illustris1', nsat)
        asp.plot_shape_obs_sims_normed('elvis', nsat)
        
if PrintCovarianceTables:
    for nsat in range(11,16):
        asp.print_covariance('illustris1dark', nsat)
        asp.print_covariance('illustris1', nsat)
        asp.print_covariance('elvis', nsat)
