#!/usr/bin/env python
# -*- coding: utf-8 -*-

from basico import *
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import pandas as pd
import pickle
from tqdm import tqdm


"""Create two kinds of plots: (a) comparisons of predicted time courses (model simulations) and measured time courses of validation experiments and (b) comparisons of the process performance with initial enzyme and substrate concentrations set to the initial concentrations of the baseline experiment and to those predicted by the selected optimization result."""


# DEFINE PLOTTING FUNCTIONS

# helper function to preprocess experimental data the same way as it was done for the parameter estimations
def preprocess_exp_data(exp_data_df, data_type='mean'):
    """ Change the following things for each data frame: (a) add up the measurements of adenosine tetra- and penta- phosphates (ATPP and ATPPP); (b) scale PolyP to get the amount of available phosphate bonds; (c) replace all zeros with 1e-12 to avoid incompatibilities with certain rate laws like the reversible Hill rate equation.

    :param exp_data_df: data frame of experimental data values (before preprocessing)
    :type exp_data_df: pandas.core.frame.DataFrame
    :param data_type: either 'mean' or 'SD' (this changes which equation is used to group ATPP and ATPP values)
    :type data_type: string

    :return exp_data_df: data frame of experimental data values (after preprocessing)
    :type exp_data_df: pandas.core.frame.DataFrame
    """

    if data_type == 'mean':
        # (a) add up the measurements of adenosine tetra- and penta- phosphates (ATPP, and ATPPP) in the ATPP column
        # mean values can just be added together to create the mean values of the new grouped
        exp_data_df['temp_ATPP_total'] = exp_data_df[['[ATPP]', '[ATPPP]']].sum(axis=1)
        # remove the old individual columns that were grouped
        exp_data_df = exp_data_df.drop(columns=['[ATPP]', '[ATPPP]'])
        # change the temporary name to the proper name that Copasi can recognize as dependent species
        exp_data_df.rename(columns={'temp_ATPP_total': '[ATPP]'}, inplace=True)
        # reorder the columns so that the new ATPP column is inserted after the ATP column
        exp_data_df = exp_data_df[['Time', '[Uri]', '[UMP]', '[UDP_GalNAc]', '[UDP]', '[UTP]', '[AMP]', '[ADP]', '[ATP]', '[ATPP]', '[GalNAc]_0', '[ATP]_0', '[ADP]_0', '[Uri]_0', '[E_PPA]_0', '[E_PPK3]_0', '[E_UDK]_0', '[E_NAHK]_0', '[E_GLMU]_0', '[E_UMPK]_0', '[PolyP]_0']]
        
        # (b) scale PolyP by *8 to get the amount of available phosphate bonds
        exp_data_df.at[0, '[PolyP]_0'] = exp_data_df.loc[0, '[PolyP]_0']*8

    elif data_type == 'SD':
        # (a) add up the measurements of adenosine tetra- and penta- phosphates (ATPP, and ATPPP) in the ATPP column
        # the standard deviations of the new group can be calculated as follows: sqrt(SD1^2 + SD2^2 + 2 * SD1 * SD2 * rho) where rho is the correlation coefficient; for now: assume that the two vriables are independent: then rho becomes 0 and we don't have to dig up the raw data points of all replicates to calculate it
        exp_data_df['temp_ATPP_total'] = np.sqrt(exp_data_df['[ATPP]']**2 + exp_data_df['[ATPPP]']**2)
        # remove the old individual columns that were grouped
        exp_data_df = exp_data_df.drop(columns=['[ATPP]', '[ATPPP]'])
        # change the temporary name to the proper name that Copasi can recognize as dependent species
        exp_data_df.rename(columns={'temp_ATPP_total': '[ATPP]'}, inplace=True)
        # reorder the columns so that the new ATPP column is inserted after the ATP column
        exp_data_df = exp_data_df[['Time', '[Uri]', '[UMP]', '[UDP_GalNAc]', '[UDP]', '[UTP]', '[AMP]', '[ADP]', '[ATP]', '[ATPP]']]

        # b) no PolyP column inside SD data frames so it does not need to be scaled

    # (c) replace all zeros with 1e-12 to avoid incompatibilities with certain rate laws like the reversible Hill rate equation
    for col in exp_data_df.columns:
        if col != 'Time':
            exp_data_df.replace({col: 0}, 1e-12, inplace=True)

    # return the modified data frame
    return exp_data_df

def visualize_predTC_withValData(model_path, fits_data_path, best_opt_result_path, exp_val_data_path, exp_val_data_SD_path, result_file_path):
    """ Visualize repeated time course simulations across an ensemble of parameter sets for a selected optimization result (prediction) together with the data of the associated validation experiment.

    :param model_path: path to the Copasi model object (.cps file)
    :type model_path: string
    :param fits_data_path: path to the parameter ensemble (.csv file)
    :type fits_data_path: string
    :param best_opt_result_path: path of the selected best optimization result (.pkl file)
    :type best_opt_result_path: string
    :param exp_val_data_path: path to the data (means over three replicates) of the validation experiment (.txt file)
    :type exp_val_data_path: string
    :param exp_val_data_SD_path: path to the data (standard deviation of the means) of the validation experiment (.txt file)
    :type exp_val_data_SD_path: string
    :param result_file_path: path for the generated files
    :type result_file_path: string
    """

    # LOAD DATA
    # ---------

    # load model file and get list of non-enzyme species names
    model_obj = load_model(model_path)
    non_E_species_names = [species_name for species_name in list(get_species(model=model_obj).index) if not species_name.startswith('E_')]

    # load the data of the parameter ensemble
    param_sets = pd.read_csv(fits_data_path)

    # load selected optimization result and get predicted initial concentrations
    with open(best_opt_result_path, 'rb') as pickle_file:
        opt_result = pickle.load(pickle_file)
    opt_result_pred_init_conc = opt_result['opt_result'][0]['sol']

    # load data of validation experiment (means and standard deviations)
    val_data_df_raw = pd.read_csv(exp_val_data_path, sep='\t')
    val_data_SD_df_raw = pd.read_csv(exp_val_data_SD_path, sep='\t')
    # preproces data of validation experiments: in the model, ATPP and ATPPP are grouped as one variable -> therefore, the data needs to be processed in the same way
    # 1) mean values can just be added together to create the mean values of the new grouped
    # 2) the standard deviations of the new group can be calculated as follows: sqrt(SD1^2 + SD2^2 + 2 * SD1 * SD2 * rho) where rho is the correlation coefficient; for now: assume that the two vriables are independent: then rho becomes 0 and we don't have to dig up the raw data points of all replicates to calculate it
    val_data_df = preprocess_exp_data(val_data_df_raw, data_type='mean')
    val_data_SD_df = preprocess_exp_data(val_data_SD_df_raw, data_type='SD')
    # rename columns so that they share the same format as the output data frames of Copasi time course simulations (remove square brackets around species names)
    val_data_col_rename_dict = dict()
    val_data_SD_col_rename_dict = dict()
    for col_name in val_data_df:
        if '[' and ']' in col_name:
            # use regex to find and remove the square brackets
            new_col_name = re.sub(r"[\[\]]", "", col_name)
            val_data_col_rename_dict.update({col_name: new_col_name})
    for col_name in val_data_SD_df:
        if '[' and ']' in col_name:
            # use regex to find and remove the square brackets
            new_col_name = re.sub(r"[\[\]]", "", col_name)
            val_data_SD_col_rename_dict.update({col_name: new_col_name})
    val_data_df.rename(columns=val_data_col_rename_dict, inplace=True)
    val_data_SD_df.rename(columns=val_data_SD_col_rename_dict, inplace=True)
    # re-index validation data frames (use 'Time' column as index)
    val_data_df.set_index('Time', inplace=True)
    val_data_SD_df.set_index('Time', inplace=True)


    # RUN SIMULATIONS
    # ---------------

    # only run the simulation if necessary
    if Path(f'{result_file_path}_sim_data.pkl').is_file():
        with open(f'{result_file_path}_sim_data.pkl', 'rb') as data_file:
            sim_tc_df_list = pickle.load(data_file)
            print(f'Simulation data taken from {result_file_path}_sim_data.pkl')
    else:
        # set  all species concentrations that were not used as optimization variables to their respective constant values
        set_species('ADP',          initial_concentration=1e-12, model=model_obj)
        set_species('AMP',          initial_concentration=1e-12, model=model_obj)
        set_species('ATPP',         initial_concentration=1e-12, model=model_obj)
        set_species('GalNAc1P',     initial_concentration=1e-12, model=model_obj)
        set_species('P',            initial_concentration=1e-12, model=model_obj)
        set_species('PolyP',        initial_concentration=256.0, model=model_obj)
        set_species('PP',           initial_concentration=1e-12, model=model_obj)
        set_species('UDP',          initial_concentration=1e-12, model=model_obj)
        set_species('UDP_GalNAc',   initial_concentration=1e-12, model=model_obj)
        set_species('UMP',          initial_concentration=1e-12, model=model_obj)
        set_species('UTP',          initial_concentration=1e-12, model=model_obj)
        if 'ATP' not in list(opt_result['opt_result'][0].index):
            # the b variant of the optimizations don't use ATP as optimization variable and instead fix it at 0.5 mM
            set_species('ATP',      initial_concentration=0.5,   model=model_obj)

        # loop over all i parameter sets and j optimization results
        sim_tc_df_list = list()
        for index, param_set in tqdm(param_sets.iterrows()):

            # 1) override model kinetic parameters with current parameter set
            # current list of parameter values
            curr_param_set_vals = list(param_set)
            # current list of parameter names
            curr_param_set_names = list(param_set.keys())
            # list() to only get content (values or names) of pandas series param_set
            # .keys() to get names of values in pandas series param_set
            # trim both lists (get rid of first entries = name and last entries = obj_val of fit)
            curr_param_set_vals = curr_param_set_vals[1:-1]
            curr_param_set_names = curr_param_set_names[1:-1]
            # set model parameter values to current estimated parameters from rand sampling result
            for i in range(len(list(curr_param_set_vals))): 
                set_reaction_parameters(name = curr_param_set_names[i], value = curr_param_set_vals[i])

            # 2) override model initial enzyme concentrations with those of the selected optimization result
            # first set_species needs to be called twice so that the species is actually set ... no idea why that is ...
            set_species('E_GLMU', initial_concentration=opt_result_pred_init_conc.iloc[0], model=model_obj)
            set_species('E_GLMU', initial_concentration=opt_result_pred_init_conc.iloc[0], model=model_obj)
            set_species('E_NAHK', initial_concentration=opt_result_pred_init_conc.iloc[1], model=model_obj)
            set_species('E_PPA',  initial_concentration=opt_result_pred_init_conc.iloc[2], model=model_obj)
            set_species('E_PPK3', initial_concentration=opt_result_pred_init_conc.iloc[3], model=model_obj)
            set_species('E_UDK',  initial_concentration=opt_result_pred_init_conc.iloc[4], model=model_obj)
            set_species('E_UMPK', initial_concentration=opt_result_pred_init_conc.iloc[5], model=model_obj)
            set_species('GalNAc', initial_concentration=opt_result_pred_init_conc.iloc[6], model=model_obj)
            if 'ATP' in list(list(opt_result['opt_result'][0].index)):
                # ATP is being used as an optimization variable so we set it here and then move to Uri
                set_species('ATP',    initial_concentration=opt_result_pred_init_conc.iloc[7], model=model_obj)
                set_species('Uri',    initial_concentration=opt_result_pred_init_conc.iloc[8], model=model_obj)
            else:
                # ATP is not being used as an optimization variable so we directly move to the next species (Uri)
                set_species('Uri',    initial_concentration=opt_result_pred_init_conc.iloc[7], model=model_obj)

            # 3) simulate the model
            sim_result = run_time_course(duration=24, automatic=False, stepsize=0.1, model=model_obj)
            # only keep non-enzyme species columns
            sim_result_species = sim_result.loc[:, non_E_species_names]
            # save result
            sim_tc_df_list.append(sim_result_species)

        # export the list of simulated time course data frames as python object
        open_file = open(f'{result_file_path}_sim_data.pkl', "wb")
        pickle.dump(sim_tc_df_list, open_file)
        open_file.close()


    # CALCULATE SIMULATION STATISTICS
    # -------------------------------

    # convert PolyP values back to concentration in mM (the PolyP variable in the model represents the amount of available phosphate bonds across all PolyP chains; we assume 8 usable phosphate bonds per chain so for simulating the model the initial PolyP concentration is multiplied by a factor of 8 in the pre-processing)
    for idx, df in enumerate(sim_tc_df_list):
        sim_tc_df_list[idx]['PolyP'] = df['PolyP'] / 8
    # get means and standard deviations across all data frames
    # 1) convert the list of DataFrames to a 3D NumPy array (shape: (n_time_points, n_species, n_ensemble_members))
    data_array = np.dstack([df.values for df in sim_tc_df_list])
    # 2) calculate mean and std along the ensemble axis (axis=2)
    mean_values = np.mean(data_array, axis=2)
    std_values = np.std(data_array, axis=2)
    # 3) convert back to DataFrames for easier plotting
    mean_df = pd.DataFrame(mean_values, columns=non_E_species_names, index=np.linspace(0, 24, data_array.shape[0]))
    std_df = pd.DataFrame(std_values, columns=non_E_species_names, index=np.linspace(0, 24, data_array.shape[0]))


    # VISUALIZE DATA
    # --------------

    # set up a grid of subplots (3 rows, 5 columns)
    fig, axs = plt.subplots(3, 5, figsize=(25, 18))

    # set global font size
    font_size = 20

    # flatten the axs array for easy iteration
    axs = axs.flatten()

    # add data and titles to subplots
    for ax_idx, ax_obj in enumerate(axs):
        if ax_idx < len(non_E_species_names):  # ensure the number of non-enzyme species is not exceeded
            species_name = non_E_species_names[ax_idx]
            # plot all simulation trajectories
            for df in sim_tc_df_list:
                ax_obj.plot(df.loc[:, species_name], color='tab:blue', alpha=0.1)
            # plot average simulation trajectory (with st. dev.)
            ax_obj.plot(mean_df[species_name], color='tab:red', label='Mean')
            ax_obj.fill_between(mean_df.index,
                                mean_df[species_name] - std_df[species_name],
                                mean_df[species_name] + std_df[species_name],
                                color='tab:red', alpha=0.2, label='±1 Std Dev')
            # plot validation data points with error bars
            if species_name in val_data_df.columns:
                ax_obj.errorbar(val_data_df.index,
                                val_data_df[species_name],
                                yerr=val_data_SD_df[species_name],
                                fmt='o',
                                color='black',
                                ecolor='black',
                                capsize=3,
                                label='Validation')
            # set subplot title to species name
            ax_obj.set_title(species_name, fontsize=font_size, loc='center')
            # set x-axis ticks
            ax_obj.set_xticks([0, 4, 8, 12, 16, 20, 24])
            ax_obj.tick_params(axis='x')
            # change the fontsize of the axis tick labels
            ax_obj.tick_params(axis='both', which='major', labelsize=font_size-4)
            ax_obj.tick_params(axis='both', which='minor', labelsize=font_size-4)

    # hide any unused subplots
    for ax_idx in range(len(non_E_species_names), len(axs)):
        axs[ax_idx].axis('off')

    # set overall x and y axis labels
    fig.supxlabel('Time [h]', y=0.05, fontsize=font_size)
    fig.supylabel('Concentration [mM]', x=0.07, fontsize=font_size)

    # add a legend for the entire figure
    fig.legend(handles=[plt.Line2D([0], [0], color='tab:blue', alpha=0.3, label='Simulated Trajectories'),
                        plt.Line2D([0], [0], color='tab:red', label='Mean Trajectory'),
                        plt.Line2D([0], [0], color='tab:red', alpha=0.2, lw=8, label='±1 Std Dev'),
                        plt.Line2D([0], [0], color='black', marker='o', markersize=8, linestyle='None', label='Validation Data')],
               loc='lower center',
               ncol=4,
               fontsize=font_size)

    # save the plot
    plt.savefig(f'{result_file_path}_plot.png', bbox_inches="tight", dpi=300)

def plot_base_opt_process_comparison(baseline_exp_data_path, best_opt_result_path, exp_val_data_path, titer_and_yield_timepoint, result_file_path):
    """Plot initial enzyme and substrate concentrations as well as various relevant performance indicators (enzyme load, product titer, yield) of the unoptimized baseline experiment and the validation experiment which was performed to test the best scoring batch optimization result.

    :param baseline_exp_data_path: path to the baseline experiment data
    :type baseline_exp_data_path: string
    :param best_opt_result_path: path to the best scoring optimization result
    :type best_opt_result_path: string
    :param exp_val_data_path: path to the validation experiment data
    :type exp_val_data_path: string
    :param titer_and_yield_timepoint: the time point where titers and yields are calculated
    :type titer_and_yield_timepoint: integer
    :param result_file_path: path of the output file
    :type result_file_path: string
    """

    # LOAD DATA

    # load selected optimization result and get start (baseline) and predicted (optimized) initial concentrations
    with open(best_opt_result_path, 'rb') as pickle_file:
        opt_result = pickle.load(pickle_file)
    opt_result_start_init_conc = opt_result['opt_result'][0]['start']
    opt_result_pred_init_conc = opt_result['opt_result'][0]['sol']

    # load data of baseline and validation experiments
    baseline_data_df_raw = pd.read_csv(baseline_exp_data_path, sep='\t')
    val_data_df_raw = pd.read_csv(exp_val_data_path, sep='\t')

    # PREPARE DATA

#    Uri_MW = 244.20         # [g/mol]
#    GalNAc_MW = 221.21      # [g/mol]
#    UDP_GalNAc_MW = 607.36  # [g/mol]
#    reaction_vol = 0.5      # [ml]

    # 1) get enzyme distributions (convert concentrations from mM to g/l)
    baseline_E_distr = dict(opt_result_start_init_conc.loc[[name for name in opt_result_start_init_conc.index if
                                                            name.startswith('E_')]])
    opt_E_distr = dict(opt_result_pred_init_conc.loc[[name for name in opt_result_pred_init_conc.index if
                                                      name.startswith('E_')]])
    E_distr_df = pd.DataFrame(zip(baseline_E_distr.keys(),
                                  baseline_E_distr.values(),
                                  opt_E_distr.values()),
                              index=baseline_E_distr.keys(),
                              columns=['enzyme', 'baseline', 'optimized'])
    # conversion factors are based on molecular weight in g/mol divided by 1000 mmol/mol -> [g/mol]/[mmol/mol] = [g/mol]*[mol/mmol] = [g/mmol]
    conversion_factors = [49190/1000,  # E_GLMU
                          39903/1000,  # E_NAHK
                          19313/1000,  # E_PPA
                          34740/1000,  # E_PPK3
                          24353/1000,  # E_UDK
                          22405/1000]  # E_UMPK
    E_distr_df_gram = E_distr_df.copy()
    E_distr_df_gram.loc[:, 'baseline'] = E_distr_df_gram.loc[:, 'baseline'] * np.array(conversion_factors)
    E_distr_df_gram.loc[:, 'optimized'] = E_distr_df_gram.loc[:, 'optimized'] * np.array(conversion_factors)

    # 2) calculate initial enzyme loads [g/l]
    baseline_E_tot = {'gram_enzyme_load': sum(E_distr_df_gram.loc[:, 'baseline'])}
    opt_E_tot = {'gram_enzyme_load': sum(E_distr_df_gram.loc[:, 'optimized'])}

    # 3) get initial substrate concentrations [mmol/l]
    if 'ATP' in list(opt_result['opt_result'][0].index):
        # ATP was used as optimization variable
        baseline_S_init_conc = dict(opt_result_start_init_conc.loc[[name for name in opt_result_start_init_conc.index if
                                    name in ['Uri', 'GalNAc', 'ATP']]])
        opt_S_init_conc = dict(opt_result_pred_init_conc.loc[[name for name in opt_result_pred_init_conc.index if
                               name in ['Uri', 'GalNAc', 'ATP']]])
        S_init_conc_df = pd.DataFrame(zip(baseline_S_init_conc.keys(),
                                          baseline_S_init_conc.values(),
                                          opt_S_init_conc.values()),
                                      index=baseline_S_init_conc.keys(),
                                      columns=['substrate', 'baseline', 'optimized'])
    else:
        # ATP was not used as optimization variable and instead fixed at 0.5 mM
        baseline_S_init_conc = dict(opt_result_start_init_conc.loc[[name for name in opt_result_start_init_conc.index if
                                    name in ['Uri', 'GalNAc']]])
        baseline_S_init_conc.update({'ATP': 0.5})
        opt_S_init_conc = dict(opt_result_pred_init_conc.loc[[name for name in opt_result_pred_init_conc.index if
                               name in ['Uri', 'GalNAc']]])
        opt_S_init_conc.update({'ATP': 0.5})
        S_init_conc_df = pd.DataFrame(zip(baseline_S_init_conc.keys(),
                                          baseline_S_init_conc.values(),
                                          opt_S_init_conc.values()),
                                      index=baseline_S_init_conc.keys(),
                                      columns=['substrate', 'baseline', 'optimized'])

    # 4) get (absolute) product titers at the specified process time point [mmol_product/l]
    baseline_prod_titer = baseline_data_df_raw.loc[baseline_data_df_raw['Time'] == titer_and_yield_timepoint, '[UDP_GalNAc]'].iloc[0]
    opt_prod_titer = val_data_df_raw.loc[val_data_df_raw['Time'] == titer_and_yield_timepoint, '[UDP_GalNAc]'].iloc[0]

    # 5) calculate (absolute) product yields [mmol product/l / mmol substrate/l = mmol product / mmol substrate]
    baseline_yield = (2*(baseline_prod_titer)) / (baseline_S_init_conc['Uri'] + baseline_S_init_conc['GalNAc'])
    opt_yield = (2*(opt_prod_titer)) / (opt_S_init_conc['Uri'] + opt_S_init_conc['GalNAc'])

    # VISUALIZE DATA
    fig, axes = plt.subplot_mosaic([['E_distr', 'E_distr', 'E_distr', 'E_distr', 'E_distr', 'E_tot',],
                                    ['S_main_distr', 'S_ATP_distr', 'titer_abs', 'yield_abs', 'titer_rel', 'yield_rel']],
                                   constrained_layout=True, figsize=(12, 7))
    colors = ['tab:gray',  # unoptimized baseline
              'tab:red']   # optimized validation exp.

    # set global font size
    font_size = 12

    # plot enzyme distributions
    E_distr_df_gram.plot(ax=axes['E_distr'], kind='bar', stacked=False,
                         color=colors)
    axes['E_distr'].set_xticklabels([name.split('E_')[1] for name in E_distr_df_gram['enzyme']], rotation=0, ha='center', fontsize=font_size)
    axes['E_distr'].set_ylabel('Concentration [g/l]', fontsize=font_size)
    axes['E_distr'].get_legend().remove()

    # plot initial gram enzyme loads [g/l]
    axes['E_tot'].bar(['baseline', 'optimized'],
                     [baseline_E_tot['gram_enzyme_load'], opt_E_tot['gram_enzyme_load']],
                     color=colors)
    axes['E_tot'].set_xticklabels([])
    axes['E_tot'].set_xticks([])
    axes['E_tot'].set_ylabel('Enzyme Load [g/l]', labelpad=7, fontsize=font_size)

    # plot initial substrate concentrations in separate plots to account for different value ranges; plot title and legends of the enzyme distribution subplot are sufficient - no need for additional titles and legends
    S_init_conc_df.loc[['Uri', 'GalNAc'], :].plot(ax=axes['S_main_distr'], kind='bar', stacked=False,
                                                  color=colors)
    axes['S_main_distr'].set_xticklabels(S_init_conc_df['substrate'].loc[['Uri', 'GalNAc']], rotation=0, ha='center', fontsize=font_size)
    axes['S_main_distr'].set_ylabel('Concentration [mM]', fontsize=font_size)
    axes['S_main_distr'].get_legend().remove()
    S_init_conc_df.loc[['ATP'], :].plot(ax=axes['S_ATP_distr'], kind='bar', stacked=False,
                                        color=colors)
    axes['S_ATP_distr'].set_xticklabels(S_init_conc_df['substrate'].loc[['ATP']], rotation=0, ha='center', fontsize=font_size)
    axes['S_ATP_distr'].set_ylabel('Concentration [mM]', fontsize=font_size)
    axes['S_ATP_distr'].get_legend().remove()

    # plot absolute values of product titers [mmol product/l] and product yields [mmol product/mmol_substrate]
    axes['titer_abs'].bar(['baseline_titer_abs', 'optimized_titer_abs'],
                [baseline_prod_titer, opt_prod_titer],
                color=colors)
    axes['titer_abs'].set_xticklabels([])
    axes['titer_abs'].set_xticks([])
    axes['titer_abs'].set_ylabel('$Titer_{' + str(titer_and_yield_timepoint) + 'h}$ [$mmol_{product}/l$]', labelpad=7, fontsize=font_size)
    axes['yield_abs'].bar(['baseline_yield_abs', 'optimized_yield_abs'],
                [baseline_yield, opt_yield],
                color=colors)
    axes['yield_abs'].set_xticklabels([])
    axes['yield_abs'].set_xticks([])
    axes['yield_abs'].set_ylabel('$Yield_{' + str(titer_and_yield_timepoint) + 'h}$ [$mmol_{product}/mmol_{substrate}$]', labelpad=7, fontsize=font_size)

    # only plot relative values of product titers and product yields if it makes sense, i.e., if the enzyme load is different between the baseline and the optimized scenario (which in turn will lead to a difference between absolute and relative titers/yields)
    if np.round(baseline_E_tot['gram_enzyme_load'], 2) != np.round(opt_E_tot['gram_enzyme_load'], 2):
        # relative titer: (absolute titer [mmol product/l] ) / enzyme load [g enzyme/l] => (mmol product)/(g enzyme)
        axes['titer_rel'].bar(['baseline_titer_rel', 'optimized_titer_rel'],
                    [baseline_prod_titer / baseline_E_tot['gram_enzyme_load'],
                     opt_prod_titer / opt_E_tot['gram_enzyme_load']],
                    color=colors)
        axes['titer_rel'].set_xticklabels([])
        axes['titer_rel'].set_xticks([])
        axes['titer_rel'].set_ylabel('$rel. Titer_{' + str(titer_and_yield_timepoint) + 'h}$ [$mmol_{product}/g_{enzyme}$]', labelpad=7, fontsize=font_size)
        # relative yield = absolute yield [mmol product/ mmol substrate] / enzyme load [g enzyme/l] => (mmol product)/(mmol substrat * (g enzyme/l))
        axes['yield_rel'].bar(['baseline_yield_rel', 'optimized_yield_rel'],
                    [baseline_yield / baseline_E_tot['gram_enzyme_load'], 
                    opt_yield / opt_E_tot['gram_enzyme_load']],
                    color=colors)
        axes['yield_rel'].set_xticklabels([])
        axes['yield_rel'].set_xticks([])
        axes['yield_rel'].set_ylabel('$rel. Yield_{' + str(titer_and_yield_timepoint) + 'h}$ [$mmol_{product}/(mmol_{substrate} \cdot (g_{enzyme}/l))$]', labelpad=7, fontsize=font_size-2)
    else:
        axes['titer_rel'].axis("off")
        axes['yield_rel'].axis("off")

    # save figure and close it to free memory
    fig.savefig(f'{result_file_path}.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    
    # also gather all data in a dictionary and store it in a text file
    data_dict = {'enzyme_init': E_distr_df_gram,
                 'gram_enzyme_loads': {'baseline': baseline_E_tot['gram_enzyme_load'], 
                                       'optimized': opt_E_tot['gram_enzyme_load']},
                 'substrates_init': S_init_conc_df,
                 'titers_abs': {'baseline': baseline_prod_titer, 
                                'optimized': opt_prod_titer},
                 'titers_rel': {'baseline': baseline_prod_titer / baseline_E_tot['gram_enzyme_load'],
                                'optimized': opt_prod_titer / opt_E_tot['gram_enzyme_load']},
                 'yields_abs': {'baseline': baseline_yield,
                                'optimized': opt_yield},
                 'yields_rel': {'baseline': baseline_yield / baseline_E_tot['gram_enzyme_load'],
                                'optimized': opt_yield / opt_E_tot['gram_enzyme_load']}}
    list_of_strings = [f'{key} : {data_dict[key]}' for key in data_dict]
    with open(f'{result_file_path}_data.txt', 'w') as file:
        [file.write(f'{string}\n') for string in list_of_strings]


# LOAD RELEVANT DATA AND CREATE PLOTS

# get path to the experimental data directory and the data of the baseline experiment Exp36_50
cwd = Path.cwd()
exp_val_data_dir = cwd / 'exp_val_data'
baseline_exp_data_path = f'{str(exp_val_data_dir)}\\UDP-GalNAc36_50_with_initConcColumns.txt'

# --------------------------------------------------------------------------------------------------------#

# Optimization Variant: a_constPolyP32mM_ubATP2.5mM
opt_var_a_path = cwd / 'a_constPolyP32mM_ubATP2.5mM'

# 1) Titer Maximization at 7 h s.t. Maximal Enzyme Load

# Model Variant: rep2_MV5, Optimization Result: OP01withS_7h[O11], Validation Experiment: Exp40_TO7h
model_path = f'{opt_var_a_path}\\rep2_MV5\\ImpExtSearch_v23b_rep2_MV5.cps'
fits_data_path = f'{opt_var_a_path}\\rep2_MV5\\sampling_output_Particle_Swarm50runs_ImpExtSearch_v23b_rep2_ModelVar5.csv'
best_opt_result_path = f'{opt_var_a_path}\\rep2_MV5\\ImpExtSearch_v23b_rep2_MV5_PartSwarm50x_OP01withS_7h_CrossVal_Scoring_Result_O11.pkl'
exp_val_data_path = f'{str(exp_val_data_dir)}\\UDP-GalNAc40_TO7h_with_initConcColumns.txt'
exp_val_data_SD_path = f'{str(exp_val_data_dir)}\\UDP-GalNAc40_TO7h_SD.txt'
tc_plot_result_file_path = f'{opt_var_a_path}\\rep2_MV5\\ImpExtSearch_v23b_rep2_MV5_OP01withS_7h_O11'
visualize_predTC_withValData(model_path, fits_data_path, best_opt_result_path, exp_val_data_path, exp_val_data_SD_path, tc_plot_result_file_path)
titer_and_yield_timepoint = 7  # [h]
perf_plot_result_file_path = f'{opt_var_a_path}\\rep2_MV5\\ImpExtSearch_v23b_rep2_MV5_OP01withS_7h_O11_base_opt_compare'
plot_base_opt_process_comparison(baseline_exp_data_path, best_opt_result_path, exp_val_data_path, titer_and_yield_timepoint, perf_plot_result_file_path)

# --------------------------------------------------------------------------------------------------------#

# Optimization Variant: b_constPolyP32mMATP0.5mM
opt_var_b_path = cwd / 'b_constPolyP32mMATP0.5mM'


# 2) Enzyme Load Minimization s.t. Minimal Titer and Yield

# Model Variant: rep2_MV5, Optimization Result: OP05withS[O19], Validation Experiment: Exp40_LO
model_path = f'{opt_var_b_path}\\rep2_MV5\\ImpExtSearch_v23b_rep2_MV5.cps'
fits_data_path = f'{opt_var_b_path}\\rep2_MV5\\sampling_output_Particle_Swarm50runs_ImpExtSearch_v23b_rep2_ModelVar5.csv'
best_opt_result_path = f'{opt_var_b_path}\\rep2_MV5\\ImpExtSearch_v23b_rep2_MV5_PartSwarm50x_OP05withS_CrossVal_Scoring_Result_O19.pkl'
exp_val_data_path = f'{str(exp_val_data_dir)}\\UDP-GalNAc40_LO_with_initConcColumns.txt'
exp_val_data_SD_path = f'{str(exp_val_data_dir)}\\UDP-GalNAc40_LO_SD.txt'
result_file_name = f'{opt_var_b_path}\\rep2_MV5\\ImpExtSearch_v23b_rep2_MV5_OP05withS_O19'
visualize_predTC_withValData(model_path, fits_data_path, best_opt_result_path, exp_val_data_path, exp_val_data_SD_path, result_file_name)
titer_and_yield_timepoint = 24  # [h]
perf_plot_result_file_path = f'{opt_var_b_path}\\rep2_MV5\\ImpExtSearch_v23b_rep2_MV5_OP05withS_O19_base_opt_compare'
plot_base_opt_process_comparison(baseline_exp_data_path, best_opt_result_path, exp_val_data_path, titer_and_yield_timepoint, perf_plot_result_file_path)
