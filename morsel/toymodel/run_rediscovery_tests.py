#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from matplotlib.lines import Line2D
import time
import shutil
import sys

sys.path.append("..")
from func_lib import *


# IDEA

""" # Idea of the rediscovery test:

# 1) create a small and simple toy model

# 2) add some defined model changes to create the 'true' model version (that is what we're trying to rediscover)

# 3) simulate the true model structure with defined parameter values to generate 'fake experimental data'

# 4) run the whole model selection (NegCtrl Ensemble + ImpExtSearch + Final Selection) with the unchanged toymodel as start variant + a list of potential model changes that includes the ones that were actually added to create the gold standard model version; use the generated 'fake experimental data'

# 5) check if at the end the model variant that comes out of the final selection is close or even identical to the true model (if that is the case then the true standard model was successfully rediscovered from the data)
"""

""" Definition of a small and simple toymodel:
- substrate S, product P
- main reaction r1: S = P (catalyzed by an enzyme E which is not part of the model description)

Uncertainties to address:
1) measurements also show a compound X => where does it come from? maybe a side reaction also catalyzed by E?
2) does X act as a regulator of the main reaction r1? 

Library of potential model changes:
- 3 stoichiometric changes (different hypotheses for where X comes from):
    - r2: S = X
    - r3: P = X
    - r4: S + P = X
    => different potential unknown side reactions catalyzed by enzyme E
- 2 kinetic changes (different hypotheses for how X affects r1):
    - r1.X.act
    - r1.X_inhib

Assumptions:
- all reactions are reversible
- all reaction can adequatly be described by mass action rate laws
- for the sake of simplicity's any competitive effects (between metabolites that bind to the same enzyme) are ignored

Model units:
- time: [s]
- volume: [l]
- substance amount: [mmol]
=> substance concentration: [mmol/l]
"""


# SETUP

# set a random seed for reproducibility
np.random.seed(42)

def calc_and_vis_rediscovery_test(noise_std, vari_terms):
    """Set up and run the rediscovery test using the toy model which is defined in the function library.
    
    :param noise_std: standard deviation of the Gaussian noise that is added to the data generated from simulating the true model (float)
    :param vari_terms: list of candidate model changes (each model change is represented by a list of tuples that specify where and how the model structure should be changed)
    """

    # SETUP

    num_candidates = len(vari_terms) 
    test_scenario = f'noise_std_{noise_std}_num_candidates_{num_candidates}'
    # set model name
    model_name = f'toymodel_{test_scenario}'
    # define fit setup for all parameters that are part of the NegCtrl model structure (all mass action rate constants)
    fit_params_info = {'(r1).k_f': {'start': 1, 'lower': 1e-1, 'upper': 1e1},
                       '(r1).k_r': {'start': 1, 'lower': 1e-1, 'upper': 1e1},
                       '(r2).k_f': {'start': 1, 'lower': 1e-1, 'upper': 1e1},
                       '(r2).k_r': {'start': 1, 'lower': 1e-1, 'upper': 1e1},
                       '(r3).k_f': {'start': 1, 'lower': 1e-1, 'upper': 1e1},
                       '(r3).k_r': {'start': 1, 'lower': 1e-1, 'upper': 1e1},
                       '(r4).k_f': {'start': 1, 'lower': 1e-1, 'upper': 1e1},
                       '(r4).k_r': {'start': 1, 'lower': 1e-1, 'upper': 1e1}}
    # define global settings for all parameter estimations
    PE_algorithm_info = {'name': 'Particle Swarm', 
                     'settings': {'method': {'Iteration Limit': 1000,
                                             'Swarm Size': 50,
                                             'Std. Deviation': 1e-06}}}
    # set how many parameter sets are calculated for each ensemble
    n_runs = 20
    # set how often the imp. ext. search is repeated
    morsel_runs = 5

    # ------------------------------------------------------------------------------------------------------ #

    # DEFINE TRUE MODEL STRUCTURE, REGULATIONS AND KINETIC PARAMETERS
    # create a specific model variant ('true model') and simulate it with specific kinetic parameter values

    #           base term1
    true_var = [1,   6,    # r1: S = P with r1.X_inhib (ID: 6)
                0,   0,    # r2: S = X
                0,   0,    # r3: P = X
                1,   0]    # r4: S + P = X
    true_var = pd.DataFrame(np.array(true_var).reshape(4, 2),
                            columns=['base', 'term1'],
                            index=['r1', 'r2', 'r3', 'r4'])
    true_model_object, true_reaction_schemes_dict = create_model_obj_from_struct_var(true_var, model_data=toymodel_data, term_libs=toymodel_term_libs)
    true_params = {'(r1).k_f': 1,
                   '(r1).k_r': 1,
                   '(r1).ki_X': 1,
                   '(r4).k_f': 1,
                   '(r4).k_r': 1}
    init_conc_dict = {'S': 10,
                      'P': 0,
                      'X': 0}
    true_tc_result = sim_model_obj_with_estim_params(true_model_object, true_params, init_conc_dict, duration=20, step_number=10)
    true_tc_result.columns = ['[P]', '[S]', '[X]']
    # use time cours result of the true model simulation to generate 'fake experimental data' for the model selection
    # generate Gaussian noise for the selected columns (excluding the first row) and store it in a DataFrame
    noisy_tc_result = true_tc_result.copy(deep=True)
    # rename columns and add extra columns for initial values
    noisy_tc_result['[P]_0'] = pd.Series([0])
    noisy_tc_result['[S]_0'] = pd.Series([10])
    noisy_tc_result['[X]_0'] = pd.Series([0])
    noise_df = pd.DataFrame(np.random.normal(0, noise_std, size=noisy_tc_result.loc[2.0:, ['[P]', '[S]', '[X]']].shape),
                            index=noisy_tc_result.loc[2.0:].index,
                            columns=['[P]', '[S]', '[X]'])
    # add the noise to the copied DataFrame (excluding the first row)
    noisy_tc_result.loc[2.0:, ['[P]', '[S]', '[X]']] += noise_df
    # export the modified data
    noisy_tc_result.to_csv(f'true_model_data_noise_std_{noise_std}.txt', sep='\t')

    # plot data generated from true model (with and without noise)
    plot_colors = {'[S]': '#17becf',  # c
                   '[P]': '#e377c2',  # m
                   '[X]': '#bcbd22'}  # y
    fig, ax = plt.subplots(figsize=(4,3))
    for species_name, line_color in plot_colors.items():
        ax.plot(noisy_tc_result.index, noisy_tc_result.loc[:, f'{species_name}'], 
                label=species_name, marker='.', linestyle='dotted', linewidth=1, color=line_color)
    ax.set(xlabel='Time [h]', ylabel='Concentration [mM]')
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                 ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(12)
    # create custom legend
    species_handles = [Line2D([], [], color=plot_colors['[S]'], marker='.', linestyle='-', linewidth=1, label='S'),
                       Line2D([], [], color=plot_colors['[P]'], marker='.', linestyle='-', linewidth=1, label='P'),
                       Line2D([], [], color=plot_colors['[X]'], marker='.', linestyle='-', linewidth=1, label='X')]
    all_handles = species_handles
    all_labels = [h.get_label() for h in all_handles]
    ax.legend(handles=all_handles, labels=all_labels, loc='upper right', fontsize=12)
    fig.savefig(f'true_model_tc_noise_std_{noise_std}.png', dpi=200, bbox_inches='tight')
    plt.close()

    # ------------------------------------------------------------------------------------------------------ #

    # NEGATIVE CONTROL VARIANT PARAMETER ENSEMBLE
    # -------------------------------------------
    # Idea: as all methods in this script file are fundamentally based on testing the inclusion of various model modifications (e.g., regulation terms) it is necessary to first establish a baseline as a negative control: it only includes reactions that are part of the cascade design (and implements them with basic convenience kinetics rate laws without any extra modifications and/or regulations)
    # define negative control baseline variant
    # set model name
    neg_ctrl_model_name = 'NegCtrl'
    # calculate a negative control ensemble
    #               base
    neg_ctrl_var = [1,    # r1: S = P
                    0,    # r2: S = X
                    0,    # r3: P = X
                    0]    # r4: S + P = X
    neg_ctrl_var = pd.DataFrame(np.array(neg_ctrl_var).reshape(4, 1),
                                columns=['base'],
                                index=['r1', 'r2', 'r3', 'r4'])
    neg_ctrl_model_object, neg_ctrl_reaction_schemes_dict = create_model_obj_from_struct_var(neg_ctrl_var, model_data=toymodel_data, term_libs=toymodel_term_libs)
    # load experimental data
    exp_data_file_names = [f'true_model_data_noise_std_{noise_std}.txt']
    exp_data_dataframes = [pd.read_csv(file_name, sep='\t') for file_name in exp_data_file_names]
    # repeat the evaluation of the same model variant n times (n replicates) and store the results
    evaluated_vars_log = create_parameter_ensemble(neg_ctrl_var, n_runs, exp_data_file_names, exp_data_dataframes, neg_ctrl_model_name, 
                                                   fit_params_info=fit_params_info, model_data=toymodel_data, term_libs=toymodel_term_libs)
    # remove 'dump' text files that are created by the 'add_experiment' function that is called by the eval_struct_var function; windows and unix systems use different file path separators for name in exp_data_names
    clean_up_exp_dump_files(exp_data_file_names)

    # ------------------------------------------------------------------------------------------------------ #

    # IMPROVED EXTENSION SEARCH
    # -------------------------
    # Idea: similar to reduction search but with two differences:
    #       (1) 'extension', meaning that the search starts with the negative control and then always adds the 
    #                        best terms instead of removing them (to identify a solution that is as sparse as 
    #                        possible) and
    #       (2) 'improved',  meaning that the estimated model parameters of the best replicate of the term 
    #                        that is selected for permanent addition are saved and used as the starting 
    #                        parameters of the estimations in the next iteration of the overall search (this 
    #                        creates an upper bound for the RSS making sure that the squared errors don't get 
    #                        worse in larger models compared to smaller ones since any fit that is found for 
    #                        a smaller model also needs to be possible in a larger model)

    # define start variant: negative control
    #            base term1
    start_var = [1,   0,    # r1: S = P
                 0,   0,    # r2: S = X
                 0,   0,    # r3: P = X
                 0,   0]    # r4: S + P = X
    start_var = pd.DataFrame(np.array(start_var).reshape(4, 2),
                             columns=['base', 'term1'],
                             index=['r1', 'r2', 'r3', 'r4'])
    # load experimental data
    exp_data_file_names = [f'true_model_data_noise_std_{noise_std}.txt']
    exp_data_dataframes = [pd.read_csv(file_name, sep='\t') for file_name in exp_data_file_names]
    # define which model selection criteria are to be used
    selected_MSC = ['AIC', 'AICc', 'BIC', 'CIC1', 'CIC2', 'CIC3']
    # use the existing parameter ensemble of the initial core model (negative control; here the start variant) instead of evaluating the start variant inside the improved_extension_search function
    with open(f'NegCtrl_{n_runs}runs_evaluated_vars_log', 'rb') as f:
        start_var_ensemble = pickle.load(f)

    # run improved extension search and store result (output dictionary that contains list of constructed 
    # model variants and analysis result of the final model selection step)

    # repeat the imp. ext. search n times
    for runID in range(1,morsel_runs+1):

        output_dict = improved_extension_search(start_var, vari_terms,
                                                exp_data_file_names, exp_data_dataframes,
                                                PE_algorithm_info=PE_algorithm_info,
                                                selected_MSC=selected_MSC,
                                                start_var_ensemble=start_var_ensemble,
                                                model_data=toymodel_data,
                                                term_libs=toymodel_term_libs)

        # store result as pickle file
        model_selection_result_name = f'ToyModelReDisTest_{test_scenario}_run{runID}'
        storage_file = open(f'{model_selection_result_name}_output_dict', 'wb')
        pickle.dump(output_dict, storage_file)
        storage_file.close()

        # remove 'dump' text files that are created by the 'add_experiment' function that is called by the 
        # eval_struct_var function; windows and unix systems use different file path separators for name in 
        # exp_data_names
        clean_up_exp_dump_files(exp_data_file_names)

    # ------------------------------------------------------------------------------------------------------ #

    # FINAL SELECTION AND RESULT ENSEMBLE

    for runID in range(1,morsel_runs+1):

        # load an ImpExtSearch result
        model_selection_result_name = f'ToyModelReDisTest_{test_scenario}_run{runID}'
        reading_file = open(f'{model_selection_result_name}_output_dict', 'rb')
        output_dict = pickle.load(reading_file)
        reading_file.close()

        # analyze the output log
        analysis_result = output_dict['analysis_result']
        # get the index of the model variant with the best median rank and the fewest additions
        best_median_rank_model_idx = analysis_result["fitness_dataframe"]['median_rank'].idxmin()
        # get the index of the model variant with the best mean rank and the fewest additions
        best_mean_rank_model_idx = analysis_result["fitness_dataframe"]['mean_rank'].idxmin()
        # compare both and select the model with the lowest index (= lowest number of additions)
        selected_idx = 0
        if best_median_rank_model_idx < best_mean_rank_model_idx:
            best_ranking_var_log_dict = output_dict['evaluated_vars_log'][best_median_rank_model_idx]
            print(f"Model variant {best_median_rank_model_idx} with best median rank {analysis_result['fitness_dataframe'].iloc[best_median_rank_model_idx, :]['median_rank']} was selected.")
            selected_idx = best_median_rank_model_idx
        elif best_median_rank_model_idx > best_mean_rank_model_idx:
            best_ranking_var_log_dict = output_dict['evaluated_vars_log'][best_mean_rank_model_idx]
            print(f"Model variant {best_mean_rank_model_idx} with best mean rank {analysis_result['fitness_dataframe'].iloc[best_mean_rank_model_idx, :]['mean_rank']} was selected.")
            selected_idx = best_mean_rank_model_idx
        elif best_median_rank_model_idx == best_mean_rank_model_idx:
            # both indices are the same so they both point to the same model variant -> therefore it doesn't matter which one is used to look it up in the ImpExtSearch result
            best_ranking_var_log_dict = output_dict['evaluated_vars_log'][best_median_rank_model_idx]
            print(f"Model variant {best_median_rank_model_idx} with best median rank {analysis_result['fitness_dataframe'].iloc[best_median_rank_model_idx, :]['median_rank']} and best mean rank {analysis_result['fitness_dataframe'].iloc[best_mean_rank_model_idx, :]['mean_rank']} was selected.")
            selected_idx = best_median_rank_model_idx

        # create parameter ensemble of the selected model variant 
        model_variant_name = f'{model_selection_result_name}_ModelVar{selected_idx}'
        evaluated_vars_log = create_parameter_ensemble(best_ranking_var_log_dict['variant'], n_runs, exp_data_file_names, exp_data_dataframes, model_variant_name,
                                                       fit_params_info=fit_params_info, model_data=toymodel_data, term_libs=toymodel_term_libs)
        # get best parameter set of the ensemble
        final_ensemble = pd.read_csv(f'sampling_output_Particle_Swarm{n_runs}runs_{model_selection_result_name}_ModelVar{selected_idx}.csv', sep=',')
        final_ensemble_best_params = final_ensemble.loc[final_ensemble['obj'] == final_ensemble['obj'].min()]

        # save result log
        result_log = {'model_selection': model_selection_result_name,
                      'runID': runID,
                      'selected_variant_name': model_variant_name,
                      'selected_variant_idx': selected_idx,
                      'full_analysis_result': analysis_result,
                      'best_ranking_model_variant_idx': selected_idx,
                      'best_ranking_model_variant_structure': best_ranking_var_log_dict['variant'],
                      'full_param_ensemble': final_ensemble,
                      'best_param_set': final_ensemble_best_params}
        with open(f'rediscovery_test_run{runID}_result.pkl', 'wb') as file:
            pickle.dump(result_log, file)
        with open(f'rediscovery_test_run{runID}_result.txt', 'w') as file:
            for key,value in result_log.items():
                file.write("\n%s:\n%s\n" % (key,value))

    # select the best result across all MoRSel runs
    result_logs_list = []
    for runID in range(1, morsel_runs+1):
        with open(f'rediscovery_test_run{runID}_result.pkl', 'rb') as file:
            result_log = pickle.load(file)
        result_logs_list.append(result_log)

    # compare the structure of the best ranking variants of each run to the true model structure and return the index of the one with the fewest differences; if multiple model structures have the same minimum number of differences, return all indices
    best_ranking_vars_struct_diff_values = [log['best_ranking_model_variant_structure'].compare(true_var).size for log in result_logs_list]
    lowest_diff_vals_indices = [idx for idx, diff_val in enumerate(best_ranking_vars_struct_diff_values) if diff_val == min(best_ranking_vars_struct_diff_values)]
    # select the model structure that is the best match (that got closest to the true structure)
    struct_vars_idx_best_match = None
    if len(lowest_diff_vals_indices) == 1:
        # if only one model structure reaches the lowest number of differences compared to the true structure then it is the closest match -> select that one
        struct_vars_idx_best_match = lowest_diff_vals_indices[0]
    else:
        # if multiple model structures reach the lowest number of differences compared to the true structure then the best one is determined based on how close the estimated parameters are to the true parameters (using the objective value of the fit as a proxy); implementation: create a list of tuples (idx, obj_val), then find the tuple with the smallest objective value and extract its index
        best_tuple = min([(idx, result_logs_list[idx]['best_param_set']['obj'].iloc[0]) for idx in lowest_diff_vals_indices], key=lambda x: x[1])
        struct_vars_idx_best_match = best_tuple[0]

    # save all information of the overall best result across all MoRSel runs
    overall_best_result_log = {'model_selection': result_logs_list[struct_vars_idx_best_match]['model_selection'],
                               'runID': struct_vars_idx_best_match+1,  # to mach the way runs are counted (starting from 1)
                               'selected_variant_name': result_logs_list[struct_vars_idx_best_match]['selected_variant_name'],
                               'selected_variant_idx': result_logs_list[struct_vars_idx_best_match]['selected_variant_idx'],
                               'full_analysis_result': result_logs_list[struct_vars_idx_best_match]['full_analysis_result'],
                               'best_ranking_model_variant_idx': result_logs_list[struct_vars_idx_best_match]['best_ranking_model_variant_idx'],
                               'best_ranking_model_variant_structure': result_logs_list[struct_vars_idx_best_match]['best_ranking_model_variant_structure'],
                               'full_param_ensemble': result_logs_list[struct_vars_idx_best_match]['full_param_ensemble'],
                               'best_param_set': result_logs_list[struct_vars_idx_best_match]['best_param_set']}
    with open('rediscovery_test_overall_best_result.pkl', 'wb') as file:
        pickle.dump(overall_best_result_log, file)
    with open('rediscovery_test_overall_best_result.txt', 'w') as file:
        for key,value in overall_best_result_log.items():
            file.write("\n%s:\n%s\n" % (key,value))

    # ------------------------------------------------------------------------------------------------------ #

    # VISUALIZATION

    # 0) load the overall best result log (across all MoRSel runs) and set all variables that are requird for the visualization
    with open('rediscovery_test_overall_best_result.pkl', 'rb') as file:
        overall_best_result_log = pickle.load(file)
    model_selection_result_name = overall_best_result_log['model_selection']
    runID = overall_best_result_log['runID']
    selected_idx = overall_best_result_log['selected_variant_idx']
    final_ensemble = overall_best_result_log['full_param_ensemble']
    final_ensemble_best_params = overall_best_result_log['best_param_set']

    # 1) visualize parameter distributions of the selected best ranking model ensemble
    num_subplots = len(final_ensemble.columns[1:-1])
    num_rows = len(final_ensemble.columns[1:-1])
    num_cols = 1
    figsize_height = num_rows
    figsize_width = int(figsize_height*0.75)
    fig, axes = plt.subplots(num_rows, num_cols, figsize=(figsize_width, figsize_height))
    for col, ax in zip(final_ensemble.columns[1:-1], axes.flat):
        # create histogram
        sns.histplot(final_ensemble[col], ax=ax)
        # plot vertical line for the parameter value of the best parameter set
        ax.axvline(final_ensemble_best_params[col].iloc[0], color='tab:red')
        # plot vertical line for the true parameter value (if the parameter is part of the selected model structure)
        if col in true_params.keys():
            ax.axvline(true_params[col], color='black', linestyle='dashed')
        # set x and y label font sizes
        for item in ([ax.xaxis.label, ax.yaxis.label]):
            item.set_fontsize(14)
        # set x and y tick label font
        for item in (ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(10)
        # hide y-axis label
        ax.set_ylabel('')
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(f'final_ensemble_param_dist_run{runID}.png', bbox_inches='tight', dpi=200)

    # 2) plot time courses of the selected best ranking model ensemble and compare it to the time courses of the NegCtrl model ensemble
    with open(f'{neg_ctrl_model_name}_{n_runs}runs_evaluated_vars_log', 'rb') as f:
        neg_ctrl_log = pickle.load(f)
    with open(f'{model_selection_result_name}_ModelVar{selected_idx}_{n_runs}runs_evaluated_vars_log', 'rb') as f:
        selected_model_var = pickle.load(f)
    var1 = neg_ctrl_log[0]
    var2 = selected_model_var[0]
    plot_name = f'NegCtrl_{model_name}'
    exp_data_file_names = [f'true_model_data_noise_std_{noise_std}.txt']
    exp_data_dataframes = [pd.read_csv(file_name, sep='\t') for file_name in exp_data_file_names]
    exp_name = f'true_model_data_noise_std_{noise_std}.txt'
    exp_ID = 0  # set it to zero as there is only one exp. data file
    # create lists of dictionaries containing the initial concentrations of all experiments (keys: species names, values: initial concentrations)
    exp_init_conc_dicts = list()
    for df in exp_data_dataframes:
        # each data frame has special columns that contain the initial concentrations in row 0 for all relevant species, and they can be identified by their '_0' suffix in the column name
        init_conc_series = df.filter(regex="_0").iloc[0, :]
        # convert series into dictionary
        init_conc_dict = dict()
        for (idx_string, val) in init_conc_series.items():
            # remove 'Exp_' prefix, '_0' suffix and square brackets from index string
            idx_string = idx_string.replace('Exp_', '')
            idx_string = idx_string.replace('_0', '')
            idx_string = idx_string.replace('[', '')
            idx_string = idx_string.replace(']', '')
            # add a new entry to the dictionary using the index string as key and the initial concentration as value
            init_conc_dict.update({idx_string: val})
        # append dictionary to list
        exp_init_conc_dicts.append(init_conc_dict)
    # VARIANT 1 SIMULATION
    # create COPASI model object based on variant 1
    var1_model_obj, _ = create_model_obj_from_struct_var(var1['variant'], model_data=toymodel_data, term_libs=toymodel_term_libs)
    # run time course simulation with model object set to all sets of estimated parameters
    var1_TC_results_list = []
    PE_replicate_idx = 0
    for PE_replicate in tqdm(var1['estimation_results']):
        PE_replicate_tc = sim_model_obj_with_estim_params(var1_model_obj, PE_replicate['estimated_parameters'],
                                                          exp_init_conc_dicts[exp_ID], duration=20)
        # add a 'Time' column (generated from data frame index)
        PE_replicate_tc['Time'] = PE_replicate_tc.index
        # add a 'PE_replicate' column
        PE_replicate_tc['PE_replicate'] = [PE_replicate_idx for i in range(len(list(PE_replicate_tc.index)))]
        # remove 'Values[X]' columns (global quantities that I don't need here)
        PE_replicate_tc.drop(list(PE_replicate_tc.filter(regex='Values')),
                             axis=1, inplace=True)
        # collect simulation data frame
        var1_TC_results_list.append(PE_replicate_tc)
        PE_replicate_idx = PE_replicate_idx + 1
    # concatenate all data frames
    var1_TC_results_merged_df = pd.concat(var1_TC_results_list)
    # VARIANT 2 SIMULATION
    # create COPASI model object based on variant 2
    var2_model_obj, _ = create_model_obj_from_struct_var(var2['variant'], model_data=toymodel_data, term_libs=toymodel_term_libs)
    # run time course simulation with model object set to all sets of estimated parameters
    var2_TC_results_list = []
    PE_replicate_idx = 0
    for PE_replicate in tqdm(var2['estimation_results']):
        PE_replicate_tc = sim_model_obj_with_estim_params(var2_model_obj, PE_replicate['estimated_parameters'],
                                                          exp_init_conc_dicts[exp_ID], duration=20)
        # add a 'Time' column (generated from data frame index)
        PE_replicate_tc['Time'] = PE_replicate_tc.index
        # add a 'PE_replicate' column
        PE_replicate_tc['PE_replicate'] = [PE_replicate_idx for i in range(len(list(PE_replicate_tc.index)))]
        # remove 'Values[X]' columns (global quantities that I don't need here)
        PE_replicate_tc.drop(list(PE_replicate_tc.filter(regex='Values')),
                             axis=1, inplace=True)
        # collect simulation data frame
        var2_TC_results_list.append(PE_replicate_tc)
        PE_replicate_idx = PE_replicate_idx + 1
    # concatenate all data frames
    var2_TC_results_merged_df = pd.concat(var2_TC_results_list)
    # MULTIPLOT
    # create seaborn line plots from simulated result data frame for variant 1 (multi plot row index 0, line style dashed) and variant 2 (multi plot row index 0, line style solid)
    fig, ax = plt.subplots(1, 1, figsize=(4,3), layout='constrained')
    scatterplot_marker_size = 20
    colors = ['#e377c2',  # P: m
              '#17becf',  # S: c
              '#bcbd22']  # X: y
    # 1) plots of variant 1 (multiplot row index 0); reshape data frame from wide to long format for plotting of simulated results with error bands
    var1_TC_results_list_merged_longform = pd.DataFrame(columns=['PE_replicate', 'Time', 'Species', 'Concentration'])
    var1_TC_results_list_merged_longform['PE_replicate'] = pd.concat(
        [pd.Series(var1_TC_results_merged_df['PE_replicate']).repeat(len(var1_TC_results_merged_df.columns) - 2)],
        axis=0)
    var1_TC_results_list_merged_longform = var1_TC_results_list_merged_longform.reset_index(drop=True)
    # since the merged data frame in wide format has two extra columns ('Time' and 'PE_replicate'): repeat(n-2), we only want the length of the species columns
    var1_TC_results_list_merged_longform['Time'] = pd.concat(
        [pd.Series(var1_TC_results_merged_df.index).repeat(len(var1_TC_results_merged_df.columns) - 2)], axis=0,
        ignore_index=True)
    var1_TC_results_list_merged_longform['Species'] = pd.concat(
        [pd.Series(list(var1_TC_results_merged_df.columns[0:-2]) * len(var1_TC_results_merged_df.index))], axis=0,
        ignore_index=True)
    var1_TC_results_list_merged_longform['Concentration'] = pd.concat([var1_TC_results_merged_df.iloc[i, 0:-2]
                                                                       for i in
                                                                       range(len(var1_TC_results_merged_df.index))],
                                                                      axis=0, ignore_index=True)
    # plot the mean and 95% confidence interval by aggregating over PE replicates (at each time point)
    sns.lineplot(data=var1_TC_results_list_merged_longform[var1_TC_results_list_merged_longform.Species.isin(['P', 'S', 'X'])],
                 x="Time", y="Concentration", hue="Species", ax=ax, linestyle='dashed', 
                 palette=[sns.color_palette(colors)[0], sns.color_palette(colors)[1], sns.color_palette(colors)[2]])
    # 2) plots of variant 2 (same row as variant 1 [multiplot row index 0] but different line style); reshape data frame from wide to long format for plotting of simulated results with error bands
    var2_TC_results_list_merged_longform = pd.DataFrame(columns=['PE_replicate', 'Time', 'Species', 'Concentration'])
    var2_TC_results_list_merged_longform['PE_replicate'] = pd.concat(
        [pd.Series(var2_TC_results_merged_df['PE_replicate']).repeat(len(var2_TC_results_merged_df.columns) - 2)],
        axis=0)
    var2_TC_results_list_merged_longform = var2_TC_results_list_merged_longform.reset_index(drop=True)
    # since the merged data frame in wide format has two extra columns ('Time' and 'PE_replicate'): repeat(n-2), we only want the length of the species columns
    var2_TC_results_list_merged_longform['Time'] = pd.concat(
        [pd.Series(var2_TC_results_merged_df.index).repeat(len(var2_TC_results_merged_df.columns) - 2)], axis=0,
        ignore_index=True)
    var2_TC_results_list_merged_longform['Species'] = pd.concat(
        [pd.Series(list(var2_TC_results_merged_df.columns[0:-2]) * len(var2_TC_results_merged_df.index))], axis=0,
        ignore_index=True)
    var2_TC_results_list_merged_longform['Concentration'] = pd.concat([var2_TC_results_merged_df.iloc[i, 0:-2]
                                                                       for i in
                                                                       range(len(var2_TC_results_merged_df.index))],
                                                                      axis=0, ignore_index=True)
    # plot the mean and 95% confidence interval by aggregating over PE replicates (at each time point)
    sns.lineplot(data=var2_TC_results_list_merged_longform[var2_TC_results_list_merged_longform.Species.isin(['P', 'S', 'X'])],
                 x="Time", y="Concentration", hue="Species", ax=ax, linestyle='solid',
                 palette=[sns.color_palette(colors)[0], sns.color_palette(colors)[1], sns.color_palette(colors)[2]])
    # add experimental data points
    sns.scatterplot(data=exp_data_dataframes[exp_ID], x='Time', y='[P]', ax=ax,
                    color=sns.color_palette(colors)[0], s=scatterplot_marker_size)
    sns.scatterplot(data=exp_data_dataframes[exp_ID], x='Time', y='[S]', ax=ax,
                    color=sns.color_palette(colors)[1], s=scatterplot_marker_size)
    sns.scatterplot(data=exp_data_dataframes[exp_ID], x='Time', y='[X]', ax=ax,
                    color=sns.color_palette(colors)[2], s=scatterplot_marker_size)
    # overwrite axis labels which were automatically generated from data frame column names
    ax.set(xlabel='Time [h]', ylabel='Concentration [mM]')
    # create legend
    ax.legend(handles=[Line2D([0], [0], color='black', linewidth=1, linestyle='dashed'),
                       Line2D([0], [0], color='black', linewidth=1, linestyle='solid'),
                       Line2D([0], [0], color='black', marker='.', linestyle='None')],
              labels=['Initial core model', 'Selected model', 'Data'])
    # save the plot and the associated time course simulation data 
    plt.savefig(f'{plot_name}_comparison_merged.png', dpi=200)
    list_of_TC_results = [var1_TC_results_merged_df, var2_TC_results_merged_df]
    with open(f'{plot_name}_TC_results_merged.pkl', 'wb') as f:
        pickle.dump(list_of_TC_results, f)

    # ------------------------------------------------------------------------------------------------------ #

    # FILE MANAGEMENT

    # move all files in my current working directory (except the script file itself) into a new result directory
    os.makedirs(test_scenario, exist_ok=True)
    for filename in os.listdir():
        # skip directories and the script file
        if os.path.isfile(filename) and filename != "run_rediscovery_tests.py":
            # move the file to the subdirectory
            shutil.move(filename, os.path.join(test_scenario, filename))

# RUN TEST

# set standard deviation of Gaussian noise (that is added to the synthetic data generated from the true model)
noise_std = 0.2  # if set to 0 the noise is turned off
# define variable terms (those that will be successively added) where the first two elements are the row and column coordinates pointing to elements in the start variant data frame and the third element is the index of the base or regulation term; tuples are placed in inner lists because combinations of tuples are also valid elements of the vari_terms list
vari_terms = [[(1,0,1)],                                    # r2.on
              [(2,0,1)],                                    # r3.on
              [(3,0,1)],                                    # r4.on
              [(0,1,5)],                                    # r1.X_Act
              [(0,1,6)]]                                    # r1.X_Inhib
# run rediscovery test and visualize the results (results are stored in a subdirectory)
calc_and_vis_rediscovery_test(noise_std, vari_terms)
