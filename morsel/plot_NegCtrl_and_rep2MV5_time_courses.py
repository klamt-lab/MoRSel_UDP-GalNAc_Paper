#!/usr/bin/env python
# -*- coding: utf-8 -*-

from func_lib import *

""" Script to analyze and visualize the results of a super structure optimization in various ways. The 
structure of the logs of evaluated variants is as follows: 

-> The outer layer is a list and each entry is a dictionary (keys: ['variant', 'estimation_results', 
'fitness']) containing results from a different structural variant:

    --> The 'variant' key points to a pandas DataFrame representation of the kinetic matrix of the 
    structural variant.

    --> The 'fitness' key points to a dictionary (keys: ['min_AIC', 'min_AICc', 'min_BIC', 'min_CIC1', 
    'min_CIC2', 'min_CIC3']) of the overall fitness measures, i.e., the model selection criteria minima 
    across all replicates.

    --> The 'estimation_results' key points to a list of parameter estimation replicates where each entry 
    is a dictionary (keys: ['structural_variant', 'reaction_schemes', 'fit_setup', 'fit_algorithm', 
    'fit_algorithm_settings', 'estimated_parameters', 'fit_statistics', 'information_criteria'])
"""

# LOAD CALCULATION RESULT LOGS
# ----------------------------
# for improved extension search load output dictionary that contains evaluated_vars_log and the analysis_result of the final model selection step:
# evaluated_vars_log: log of evaluation results (each of which is a return object from eval_struct_var()) [list]
# analysis_result: a dictionary of the fitness data frame and the identified model variants with the best mean and median fitness [dictionary]

# load parameter ensemble of the new v23 initial core model (negative control; only core reactions activated and P_Inhib for all enzyme-catalyzed reactions)
reading_file = open('v23bNegControl_Exp36.37.39_50runs_evaluated_vars_log', 'rb')
neg_ctrl_log = pickle.load(reading_file)
reading_file.close()

# load parameter ensemble of the best model variant that was selected from an ImpExtSearch result
reading_file = open('rep2//ImpExtSearch_v23b_rep2_ModelVar5_50runs_evaluated_vars_log', 'rb')
selected_model_var = pickle.load(reading_file)
reading_file.close()

# VISUALIZATION
# -------------
# Visualize the performance of two selected structural variants
var1 = neg_ctrl_log[0] # v23b negative control model (P_Inhib for all reactions)
var2 = selected_model_var[0] # param. ensemble n=50 of best ranking model variant from ImpExtSearch
# load experimental data (measured by Tuan in the UDP-GalNAc experiments 36, 37 and 39)
exp_data_file_names = ['UDP-GalNAc36_5_with_initConcColumns.txt',
                       'UDP-GalNAc36_10_with_initConcColumns.txt',
                       'UDP-GalNAc36_20_with_initConcColumns.txt',
                       'UDP-GalNAc36_50_with_initConcColumns.txt',
                       'UDP-GalNAc37_20A_with_initConcColumns.txt',
                       'UDP-GalNAc37_20B_with_initConcColumns.txt',
                       'UDP-GalNAc37_50A_with_initConcColumns.txt',
                       'UDP-GalNAc37_50B_with_initConcColumns.txt',
                       'UDP-GalNAc39_LB_with_initConcColumns.txt',
                       'UDP-GalNAc39_LO_with_initConcColumns.txt',
                       'UDP-GalNAc39_TB_with_initConcColumns.txt',
                       'UDP-GalNAc39_TO_with_initConcColumns.txt']
exp_data_dataframes = load_and_process_exp_data(exp_data_file_names)
# select the experimental data which will be used for plotting and to initialize 
# the model simulations
plot_name = 'v23bNegCtrldashed_v23IESrep2MV5solid'
exp_names = ['UDP-GalNAc36_5',
             'UDP-GalNAc36_10',
             'UDP-GalNAc36_20',
             'UDP-GalNAc36_50',
             'UDP-GalNAc37_20A',
             'UDP-GalNAc37_20B',
             'UDP-GalNAc37_50A',
             'UDP-GalNAc37_50B',
             'UDP-GalNAc39_LB',
             'UDP-GalNAc39_LO',
             'UDP-GalNAc39_TB',
             'UDP-GalNAc39_TO']
# global plotting options
sns.set(font_scale = 1.25)
sns.set_style('whitegrid')
for exp_ID, exp_name in enumerate(exp_names):
    # call plotting function which will save the plot and the associated simulation data
    compare_struct_var_tcs_merged(var1, var2, exp_data_dataframes, exp_ID, exp_name, plot_name)
