#!/usr/bin/env python
# -*- coding: utf-8 -*-

from func_lib import *
import numpy as np
import pandas as pd
import pickle

""" Create a summary of the superstructure optimization results of v23b (MV5, rep1_MV5, rep2_MV5, rep3_MV6, rep4_MV3) and compare them to the base model that was used as starting point for the Improved Extension Search (NegCtrl).
"""


# LOAD DATA AND PREPARE RESULTS
# load ImpExtSearch result logs
reading_file = open('SSO_ImpExtSearch_v23b_output_dict', 'rb')
output_dict_0 = pickle.load(reading_file)
reading_file.close()
reading_file = open('rep1//SSO_ImpExtSearch_v23b_rep1_output_dict', 'rb')
output_dict_rep1 = pickle.load(reading_file)
reading_file.close()
reading_file = open('rep2//SSO_ImpExtSearch_v23b_rep2_output_dict', 'rb')
output_dict_rep2 = pickle.load(reading_file)
reading_file.close()
reading_file = open('rep3//SSO_ImpExtSearch_v23b_rep3_output_dict', 'rb')
output_dict_rep3 = pickle.load(reading_file)
reading_file.close()
reading_file = open('rep4//SSO_ImpExtSearch_v23b_rep4_output_dict', 'rb')
output_dict_rep4 = pickle.load(reading_file)
reading_file.close()
dict_of_ImpExtSearch_logs = {'rep0': output_dict_0, 
                             'rep1': output_dict_rep1, 
                             'rep2': output_dict_rep2, 
                             'rep3': output_dict_rep3, 
                             'rep4': output_dict_rep4}
# load parameter ensembles of the base model (NegCtrl) and of all selected best ranking model variants from each ImpExtSearch run
param_sets_NegCtrl = pd.read_csv('sampling_output_Particle_Swarm50runs_v23bNegControl_Exp36.37.39.csv')
param_sets_rep0 = pd.read_csv('sampling_output_Particle_Swarm50runs_ImpExtSearch_v23b_ModelVar5.csv')
param_sets_rep1 = pd.read_csv('rep1//sampling_output_Particle_Swarm50runs_ImpExtSearch_v23b_rep1_ModelVar6.csv')
param_sets_rep2 = pd.read_csv('rep2//sampling_output_Particle_Swarm50runs_ImpExtSearch_v23b_rep2_ModelVar5.csv')
param_sets_rep3 = pd.read_csv('rep3//sampling_output_Particle_Swarm50runs_ImpExtSearch_v23b_rep3_ModelVar4.csv')
param_sets_rep4 = pd.read_csv('rep4//sampling_output_Particle_Swarm50runs_ImpExtSearch_v23b_rep4_ModelVar6.csv')
dict_of_param_ensembles = {'NegCtrl': param_sets_NegCtrl, 
                           'MV5': param_sets_rep0, 
                           'rep1_MV6': param_sets_rep1, 
                           'rep2_MV5': param_sets_rep2, 
                           'rep3_MV4': param_sets_rep3, 
                           'rep4_MV6': param_sets_rep4}
# for each ImpExtSearch result get the 'trace' from the starting point to the selected best ranking model variant (the 'trace' is the order of the terms that were added from the start to the selected model variant)
best_reps_and_MVs = [key.split('_') for key in dict_of_param_ensembles.keys() if 'MV' in key]
best_MVs = [item for item in sum(best_reps_and_MVs, []) if 'MV' in item]
best_MVs_indices = [int(x.split('MV')[1]) for x in best_MVs]
dict_of_traces = {'NegCtrl': 'no trace calculation necessary',
                  'rep0': None,
                  'rep1': None,
                  'rep2': None,
                  'rep3': None,
                  'rep4': None}
for (run_ID, output_dict), best_MV_idx in zip(dict_of_ImpExtSearch_logs.items(), best_MVs_indices):
    # set starting point of the trace calculation (= unchanged start variant of the ImpExtSearch result)
    start_var = output_dict['evaluated_vars_log'][0]['variant']
    # set end point of the trace calculation (= best ranking model variant of the ImpExtSearch result)
    selected_var = output_dict['evaluated_vars_log'][best_MV_idx]['variant']
    # get trace
    selected_var_trace = get_imp_ext_search_trace_for_selected_var(start_var, 
                                                                   selected_var,
                                                                   output_dict['evaluated_vars_log'])
    dict_of_traces.update({run_ID: selected_var_trace})

# CREATE REPORT
# collect related data of all optimization problems in separate new dictionaries
# each value of the data dictionary is a list of 6 elements (for each of the 6 model variants)
data_dict = {'Model': [run_ID for run_ID in dict_of_param_ensembles.keys()], 
             'num_estim_params': [len(param_ensemble.columns)-2 for param_ensemble in dict_of_param_ensembles.values()], # subtract 2 because of index and obj.val columns
             'best_fit': [np.round(param_ensemble.min()['obj'], 5) for param_ensemble in dict_of_param_ensembles.values()], 
             'AIC': [np.round(1080*np.log(param_ensemble.min()['obj']/1080)+2*((len(param_ensemble.columns)-2)+1), 2) for param_ensemble in dict_of_param_ensembles.values()], # AIC = n*log(RSS/n)+2*k; k = num_estim_params+1; 1080 exp. data points
             'AICc': [np.round(1080*np.log(param_ensemble.min()['obj']/1080)+((2*((len(param_ensemble.columns)-2)+1)*(((len(param_ensemble.columns)-2)+1)+1))/(1080-((len(param_ensemble.columns)-2)+1)-1)), 2) for param_ensemble in dict_of_param_ensembles.values()], # AICc = n*log(RSS/n)+((2*k*(k+1))/(n-k-1)); k = num_estim_params+1; 1080 exp. data points
             'BIC': [np.round(1080*np.log(param_ensemble.min()['obj']/1080)+np.log(1080)*((len(param_ensemble.columns)-2)+1), 2) for param_ensemble in dict_of_param_ensembles.values()], # BIC = n*log(RSS/n)+log(n)*k; k = num_estim_params+1; 1080 exp. data points
             'CIC1': [np.round(1080*np.log(param_ensemble.min()['obj']/1080)+(np.abs((1080*np.log(param_ensemble.min()['obj']/1080)))*0.01)*((len(param_ensemble.columns)-2)+1), 2) for param_ensemble in dict_of_param_ensembles.values()], # CIC1 = n*log(RSS/n)+(abs((n*log(RSS/n)))*0.01)*k; k = num_estim_params+1; 1080 exp. data points
             'CIC2': [np.round(1080*np.log(param_ensemble.min()['obj']/1080)+200*((len(param_ensemble.columns)-2)+1), 2) for param_ensemble in dict_of_param_ensembles.values()], # CIC2 = n*log(RSS/n)+200*k; k = num_estim_params+1; 1080 exp. data points
             'CIC3': [np.round(1080*np.log(param_ensemble.min()['obj']/1080)*1/((len(param_ensemble.columns)-2)+1), 2) for param_ensemble in dict_of_param_ensembles.values()], # CIC3 = (n*log(RSS/n))*1/k; k = num_estim_params+1; 1080 exp. data points
             'trace': [trace for run_ID, trace in dict_of_traces.items()]}
# create pandas data frame for report from the dictionary
report_df = pd.DataFrame(data_dict)
# create extra columns that show the average and median ranks across all 6 model selection criteria
report_df_ranked = report_df.loc[:, ['AIC', 'AICc', 'BIC', 'CIC1', 'CIC2', 'CIC3']].rank()
report_df['MSC_median_rank'] = np.round(report_df_ranked.median(axis=1), 2)
report_df['MSC_mean_rank'] = np.round(report_df_ranked.mean(axis=1), 2)
# reorder columns
report_df = report_df[['Model', 'num_estim_params', 'best_fit', 'AIC', 'AICc', 'BIC', 'CIC1', 'CIC2', 'CIC3', 'MSC_mean_rank', 'MSC_median_rank', 'trace']]

# EXPORT REPORT
# save report data frames as .csv files
report_df.to_csv('IESv23b_SSO_report.csv', index=False)
