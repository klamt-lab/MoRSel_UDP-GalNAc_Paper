#!/usr/bin/env python
# -*- coding: utf-8 -*-

from func_lib import *
import numpy as np
import pandas as pd
import pickle

""" Create a summary of the superstructure optimization results of v23b (MV5, rep1_MV5, rep2_MV5, rep3_MV6, rep4_MV3) and compare them to the base model that was used as starting point for the Improved Extension Search (NegCtrl). Also include additional test scenarios where only certain model selection critera were used for the Imp. Ext. Search (onlyAIC, onlyAICc, onlyBIC, onlyCIC1, onlyCIC2, onlyCIC3, AIC_AICc_BIC).
"""

# LOAD DATA AND PREPARE RESULTS
# load ImpExtSearch result logs
ImpExtSearch_logs_paths = {'rep0': 'SSO_ImpExtSearch_v23b_output_dict',
                           'rep1': 'rep1//SSO_ImpExtSearch_v23b_rep1_output_dict',
                           'rep2': 'rep2//SSO_ImpExtSearch_v23b_rep2_output_dict',
                           'rep3': 'rep3//SSO_ImpExtSearch_v23b_rep3_output_dict',
                           'rep4': 'rep4//SSO_ImpExtSearch_v23b_rep4_output_dict',
                           'onlyAIC': 'onlyAIC//SSO_ImpExtSearch_v23b_onlyAIC_output_dict',
                           'onlyAICc': 'onlyAICc//SSO_ImpExtSearch_v23b_onlyAICc_output_dict',
                           'onlyBIC': 'onlyBIC//SSO_ImpExtSearch_v23b_onlyBIC_output_dict',
                           'onlyCIC1': 'onlyCIC1//SSO_ImpExtSearch_v23b_onlyCIC1_output_dict',
                           'onlyCIC2': 'onlyCIC2//SSO_ImpExtSearch_v23b_onlyCIC2_output_dict',
                           'onlyCIC3': 'onlyCIC3//SSO_ImpExtSearch_v23b_onlyCIC3_output_dict',
                           'AIC_AICc_BIC': 'AIC_AICc_BIC//SSO_ImpExtSearch_v23b_AIC_AICc_BIC_output_dict'}
dict_of_ImpExtSearch_logs = {}
for name, path in ImpExtSearch_logs_paths.items():
    with open(path, 'rb') as f:
        dict_of_ImpExtSearch_logs[name] = pickle.load(f)
# load pickled evaluated vars. logs of the base model (NegCtrl) and of all selected best ranking model variants from each ImpExtSearch run
final_param_ensembles_paths = {'NegCtrl': 'v23bNegControl_Exp36.37.39_50runs_evaluated_vars_log',
                               'MV5': 'ImpExtSearch_v23b_ModelVar5_50runs_evaluated_vars_log',
                               'rep1_MV6': 'rep1//ImpExtSearch_v23b_rep1_ModelVar6_50runs_evaluated_vars_log',
                               'rep2_MV5': 'rep2//ImpExtSearch_v23b_rep2_ModelVar5_50runs_evaluated_vars_log',
                               'rep3_MV4': 'rep3//ImpExtSearch_v23b_rep3_ModelVar4_50runs_evaluated_vars_log',
                               'rep4_MV6': 'rep4//ImpExtSearch_v23b_rep4_ModelVar6_50runs_evaluated_vars_log',
                               'onlyAIC_MV17': 'onlyAIC//ImpExtSearch_v23b_onlyAIC_ModelVar17_50runs_evaluated_vars_log',
                               'onlyAICc_MV17': 'onlyAICc//ImpExtSearch_v23b_onlyAICc_ModelVar17_50runs_evaluated_vars_log',
                               'onlyBIC_MV16': 'onlyBIC//ImpExtSearch_v23b_onlyBIC_ModelVar16_50runs_evaluated_vars_log',
                               'onlyCIC1_MV2': 'onlyCIC1//ImpExtSearch_v23b_onlyCIC1_ModelVar2_50runs_evaluated_vars_log',
                               'onlyCIC2_MV2': 'onlyCIC2//ImpExtSearch_v23b_onlyCIC2_ModelVar2_50runs_evaluated_vars_log',
                               'onlyCIC3_MV2': 'onlyCIC3//ImpExtSearch_v23b_onlyCIC3_ModelVar2_50runs_evaluated_vars_log',
                               'AIC_AICc_BIC_MV11': 'AIC_AICc_BIC//ImpExtSearch_v23b_AIC_AICc_BIC_ModelVar11_50runs_evaluated_vars_log'}
dict_of_param_ensembles = {}
for name, path in final_param_ensembles_paths.items():
    with open(path, 'rb') as f:
        dict_of_param_ensembles[name] = pickle.load(f)
# save the best objective value that was reached in each parameter ensemble
param_ensembles_best_obj_vals = {}
for name, param_ensemble_info in dict_of_param_ensembles.items():
    best_obj_val = min([param_ensemble_info[0]['estimation_results'][i]['fit_statistics']['obj'] for i in range(len(param_ensemble_info[0]['estimation_results']))])
    param_ensembles_best_obj_vals[name] = best_obj_val
# for each ImpExtSearch result get the 'trace' from the starting point to the selected best ranking model variant (the 'trace' is the order of the terms that were added from the start to the selected model variant)
best_reps_and_MVs = [key.split('_') for key in dict_of_param_ensembles.keys() if 'MV' in key]
best_MVs = [item for item in sum(best_reps_and_MVs, []) if 'MV' in item]
best_MVs_indices = [int(x.split('MV')[1]) for x in best_MVs]
dict_of_traces = {'NegCtrl': 'no trace calculation necessary',
                  'rep0': None,
                  'rep1': None,
                  'rep2': None,
                  'rep3': None,
                  'rep4': None,
                  'onlyAIC': None,
                  'onlyAICc': None,
                  'onlyBIC': None,
                  'onlyCIC1': None,
                  'onlyCIC2': None,
                  'onlyCIC3': None,
                  'AIC_AICc_BIC': None}
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
             'num_estim_params': [len(entry[0]['estimation_results'][0]['estimated_parameters']) for entry in dict_of_param_ensembles.values()],
             'best_fit': np.round([obj_val for obj_val  in param_ensembles_best_obj_vals.values()], 5), 
             'AIC': np.round([entry[0]['fitness']['min_AIC'] for entry in dict_of_param_ensembles.values()], 2), 
             'AICc': np.round([entry[0]['fitness']['min_AICc'] for entry in dict_of_param_ensembles.values()], 2),
             'BIC': np.round([entry[0]['fitness']['min_BIC'] for entry in dict_of_param_ensembles.values()], 2),
             'CIC1': np.round([entry[0]['fitness']['min_CIC1'] for entry in dict_of_param_ensembles.values()], 2),
             'CIC2': np.round([entry[0]['fitness']['min_CIC2'] for entry in dict_of_param_ensembles.values()], 2),
             'CIC3': np.round([entry[0]['fitness']['min_CIC3'] for entry in dict_of_param_ensembles.values()], 2),
             'trace': [trace for run_ID, trace in dict_of_traces.items()]}
# create pandas data frame for report from the dictionary
report_df = pd.DataFrame(data_dict)
# create extra columns that show the average and median ranks across all 6 model selection criteria for the original 5 runs and for all runs
report_df_ranked = report_df.loc[0:6, ['AIC', 'AICc', 'BIC', 'CIC1', 'CIC2', 'CIC3']].rank()
report_df['original5_MSC_median_rank'] = np.round(report_df_ranked.median(axis=1), 2)
report_df['original5_MSC_mean_rank'] = np.round(report_df_ranked.mean(axis=1), 2)
report_df_ranked_all = report_df.loc[:, ['AIC', 'AICc', 'BIC', 'CIC1', 'CIC2', 'CIC3']].rank()
report_df['allRuns_MSC_median_rank'] = np.round(report_df_ranked_all.median(axis=1), 2)
report_df['allRuns_MSC_mean_rank'] = np.round(report_df_ranked_all.mean(axis=1), 2)
# reorder columns
report_df = report_df[['Model', 'num_estim_params', 'best_fit', 'AIC', 'AICc', 'BIC', 'CIC1', 'CIC2', 'CIC3', 'original5_MSC_mean_rank', 'original5_MSC_median_rank', 'allRuns_MSC_mean_rank', 'allRuns_MSC_median_rank', 'trace']]
# create an alternative version of the report DataFrame where traces are expressed not as lists of term names but instead as sequences of numbers (where each number refers to a term name)
trans_dict = {('ADP_Decay_v1', 'ADP_Decay.Active.Variant1'): '1',
              ('ADP_Decay_v2', 'ADP_Decay.Active.Variant2'): '2',
              ('GLMU', 'ADP_Inhib'):                         '3',
              ('GLMU', 'UDP_GalNAc_Inhib'):                  '4',
              ('NAHK_ATP', 'ADP_Variant'):                   '5',
              ('NAHK_ATP', 'ATPP_Variant'):                  '6',
              ('NAHK_ATP', 'ADP_and_ATPP_Variant'):          '7',
              ('NAHK_ATP', 'ADP_Inhib'):                     '8',
              ('PPA', 'ADP_Inhib'):                          '9',
              ('PPA', 'ATP_Act'):                            '10',
              ('PPA', 'PP_Act'):                             '11',
              ('PPK3_A', 'ATPP_Variant'):                    '12',
              ('PPK3_A', 'ADP_Inhib'):                       '13',
              ('PPK3_A', 'PP_Inhib'):                        '14',
              ('UDK_ATP', 'ADP_Variant'):                    '15',
              ('UDK_ATP', 'ATPP_Variant'):                   '16',
              ('UDK_ATP', 'ADP_and_ATPP_Variant'):           '17',
              ('UDK_ATP', 'ADP_Inhib'):                      '18',
              ('UDK_ATP', 'UTP_Inhib'):                      '19',
              ('UMPK_ATP', 'ADP_Variant'):                   '20',
              ('UMPK_ATP', 'ATPP_Variant'):                  '21',
              ('UMPK_ATP', 'ADP_and_ATPP_Variant'):          '22',
              ('UMPK_ATP', 'ADP_Inhib'):                     '23'}
report_df_alt = report_df.copy(deep=True)
for row_idx, trace_original in report_df_alt.loc[1::, 'trace'].items():
    trace_translated = []
    for term in trace_original:
        # the first value of the term tuple is the number denoting when it was added durin the Imp. Ext. Search so we can skip that for the look up in the translation dictionary
        trace_translated.append(trans_dict[term[1::]])
    report_df_alt.at[row_idx, 'trace'] = trace_translated

# EXPORT REPORTS
# save report data frames as .csv files
report_df.to_csv('IESv23b_SSO_report.csv', index=False)
report_df_alt.to_csv('IESv23b_SSO_report_numerical_traces.csv', index=False)
