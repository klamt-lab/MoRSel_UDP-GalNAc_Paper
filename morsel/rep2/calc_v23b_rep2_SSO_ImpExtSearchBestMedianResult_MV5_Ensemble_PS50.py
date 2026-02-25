#!/usr/bin/env python
# -*- coding: utf-8 -*-

import time
import sys

sys.path.append("..")
from func_lib import *

# load an ImpExtSearch result
reading_file = open('SSO_ImpExtSearch_v23b_rep2_output_dict', 'rb')
output_dict = pickle.load(reading_file)
reading_file.close()

# analyze single output log
analysis_result = output_dict['analysis_result']
# get the index of the model variant with the best median rank and the fewest additions
best_median_rank_model_idx = analysis_result["fitness_dataframe"]['median_rank'].idxmin()
# get the index of the model variant with the best mean rank and the fewest additions
best_mean_rank_model_idx = analysis_result["fitness_dataframe"]['mean_rank'].idxmin()
# compare both and select the model with the lowest index (= lowest number of additions)
if best_median_rank_model_idx < best_mean_rank_model_idx:
    best_ranking_var_log_dict = output_dict['evaluated_vars_log'][best_median_rank_model_idx]
    print(f"Model variant {best_median_rank_model_idx} with best median rank {analysis_result['fitness_dataframe'].iloc[best_median_rank_model_idx, :]['median_rank']} was selected.")
elif best_median_rank_model_idx > best_mean_rank_model_idx:
    best_ranking_var_log_dict = output_dict['evaluated_vars_log'][best_mean_rank_model_idx]
    print(f"Model variant {best_mean_rank_model_idx} with best mean rank {analysis_result['fitness_dataframe'].iloc[best_mean_rank_model_idx, :]['mean_rank']} was selected.")
elif best_median_rank_model_idx == best_mean_rank_model_idx:
    # both indices are the same so they both point to the same model variant -> therefore it doesn't matter which one is used to look it up in the ImpExtSearch result
    best_ranking_var_log_dict = output_dict['evaluated_vars_log'][best_median_rank_model_idx]
    print(f"Model variant {best_median_rank_model_idx} with best median rank {analysis_result['fitness_dataframe'].iloc[best_median_rank_model_idx, :]['median_rank']} and best mean rank {analysis_result['fitness_dataframe'].iloc[best_mean_rank_model_idx, :]['mean_rank']} was selected.")
# output: "Model variant 5 with best median rank 6.0 and best mean rank 9.333333333333334 was selected."

# create parameter ensemble of the selected model variant 
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
exp_data_dataframes = load_and_process_exp_data([f'../{name}' for name in exp_data_file_names])  # look for data in parent directory
model_variant_name = 'ImpExtSearch_v23b_rep2_ModelVar5'
evaluated_vars_log = create_parameter_ensemble(best_ranking_var_log_dict['variant'], 50, 
                                               exp_data_file_names, exp_data_dataframes, model_variant_name)
