#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

sys.path.append("..")
from batch_opt_func_lib import *

""" Select the best optimization result (best prediction) by calculating the objective value for all combinations of all parameter sets p with all optimization results O.

    Approach: Load the list of parameter sets from the provided ensemble and load the list of optimization results (= sets of optimal initial enzyme and substrate concentrations) from 'opt_calc_and_vis.py'. Both lists are of equal size (usually 100). Then, build a two level loop to iterate over all possible combinations of i parameter sets p and j optimization results O. Set the respective parameters and initial enzyme and substrate concentrations in the model and run a simulation. Store the results in a i x j matrix (rows: parameter sets, columns: optimization results). Then, calculate different measures (median, minimum, sum) for each column (= each optimization result) and combine them to calculate a score for each column. The optimization result with the highest score is selected for experimental testing.

    Main package: basiCO - simplified Copasi Python API
                  <https://github.com/copasi/basico>"""

# SETUP
model_path = 'ImpExtSearch_v23b_rep2_MV5.cps'
fits_data_path = 'sampling_output_Particle_Swarm50runs_ImpExtSearch_v23b_rep2_ModelVar5.csv'

# OP05withS
""" Scoring function: '(col_enyzme_load_value-0.53)/(overall_enzyme_load_min-0.53) + min(1, (titers_col_median/40.37)) + min(1, (yields_col_median/0.8274))'
"""
opt_results_path = 'ImpExtSearch_v23b_rep2_MV5_PartSwarm50x_OP05withS.pkl'
result_file_name = "ImpExtSearch_v23b_rep2_MV5_PartSwarm50x_OP05withS_CrossVal"
cross_val_table_E_tot, cross_val_table_titers, cross_val_table_yields, best_opt_result = cross_val_and_score_OP05withS(model_path, fits_data_path, opt_results_path, result_file_name)
