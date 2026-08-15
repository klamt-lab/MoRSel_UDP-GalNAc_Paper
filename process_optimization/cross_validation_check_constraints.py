#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import os
import pandas as pd
import pickle
from pprint import pformat

working_dir = os.getcwd()

def cross_validation_check_constraints(model_variant_name, opt_variant_name, path_to_files):
    """ Check how often a relevant constraint is violated in the corresponding cross-validation table.
    
    :param model_variant_name: name of the model variant that was used to calculate the optimizations
    :type model_variant_name: string
    :param opt_variant_name: optimization variant ('a' or 'b')
    :type opt_variant_name: string
    :param path_to_files: path to the result files of the optimization (this is also where the generated figures will be stored)
    :type path_to_files:
    """

    # set working directory to the provided path so that result files are found and generated figures are saved to the correct directory
    os.chdir(path_to_files)

    # load relevant cross-validation tables:
    # some optimization problems include constraints for model variables that are also contained in a cross-validation table (for example the optimization problem could contain a minimum titer constraint and one of the cross-validation tables lists the titers for all combinations of optimized initial concentrations O1 ... On and parameter sets P1 ... Pn)
    # -> currently, this is only relevant for OP05 (while OP01 includes an E_tot constraint there is no E_tot cross-validation table)
    with open(f'{model_variant_name}_PartSwarm50x_OP05withS_CrossVal_allEtot.pkl', 'rb') as file:
        cross_val_result_OP05withS_allEtot = pickle.load(file)
    with open(f'{model_variant_name}_PartSwarm50x_OP05withS_CrossVal_allTit.pkl', 'rb') as file:
        cross_val_result_OP05withS_allTit = pickle.load(file)
    with open(f'{model_variant_name}_PartSwarm50x_OP05withS_CrossVal_allYields.pkl', 'rb') as file:
        cross_val_result_OP05withS_allYields = pickle.load(file)

    # collect pandas DataFrames of all tables in a dictionary
    df_col_names = [f'O{i}' for i in range(1, 51)]  # optimization results O1, ..., O50
    df_row_names = [f'P{i}' for i in range(1, 51)]  # parameter sets P1, ..., P50
    cross_val_tables_dict = {'OP05withS': {'titers': pd.DataFrame(cross_val_result_OP05withS_allTit, columns=df_col_names, index=df_row_names), 
                                           'yields': pd.DataFrame(cross_val_result_OP05withS_allYields, columns=df_col_names, index=df_row_names), 
                                           'E_tot': pd.DataFrame(cross_val_result_OP05withS_allEtot, columns=df_col_names, index=df_row_names)}}

    # create a dictionary with the relevant constraint values
    constraints_dict = {'OP05withS': {'titer_lower': 40.37,     # mM (t = t_end)
                                      'yield_lower': 0.8274}}   # [-] (t = t_end)

    # go through the relevant cross-validation tables and check how often the associated constraint was violated
    cross_val_tables_constraint_eval_dict = {'OP05withS': {'titers': {'total_table_size': cross_val_tables_dict['OP05withS']['titers'].size,
                                                                      'num_constraint_violated': int(((cross_val_tables_dict['OP05withS']['titers'] < constraints_dict['OP05withS']['titer_lower']).sum().sum())),
                                                                      'percentage_constraint_violated': float((((cross_val_tables_dict['OP05withS']['titers'] < constraints_dict['OP05withS']['titer_lower']).sum().sum()) / cross_val_tables_dict['OP05withS']['titers'].size) * 100)},
                                                           'yields': {'total_table_size': cross_val_tables_dict['OP05withS']['yields'].size,
                                                                      'num_constraint_violated': int(((cross_val_tables_dict['OP05withS']['yields'] < constraints_dict['OP05withS']['yield_lower']).sum().sum())),
                                                                      'percentage_constraint_violated': float((((cross_val_tables_dict['OP05withS']['yields'] < constraints_dict['OP05withS']['yield_lower']).sum().sum()) / cross_val_tables_dict['OP05withS']['yields'].size) * 100)}}}

    return cross_val_tables_constraint_eval_dict

# optimization variant b: PolyP fixed at 32 mM and ATP fixed at 0.5 mM
model_dir_info = [('ImpExtSearch_v23b_rep2_MV5', 'b', f'{working_dir}\\b_constPolyP32mMATP0.5mM\\rep2_MV5')]
opt_b_cross_val_tables_constraint_eval_results = {}
for (name, opt_var, path) in model_dir_info:
    opt_b_cross_val_tables_constraint_eval_results.update({f'{name}_{opt_var}': cross_validation_check_constraints(name, opt_var, path)})
# export results as pickle and text files
with open(os.path.join(working_dir, f'opt_b_cross_val_tables_constraint_eval_results.pkl'), 'wb') as file:
    pickle.dump(opt_b_cross_val_tables_constraint_eval_results, file)
with open(os.path.join(working_dir, f'opt_b_cross_val_tables_constraint_eval_results.txt'), 'w') as file:
    file.write(pformat(opt_b_cross_val_tables_constraint_eval_results, indent=4))
