#!/usr/bin/env python
# -*- coding: utf-8 -*-

from basico import *
import math
import numpy as np
import pandas as pd
import pickle
import re
import sys
from tqdm import tqdm

""" Collection of functions used by the batch optimization calculation, analysis and visualization scripts as well as within other functions of this library.
"""

## here, the enzyme loads (both in constraints and in target functions) are defined using grams per litre!

# ------------------------------------------------------------------------------------------------------ #

# DEFINITION AND CALCULATION OF OPTIMIZATION PROBLEMS

# 1) OP01_7h: maximize titer at 7h (s.t. maximum enzyme load; opt. vars: initial enzyme concentrations)
def solve_OP01_7h(model_path, fits_data_path, result_file_name):
    """ Solve a specific implementation of the optimization problem OP01 (maximize titer at 7h for a maximum enzyme load while reaching a certain minimum yield at 7h) with initial enzyme concentrations set as optimization variables. Repeat calculation for all parameter sets of a given parameter ensemble.

    :param model_path: path to the Copasi model object (.cps file)
    :type model_path: string
    :param fits_data_path: path to the parameter ensemble (.csv file)
    :type fits_data_path: string
    :param result_file_name: path of the calculated optimization result (.pkl file)
    :type fits_data_path: string

    :return opt_results: list of optimization results (elements: lists of result pandas data frames and dictionaries with information of the optimization calculation)
    :type opt_results: list
    """

    # LOAD MODEL OBJECT AND PARAMETER ENSEMBLE
    # the Copasi model file contains the structure of the model (stoichiometry, species, reactions, rate laws, ODEs); initial concentrations of all species, kinetic parameter values, and the specific optimization problem are set in this script
    model = load_model(model_path)
    # read the data of the parameter estimation sampling (the first column is all zeros, parameter columns start at .iloc[:,1])
    fits_data = pd.read_csv(fits_data_path)

    # SETUP OPTIMIZATION PROBLEM
    # set all species concentrations; ADP is naturally present in ATP stock solution at a rough ratio of 95% ATP/5% ADP but for the model simulation we assume an ideal 100% ATP stock solution (i.e., ADP is not present at t=0): the initial concentrations are from Tuan:Exp36_50mM (planned concentrations)
    set_species('ADP',          initial_concentration=1e-12)
    set_species('AMP',          initial_concentration=1e-12)
    set_species('ATP',          initial_concentration=0.5)
    set_species('ATPP',         initial_concentration=1e-12)
    set_species('E_GLMU',       initial_concentration=0.002032934)
    set_species('E_NAHK',       initial_concentration=0.002506077)
    set_species('E_PPA',        initial_concentration=0.001553358)
    set_species('E_PPK3',       initial_concentration=0.002878526)
    set_species('E_UDK',        initial_concentration=0.00410627)
    set_species('E_UMPK',       initial_concentration=0.004463289)
    set_species('GalNAc',       initial_concentration=50)
    set_species('GalNAc1P',     initial_concentration=1e-12)
    set_species('P',            initial_concentration=1e-12)
    set_species('PolyP',        initial_concentration=256.0)
    set_species('PP',           initial_concentration=1e-12)
    set_species('UDP',          initial_concentration=1e-12)
    set_species('UDP_GalNAc',   initial_concentration=1e-12)
    set_species('UMP',          initial_concentration=1e-12)
    set_species('UTP',          initial_concentration=1e-12)
    set_species('Uri',          initial_concentration=50)
    # create the necessary global quantities and events used for the target function (there can be no spaces in Copasi object names)
    add_parameter('UDP_GalNAc_at_7h', type='fixed', initial_value=0)
    add_event('save_UDP-GalNAc_at_7h', 'Time == 7', [('Values[UDP_GalNAc_at_7h]', '[UDP_GalNAc]')])
    # set target function
    set_objective_function(expression='Values[UDP_GalNAc_at_7h]', maximize=True)
    # setup optimization variables
    opt_vars = [{'name': '[E_GLMU]_0',
                'lower': 0.000203,
                'upper': 0.041,
                'start': 0.002032934},
                {'name': '[E_NAHK]_0',
                'lower': 0.000251,
                'upper': 0.05,
                'start': 0.002506077},
                {'name': '[E_PPA]_0',
                'lower': 0.000518,
                'upper': 0.103,
                'start': 0.001553358},
                {'name': '[E_PPK3]_0',
                'lower': 0.000288,
                'upper': 0.029,
                'start': 0.002878526},
                {'name': '[E_UDK]_0',
                'lower': 0.000411,
                'upper': 0.08,
                'start': 0.00410627},
                {'name': '[E_UMPK]_0',
                'lower': 0.000169,
                'upper': 0.034,
                'start': 0.004463289}]
    set_opt_parameters(opt_vars)
    # create the necessary global quantities used for the constraints (there can be no spaces in Copasi object names); enzyme concentrations in mM (mmol/l) are multiplied with their molecular weights (in g/mol divided by 1000 to get g/mmol) so that the resulting enzyme load is expressed in g/l
    add_parameter('E_tot_MW', type='assignment', expression='(49190/1000) * [E_GLMU]_0 + (39903/1000) * [E_NAHK]_0 + (19313/1000) * [E_PPA]_0 + (34740/1000) * [E_PPK3]_0 + (24353/1000) * [E_UDK]_0 + (22405/1000) * [E_UMPK]_0')
    # set constraints
    set_opt_constraints([
        {'name': 'Values[E_tot_MW]', 'lower': 0, 'upper': 0.53}
    ])
    # change time course settings (the optimization task calls the time course task)
    tc_settings = get_task_settings(T.TIME_COURSE)
    tc_settings['problem']['Duration'] = 7
    tc_settings['problem']['StepNumber'] = 100
    tc_settings['problem']['StepSize'] = 7/100
    set_task_settings(T.TIME_COURSE, settings=tc_settings)
    # set subtask and solver
    set_opt_settings(settings={'subtask': T.TIME_COURSE,
                               'method': {'name': PE.GENETIC_ALGORITHM}
                               })

    # SOLVE OPTIMIZATION PROBLEM FOR PARAMETER ENSEMBLE
    opt_results = list()
    # loop over result data frame (every row = one parameter set)
    for index, param_set in tqdm(fits_data.iterrows()):
        # current list of parameter values
        curr_param_set_vals = list(param_set)
        # current list of parameter names; list() to only get content (values or names) of pandas series; param_set.keys() to get names of values in pandas series param_set
        curr_param_set_names = list(param_set.keys())
        # trim both lists (get rid of first entries = name and last entries = obj_val of fit)
        curr_param_set_vals = curr_param_set_vals[1:-1]
        curr_param_set_names = curr_param_set_names[1:-1]
        # set model parameter values to current estimated parameters from rand sampling result
        for i in range(len(list(curr_param_set_vals))): 
            set_reaction_parameters(name = curr_param_set_names[i], value = curr_param_set_vals[i])
        # optimize the model with these updated parameter values
        opt_result = run_optimization() # opt_result is a pandas.core.frame.DataFrame
        opt_stats = get_opt_statistic() # opt_stats is a dictionary
        # get list of optimization variable start values
        opt_vars_start_vals = [opt_vars[0]['start'], 
                               opt_vars[1]['start'], 
                               opt_vars[2]['start'], 
                               opt_vars[3]['start'], 
                               opt_vars[4]['start'], 
                               opt_vars[5]['start']]
        # add start values of optimization variables to result data frame as column 3
        opt_result.insert(2, "start", opt_vars_start_vals, True)
        # store the optimization result data frame and the opt_stats dictionary as a list
        opt_results.append([opt_result, opt_stats])
    # export opt_results as python object
    open_file = open(result_file_name, "wb")
    pickle.dump(opt_results, open_file)
    open_file.close()
    # result structure:
    # outer list opt_results [ inner lists  [ df opt result, dict opt_stats  ]  [...]  [...] ...]
    # indexing: opt_results[0][0]           -> first opt result data frame
    #           opt_results[0][0].iloc[:,0] -> first column in first opt_result data frame
    #           opt_results[0][1]           -> first opt result statistics dictionary
    #           opt_results[0][1]['key']    -> entry of 'key' in first opt result  
    #                                          statistics dictionary
    return opt_results

# 2) OP01withS_7h: maximize titer at 7h (s.t. maximum enzyme load; opt. vars: initial enzyme and substrate concentrations [main substrates capped to baseline concentrations])
def solve_OP01withS_7h(model_path, fits_data_path, result_file_name):
    """ Solve a specific implementation of the optimization problem OP01 (maximize titer at 7h for a maximum enzyme load while reaching a certain minimum yield at 7h) with initial enzyme and substrate concentrations set as optimization variables. Repeat calculation for all parameter sets of a given parameter ensemble.

    :param model_path: path to the Copasi model object (.cps file)
    :type model_path: string
    :param fits_data_path: path to the parameter ensemble (.csv file)
    :type fits_data_path: string
    :param result_file_name: path of the calculated optimization result (.pkl file)
    :type fits_data_path: string

    :return opt_results: list of optimization results (elements: lists of result pandas data frames and dictionaries with information of the optimization calculation)
    :type opt_results: list
    """

    # LOAD MODEL OBJECT AND PARAMETER ENSEMBLE
    # the Copasi model file contains the structure of the model (stoichiometry, species, reactions, rate laws, ODEs); initial concentrations of all species, kinetic parameter values, and the specific optimization problem are set in this script
    model = load_model(model_path)
    # read the data of the parameter estimation sampling (the first column is all zeros, parameter columns start at .iloc[:,1])
    fits_data = pd.read_csv(fits_data_path)

    # SETUP OPTIMIZATION PROBLEM
    # set all species concentrations; ADP is naturally present in ATP stock solution at a rough ratio of 95% ATP/5% ADP but for the model simulation we assume an ideal 100% ATP stock solution (i.e., ADP is not present at t=0): the initial concentrations are from Tuan:Exp36_50mM (planned concentrations)
    set_species('ADP',          initial_concentration=1e-12)
    set_species('AMP',          initial_concentration=1e-12)
    set_species('ATP',          initial_concentration=0.5)
    set_species('ATPP',         initial_concentration=1e-12)
    set_species('E_GLMU',       initial_concentration=0.002032934)
    set_species('E_NAHK',       initial_concentration=0.002506077)
    set_species('E_PPA',        initial_concentration=0.001553358)
    set_species('E_PPK3',       initial_concentration=0.002878526)
    set_species('E_UDK',        initial_concentration=0.00410627)
    set_species('E_UMPK',       initial_concentration=0.004463289)
    set_species('GalNAc',       initial_concentration=50)
    set_species('GalNAc1P',     initial_concentration=1e-12)
    set_species('P',            initial_concentration=1e-12)
    set_species('PolyP',        initial_concentration=256.0)
    set_species('PP',           initial_concentration=1e-12)
    set_species('UDP',          initial_concentration=1e-12)
    set_species('UDP_GalNAc',   initial_concentration=1e-12)
    set_species('UMP',          initial_concentration=1e-12)
    set_species('UTP',          initial_concentration=1e-12)
    set_species('Uri',          initial_concentration=50)
    # create the necessary global quantities and events used for the target function (there can be no spaces in Copasi object names)
    add_parameter('UDP_GalNAc_at_7h', type='fixed', initial_value=0)
    add_event('save_UDP-GalNAc_at_7h', 'Time == 7', [('Values[UDP_GalNAc_at_7h]', '[UDP_GalNAc]')])
    # set target function
    set_objective_function(expression='Values[UDP_GalNAc_at_7h]', maximize=True)
    # setup optimization variables
    opt_vars = [{'name': '[E_GLMU]_0',
                'lower': 0.000203,
                'upper': 0.041,
                'start': 0.002032934},
                {'name': '[E_NAHK]_0',
                'lower': 0.000251,
                'upper': 0.05,
                'start': 0.002506077},
                {'name': '[E_PPA]_0',
                'lower': 0.000518,
                'upper': 0.103,
                'start': 0.001553358},
                {'name': '[E_PPK3]_0',
                'lower': 0.000288,
                'upper': 0.029,
                'start': 0.002878526},
                {'name': '[E_UDK]_0',
                'lower': 0.000411,
                'upper': 0.08,
                'start': 0.00410627},
                {'name': '[E_UMPK]_0',
                'lower': 0.000169,
                'upper': 0.034,
                'start': 0.004463289},
                {'name': '[GalNAc]_0',
                'lower': 0.1,
                'upper': 50,
                'start': 50},
                {'name': '[Uri]_0',
                'lower': 0.1,
                'upper': 50,
                'start': 50}]
    set_opt_parameters(opt_vars)
    # create the necessary global quantities used for the constraints (there can be no spaces in Copasi object names); enzyme concentrations in mM (mmol/l) are multiplied with their molecular weights (in g/mol divided by 1000 to get g/mmol) so that the resulting enzyme load is expressed in g/l
    add_parameter('E_tot_MW', type='assignment', expression='(49190/1000) * [E_GLMU]_0 + (39903/1000) * [E_NAHK]_0 + (19313/1000) * [E_PPA]_0 + (34740/1000) * [E_PPK3]_0 + (24353/1000) * [E_UDK]_0 + (22405/1000) * [E_UMPK]_0')
    # set constraints
    set_opt_constraints([
        {'name': 'Values[E_tot_MW]', 'lower': 0, 'upper': 0.53}
    ])
    # change time course settings (the optimization task calls the time course task)
    tc_settings = get_task_settings(T.TIME_COURSE)
    tc_settings['problem']['Duration'] = 7
    tc_settings['problem']['StepNumber'] = 100
    tc_settings['problem']['StepSize'] = 7/100
    set_task_settings(T.TIME_COURSE, settings=tc_settings)
    # set subtask and solver
    set_opt_settings(settings={'subtask': T.TIME_COURSE,
                               'method': {'name': PE.GENETIC_ALGORITHM}
                               })

    # SOLVE OPTIMIZATION PROBLEM FOR PARAMETER ENSEMBLE
    opt_results = list()
    # loop over result data frame (every row = one parameter set)
    for index, param_set in tqdm(fits_data.iterrows()):
        # current list of parameter values
        curr_param_set_vals = list(param_set)
        # current list of parameter names; list() to only get content (values or names) of pandas series; param_set.keys() to get names of values in pandas series param_set
        curr_param_set_names = list(param_set.keys())
        # trim both lists (get rid of first entries = name and last entries = obj_val of fit)
        curr_param_set_vals = curr_param_set_vals[1:-1]
        curr_param_set_names = curr_param_set_names[1:-1]
        # set model parameter values to current estimated parameters from rand sampling result
        for i in range(len(list(curr_param_set_vals))): 
            set_reaction_parameters(name = curr_param_set_names[i], value = curr_param_set_vals[i])
        # optimize the model with these updated parameter values
        opt_result = run_optimization() # opt_result is a pandas.core.frame.DataFrame
        opt_stats = get_opt_statistic() # opt_stats is a dictionary
        # get list of optimization variable start values
        opt_vars_start_vals = [opt_vars[0]['start'], 
                               opt_vars[1]['start'], 
                               opt_vars[2]['start'], 
                               opt_vars[3]['start'], 
                               opt_vars[4]['start'], 
                               opt_vars[5]['start'],
                               opt_vars[6]['start'],
                               opt_vars[7]['start']]
        # add start values of optimization variables to result data frame as column 3
        opt_result.insert(2, "start", opt_vars_start_vals, True)
        # store the optimization result data frame and the opt_stats dictionary as a list
        opt_results.append([opt_result, opt_stats])
    # export opt_results as python object
    open_file = open(result_file_name, "wb")
    pickle.dump(opt_results, open_file)
    open_file.close()
    # result structure:
    # outer list opt_results [ inner lists  [ df opt result, dict opt_stats  ]  [...]  [...] ...]
    # indexing: opt_results[0][0]           -> first opt result data frame
    #           opt_results[0][0].iloc[:,0] -> first column in first opt_result data frame
    #           opt_results[0][1]           -> first opt result statistics dictionary
    #           opt_results[0][1]['key']    -> entry of 'key' in first opt result  
    #                                          statistics dictionary
    return opt_results

# 3) OP01_12h: maximize titer at 12h (s.t. maximum enzyme load; opt. vars: initial enzyme concentrations)
def solve_OP01_12h(model_path, fits_data_path, result_file_name):
    """ Solve a specific implementation of the optimization problem OP01 (maximize titer at 12h for a maximum enzyme load while reaching a certain minimum yield at 12h) with initial enzyme concentrations set as optimization variables. Repeat calculation for all parameter sets of a given parameter ensemble.

    :param model_path: path to the Copasi model object (.cps file)
    :type model_path: string
    :param fits_data_path: path to the parameter ensemble (.csv file)
    :type fits_data_path: string
    :param result_file_name: path of the calculated optimization result (.pkl file)
    :type fits_data_path: string

    :return opt_results: list of optimization results (elements: lists of result pandas data frames and dictionaries with information of the optimization calculation)
    :type opt_results: list
    """

    # LOAD MODEL OBJECT AND PARAMETER ENSEMBLE
    # the Copasi model file contains the structure of the model (stoichiometry, species, reactions, rate laws, ODEs); initial concentrations of all species, kinetic parameter values, and the specific optimization problem are set in this script
    model = load_model(model_path)
    # read the data of the parameter estimation sampling (the first column is all zeros, parameter columns start at .iloc[:,1])
    fits_data = pd.read_csv(fits_data_path)

    # SETUP OPTIMIZATION PROBLEM
    # set all species concentrations; ADP is naturally present in ATP stock solution at a rough ratio of 95% ATP/5% ADP but for the model simulation we assume an ideal 100% ATP stock solution (i.e., ADP is not present at t=0): the initial concentrations are from Tuan:Exp36_50mM (planned concentrations)
    set_species('ADP',          initial_concentration=1e-12)
    set_species('AMP',          initial_concentration=1e-12)
    set_species('ATP',          initial_concentration=0.5)
    set_species('ATPP',         initial_concentration=1e-12)
    set_species('E_GLMU',       initial_concentration=0.002032934)
    set_species('E_NAHK',       initial_concentration=0.002506077)
    set_species('E_PPA',        initial_concentration=0.001553358)
    set_species('E_PPK3',       initial_concentration=0.002878526)
    set_species('E_UDK',        initial_concentration=0.00410627)
    set_species('E_UMPK',       initial_concentration=0.004463289)
    set_species('GalNAc',       initial_concentration=50)
    set_species('GalNAc1P',     initial_concentration=1e-12)
    set_species('P',            initial_concentration=1e-12)
    set_species('PolyP',        initial_concentration=256.0)
    set_species('PP',           initial_concentration=1e-12)
    set_species('UDP',          initial_concentration=1e-12)
    set_species('UDP_GalNAc',   initial_concentration=1e-12)
    set_species('UMP',          initial_concentration=1e-12)
    set_species('UTP',          initial_concentration=1e-12)
    set_species('Uri',          initial_concentration=50)
    # create the necessary global quantities and events used for the target function (there can be no spaces in Copasi object names)
    add_parameter('UDP_GalNAc_at_12h', type='fixed', initial_value=0)
    add_event('save_UDP-GalNAc_at_12h', 'Time == 12', [('Values[UDP_GalNAc_at_12h]', '[UDP_GalNAc]')])
    # set target function
    set_objective_function(expression='Values[UDP_GalNAc_at_12h]', maximize=True)
    # setup optimization variables
    opt_vars = [{'name': '[E_GLMU]_0',
                'lower': 0.000203,
                'upper': 0.041,
                'start': 0.002032934},
                {'name': '[E_NAHK]_0',
                'lower': 0.000251,
                'upper': 0.05,
                'start': 0.002506077},
                {'name': '[E_PPA]_0',
                'lower': 0.000518,
                'upper': 0.103,
                'start': 0.001553358},
                {'name': '[E_PPK3]_0',
                'lower': 0.000288,
                'upper': 0.029,
                'start': 0.002878526},
                {'name': '[E_UDK]_0',
                'lower': 0.000411,
                'upper': 0.08,
                'start': 0.00410627},
                {'name': '[E_UMPK]_0',
                'lower': 0.000169,
                'upper': 0.034,
                'start': 0.004463289}]
    set_opt_parameters(opt_vars)
    # create the necessary global quantities used for the constraints (there can be no spaces in Copasi object names); enzyme concentrations in mM (mmol/l) are multiplied with their molecular weights (in g/mol divided by 1000 to get g/mmol) so that the resulting enzyme load is expressed in g/l
    add_parameter('E_tot_MW', type='assignment', expression='(49190/1000) * [E_GLMU]_0 + (39903/1000) * [E_NAHK]_0 + (19313/1000) * [E_PPA]_0 + (34740/1000) * [E_PPK3]_0 + (24353/1000) * [E_UDK]_0 + (22405/1000) * [E_UMPK]_0')
    # set constraints
    set_opt_constraints([
        {'name': 'Values[E_tot_MW]', 'lower': 0, 'upper': 0.53}
    ])
    # change time course settings (the optimization task calls the time course task)
    tc_settings = get_task_settings(T.TIME_COURSE)
    tc_settings['problem']['Duration'] = 12
    tc_settings['problem']['StepNumber'] = 100
    tc_settings['problem']['StepSize'] = 12/100
    set_task_settings(T.TIME_COURSE, settings=tc_settings)
    # set subtask and solver
    set_opt_settings(settings={'subtask': T.TIME_COURSE,
                               'method': {'name': PE.GENETIC_ALGORITHM}
                               })

    # SOLVE OPTIMIZATION PROBLEM FOR PARAMETER ENSEMBLE
    # setup for repeated optimizations with all different estimated parameter sets from random parameter sampling result
    opt_results = list()
    # loop over result data frame (every row = one parameter set)
    for index, param_set in tqdm(fits_data.iterrows()):
        # current list of parameter values
        curr_param_set_vals = list(param_set)
        # current list of parameter names; list() to only get content (values or names) of pandas series; param_set.keys() to get names of values in pandas series param_set
        curr_param_set_names = list(param_set.keys())
        # trim both lists (get rid of first entries = name and last entries = obj_val of fit)
        curr_param_set_vals = curr_param_set_vals[1:-1]
        curr_param_set_names = curr_param_set_names[1:-1]
        # set model parameter values to current estimated parameters from rand sampling result
        for i in range(len(list(curr_param_set_vals))): 
            set_reaction_parameters(name = curr_param_set_names[i], value = curr_param_set_vals[i])
        # optimize the model with these updated parameter values
        opt_result = run_optimization() # opt_result is a pandas.core.frame.DataFrame
        opt_stats = get_opt_statistic() # opt_stats is a dictionary
        # get list of optimization variable start values
        opt_vars_start_vals = [opt_vars[0]['start'], 
                               opt_vars[1]['start'], 
                               opt_vars[2]['start'], 
                               opt_vars[3]['start'], 
                               opt_vars[4]['start'], 
                               opt_vars[5]['start']]
        # add start values of optimization variables to result data frame as column 3
        opt_result.insert(2, "start", opt_vars_start_vals, True)
        # store the optimization result data frame and the opt_stats dictionary as a list
        opt_results.append([opt_result, opt_stats])
    # export opt_results as python object
    open_file = open(result_file_name, "wb")
    pickle.dump(opt_results, open_file)
    open_file.close()
    # result structure:
    # outer list opt_results [ inner lists  [ df opt result, dict opt_stats  ]  [...]  [...] ...]
    # indexing: opt_results[0][0]           -> first opt result data frame
    #           opt_results[0][0].iloc[:,0] -> first column in first opt_result data frame
    #           opt_results[0][1]           -> first opt result statistics dictionary
    #           opt_results[0][1]['key']    -> entry of 'key' in first opt result  
    #                                          statistics dictionary
    return opt_results

# 4) OP01withS_12h: maximize titer at 12h (s.t. maximum enzyme load; opt. vars: initial enzyme and substrate concentrations [main substrates capped to baseline concentrations])
def solve_OP01withS_12h(model_path, fits_data_path, result_file_name):
    """ Solve a specific implementation of the optimization problem OP01 (maximize titer at 12h for a maximum enzyme load while reaching a certain minimum yield at 12h) with initial enzyme and substrate concentrations set as optimization variables. Repeat calculation for all parameter sets of a given parameter ensemble.

    :param model_path: path to the Copasi model object (.cps file)
    :type model_path: string
    :param fits_data_path: path to the parameter ensemble (.csv file)
    :type fits_data_path: string
    :param result_file_name: path of the calculated optimization result (.pkl file)
    :type fits_data_path: string

    :return opt_results: list of optimization results (elements: lists of result pandas data frames and dictionaries with information of the optimization calculation)
    :type opt_results: list
    """

    # LOAD MODEL OBJECT AND PARAMETER ENSEMBLE
    # the Copasi model file contains the structure of the model (stoichiometry, species, reactions, rate laws, ODEs); initial concentrations of all species, kinetic parameter values, and the specific optimization problem are set in this script
    model = load_model(model_path)
    # read the data of the parameter estimation sampling (the first column is all zeros, parameter columns start at .iloc[:,1])
    fits_data = pd.read_csv(fits_data_path)

    # SETUP OPTIMIZATION PROBLEM
    # set all species concentrations; ADP is naturally present in ATP stock solution at a rough ratio of 95% ATP/5% ADP but for the model simulation we assume an ideal 100% ATP stock solution (i.e., ADP is not present at t=0): the initial concentrations are from Tuan:Exp36_50mM (planned concentrations)
    set_species('ADP',          initial_concentration=1e-12)
    set_species('AMP',          initial_concentration=1e-12)
    set_species('ATP',          initial_concentration=0.5)
    set_species('ATPP',         initial_concentration=1e-12)
    set_species('E_GLMU',       initial_concentration=0.002032934)
    set_species('E_NAHK',       initial_concentration=0.002506077)
    set_species('E_PPA',        initial_concentration=0.001553358)
    set_species('E_PPK3',       initial_concentration=0.002878526)
    set_species('E_UDK',        initial_concentration=0.00410627)
    set_species('E_UMPK',       initial_concentration=0.004463289)
    set_species('GalNAc',       initial_concentration=50)
    set_species('GalNAc1P',     initial_concentration=1e-12)
    set_species('P',            initial_concentration=1e-12)
    set_species('PolyP',        initial_concentration=256.0)
    set_species('PP',           initial_concentration=1e-12)
    set_species('UDP',          initial_concentration=1e-12)
    set_species('UDP_GalNAc',   initial_concentration=1e-12)
    set_species('UMP',          initial_concentration=1e-12)
    set_species('UTP',          initial_concentration=1e-12)
    set_species('Uri',          initial_concentration=50)
    # create the necessary global quantities and events used for the target function (there can be no spaces in Copasi object names)
    add_parameter('UDP_GalNAc_at_12h', type='fixed', initial_value=0)
    add_event('save_UDP-GalNAc_at_12h', 'Time == 12', [('Values[UDP_GalNAc_at_12h]', '[UDP_GalNAc]')])
    # set target function
    set_objective_function(expression='Values[UDP_GalNAc_at_12h]', maximize=True)
    # setup optimization variables
    opt_vars = [{'name': '[E_GLMU]_0',
                'lower': 0.000203,
                'upper': 0.041,
                'start': 0.002032934},
                {'name': '[E_NAHK]_0',
                'lower': 0.000251,
                'upper': 0.05,
                'start': 0.002506077},
                {'name': '[E_PPA]_0',
                'lower': 0.000518,
                'upper': 0.103,
                'start': 0.001553358},
                {'name': '[E_PPK3]_0',
                'lower': 0.000288,
                'upper': 0.029,
                'start': 0.002878526},
                {'name': '[E_UDK]_0',
                'lower': 0.000411,
                'upper': 0.08,
                'start': 0.00410627},
                {'name': '[E_UMPK]_0',
                'lower': 0.000169,
                'upper': 0.034,
                'start': 0.004463289},
                {'name': '[GalNAc]_0',
                'lower': 0.1,
                'upper': 50,
                'start': 50},
                {'name': '[Uri]_0',
                'lower': 0.1,
                'upper': 50,
                'start': 50}]
    set_opt_parameters(opt_vars)
    # create the necessary global quantities used for the constraints (there can be no spaces in Copasi object names); enzyme concentrations in mM (mmol/l) are multiplied with their molecular weights (in g/mol divided by 1000 to get g/mmol) so that the resulting enzyme load is expressed in g/l
    add_parameter('E_tot_MW', type='assignment', expression='(49190/1000) * [E_GLMU]_0 + (39903/1000) * [E_NAHK]_0 + (19313/1000) * [E_PPA]_0 + (34740/1000) * [E_PPK3]_0 + (24353/1000) * [E_UDK]_0 + (22405/1000) * [E_UMPK]_0')
    # set constraints
    set_opt_constraints([
        {'name': 'Values[E_tot_MW]', 'lower': 0, 'upper': 0.53}
    ])
    # change time course settings (the optimization task calls the time course task)
    tc_settings = get_task_settings(T.TIME_COURSE)
    tc_settings['problem']['Duration'] = 12
    tc_settings['problem']['StepNumber'] = 100
    tc_settings['problem']['StepSize'] = 12/100
    set_task_settings(T.TIME_COURSE, settings=tc_settings)
    # set subtask and solver
    set_opt_settings(settings={'subtask': T.TIME_COURSE,
                               'method': {'name': PE.GENETIC_ALGORITHM}
                               })

    # SOLVE OPTIMIZATION PROBLEM FOR PARAMETER ENSEMBLE
    # setup for repeated optimizations with all different estimated parameter sets from random parameter sampling result
    opt_results = list()
    # loop over result data frame (every row = one parameter set)
    for index, param_set in tqdm(fits_data.iterrows()):
        # current list of parameter values
        curr_param_set_vals = list(param_set)
        # current list of parameter names; list() to only get content (values or names) of pandas series; param_set.keys() to get names of values in pandas series param_set
        curr_param_set_names = list(param_set.keys())
        # trim both lists (get rid of first entries = name and last entries = obj_val of fit)
        curr_param_set_vals = curr_param_set_vals[1:-1]
        curr_param_set_names = curr_param_set_names[1:-1]
        # set model parameter values to current estimated parameters from rand sampling result
        for i in range(len(list(curr_param_set_vals))): 
            set_reaction_parameters(name = curr_param_set_names[i], value = curr_param_set_vals[i])
        # optimize the model with these updated parameter values
        opt_result = run_optimization() # opt_result is a pandas.core.frame.DataFrame
        opt_stats = get_opt_statistic() # opt_stats is a dictionary
        # get list of optimization variable start values
        opt_vars_start_vals = [opt_vars[0]['start'], 
                               opt_vars[1]['start'], 
                               opt_vars[2]['start'], 
                               opt_vars[3]['start'], 
                               opt_vars[4]['start'], 
                               opt_vars[5]['start'],
                               opt_vars[6]['start'],
                               opt_vars[7]['start']]
        # add start values of optimization variables to result data frame as column 3
        opt_result.insert(2, "start", opt_vars_start_vals, True)
        # store the optimization result data frame and the opt_stats dictionary as a list
        opt_results.append([opt_result, opt_stats])
    # export opt_results as python object
    open_file = open(result_file_name, "wb")
    pickle.dump(opt_results, open_file)
    open_file.close()
    # result structure:
    # outer list opt_results [ inner lists  [ df opt result, dict opt_stats  ]  [...]  [...] ...]
    # indexing: opt_results[0][0]           -> first opt result data frame
    #           opt_results[0][0].iloc[:,0] -> first column in first opt_result data frame
    #           opt_results[0][1]           -> first opt result statistics dictionary
    #           opt_results[0][1]['key']    -> entry of 'key' in first opt result  
    #                                          statistics dictionary
    return opt_results

# 5) OP05: minimize enzyme load (s.t. minimum titer; opt. vars: initial enzyme concentrations)
def solve_OP05(model_path, fits_data_path, result_file_name):
    """ Solve a specific implementation of the optimization problem OP05 (minimize enzyme load while reaching a certain minimum titer) with initial enzyme concentrations set as optimization variables. Repeat calculation for all parameter sets of a given parameter ensemble.

    :param model_path: path to the Copasi model object (.cps file)
    :type model_path: string
    :param fits_data_path: path to the parameter ensemble (.csv file)
    :type fits_data_path: string
    :param result_file_name: path of the calculated optimization result (.pkl file)
    :type fits_data_path: string

    :return opt_results: list of optimization results (elements: lists of result pandas data frames and dictionaries with information of the optimization calculation)
    :type opt_results: list
    """
    
    # LOAD MODEL OBJECT AND PARAMETER ENSEMBLE
    # the Copasi model file contains the structure of the model (stoichiometry, species, reactions, rate laws, ODEs); initial concentrations of all species, kinetic parameter values, and the specific optimization problem are set in this script
    model = load_model(model_path)
    # read the data of the parameter estimation sampling (the first column is all zeros, parameter columns start at .iloc[:,1])
    fits_data = pd.read_csv(fits_data_path)

    # SETUP OPTIMIZATION PROBLEM
    # set all species concentrations; ADP is naturally present in ATP stock solution at a rough ratio of 95% ATP/5% ADP but for the model simulation we assume an ideal 100% ATP stock solution (i.e., ADP is not present at t=0): the initial concentrations are from Tuan:Exp36_50mM (planned concentrations)
    set_species('ADP',          initial_concentration=1e-12)
    set_species('AMP',          initial_concentration=1e-12)
    set_species('ATP',          initial_concentration=0.5)
    set_species('ATPP',         initial_concentration=1e-12)
    set_species('E_GLMU',       initial_concentration=0.002032934)
    set_species('E_NAHK',       initial_concentration=0.002506077)
    set_species('E_PPA',        initial_concentration=0.001553358)
    set_species('E_PPK3',       initial_concentration=0.002878526)
    set_species('E_UDK',        initial_concentration=0.00410627)
    set_species('E_UMPK',       initial_concentration=0.004463289)
    set_species('GalNAc',       initial_concentration=50)
    set_species('GalNAc1P',     initial_concentration=1e-12)
    set_species('P',            initial_concentration=1e-12)
    set_species('PolyP',        initial_concentration=256.0)
    set_species('PP',           initial_concentration=1e-12)
    set_species('UDP',          initial_concentration=1e-12)
    set_species('UDP_GalNAc',   initial_concentration=1e-12)
    set_species('UMP',          initial_concentration=1e-12)
    set_species('UTP',          initial_concentration=1e-12)
    set_species('Uri',          initial_concentration=50)
    # create the necessary global quantities and events used for the target function (there can be no spaces in Copasi object names); enzyme concentrations in mM (mmol/l) are multiplied with their molecular weights (in g/mol divided by 1000 to get g/mmol) so that the resulting enzyme load is expressed in g/l
    add_parameter('E_tot_MW', type='assignment', expression='(49190/1000) * [E_GLMU]_0 + (39903/1000) * [E_NAHK]_0 + (19313/1000) * [E_PPA]_0 + (34740/1000) * [E_PPK3]_0 + (24353/1000) * [E_UDK]_0 + (22405/1000) * [E_UMPK]_0')
    # set target function
    set_objective_function(expression='Values[E_tot_MW]', minimize=True)
    # setup optimization variables
    opt_vars = [{'name': '[E_GLMU]_0',
                'lower': 0.000203,
                'upper': 0.041,
                'start': 0.002032934},
                {'name': '[E_NAHK]_0',
                'lower': 0.000251,
                'upper': 0.05,
                'start': 0.002506077},
                {'name': '[E_PPA]_0',
                'lower': 0.000518,
                'upper': 0.103,
                'start': 0.001553358},
                {'name': '[E_PPK3]_0',
                'lower': 0.000288,
                'upper': 0.029,
                'start': 0.002878526},
                {'name': '[E_UDK]_0',
                'lower': 0.000411,
                'upper': 0.08,
                'start': 0.00410627},
                {'name': '[E_UMPK]_0',
                'lower': 0.000169,
                'upper': 0.034,
                'start': 0.004463289}]
    set_opt_parameters(opt_vars)
    # create the necessary global quantities used for the constraints (there can be no spaces in Copasi object names)
    add_parameter('UDP_GalNAc_at_24h', type='fixed', initial_value=0)
    add_event('save_UDP-GalNAc_at_24h', 'Time == 24', [('Values[UDP_GalNAc_at_24h]', '[UDP_GalNAc]')])
    # set constraints; with fixed substrate concentrations an additional yield constraint is not necessary
    set_opt_constraints([
        {'name': 'Values[UDP_GalNAc_at_24h]', 'lower': 40.37, 'upper': 1000} # the upper limit is just a placeholder
    ])
    # change time course settings (the optimization task calls the time course task)
    tc_settings = get_task_settings(T.TIME_COURSE)
    tc_settings['problem']['Duration'] = 24
    tc_settings['problem']['StepNumber'] = 100
    tc_settings['problem']['StepSize'] = 24/100
    set_task_settings(T.TIME_COURSE, settings=tc_settings)
    # set subtask and solver
    set_opt_settings(settings={'subtask': T.TIME_COURSE,
                               'method': {'name': PE.GENETIC_ALGORITHM}
                               })

    # SOLVE OPTIMIZATION PROBLEM FOR PARAMETER ENSEMBLE
    # setup for repeated optimizations with all different estimated parameter sets from random parameter sampling result
    opt_results = list()
    # loop over result data frame (every row = one parameter set)
    for index, param_set in tqdm(fits_data.iterrows()):
        # current list of parameter values
        curr_param_set_vals = list(param_set)
        # current list of parameter names; list() to only get content (values or names) of pandas series; param_set.keys() to get names of values in pandas series param_set
        curr_param_set_names = list(param_set.keys())
        # trim both lists (get rid of first entries = name and last entries = obj_val of fit)
        curr_param_set_vals = curr_param_set_vals[1:-1]
        curr_param_set_names = curr_param_set_names[1:-1]
        # set model parameter values to current estimated parameters from rand sampling result
        for i in range(len(list(curr_param_set_vals))): 
            set_reaction_parameters(name = curr_param_set_names[i], value = curr_param_set_vals[i])
        # optimize the model with these updated parameter values
        opt_result = run_optimization() # opt_result is a pandas.core.frame.DataFrame
        opt_stats = get_opt_statistic() # opt_stats is a dictionary
        # get list of optimization variable start values
        opt_vars_start_vals = [opt_vars[0]['start'], 
                               opt_vars[1]['start'], 
                               opt_vars[2]['start'], 
                               opt_vars[3]['start'], 
                               opt_vars[4]['start'], 
                               opt_vars[5]['start']]
        # add start values of optimization variables to result data frame as column 3
        opt_result.insert(2, "start", opt_vars_start_vals, True)
        # store the optimization result data frame and the opt_stats dictionary as a list
        opt_results.append([opt_result, opt_stats])
    # export opt_results as python object
    open_file = open(result_file_name, "wb")
    pickle.dump(opt_results, open_file)
    open_file.close()
    # result structure:
    # outer list opt_results [ inner lists  [ df opt result, dict opt_stats  ]  [...]  [...] ...]
    # indexing: opt_results[0][0]           -> first opt result data frame
    #           opt_results[0][0].iloc[:,0] -> first column in first opt_result data frame
    #           opt_results[0][1]           -> first opt result statistics dictionary
    #           opt_results[0][1]['key']    -> entry of 'key' in first opt result  
    #                                          statistics dictionary
    return opt_results

# 6) OP05withS: minimize enzyme load (s.t. minimum titer and minimum yield; opt. vars: initial enzyme and substrate concentrations)
def solve_OP05withS(model_path, fits_data_path, result_file_name):
    """ Solve a specific implementation of the optimization problem OP05withS (minimize enzyme load while reaching a certain minimum titer and minimum yield) with initial enzyme and substrate concentrations set as optimization variables. Repeat calculation for all parameter sets of a given parameter ensemble.

    :param model_path: path to the Copasi model object (.cps file)
    :type model_path: string
    :param fits_data_path: path to the parameter ensemble (.csv file)
    :type fits_data_path: string
    :param result_file_name: path of the calculated optimization result (.pkl file)
    :type fits_data_path: string

    :return opt_results: list of optimization results (elements: lists of result pandas data frames and dictionaries with information of the optimization calculation)
    :type opt_results: list
    """

    # LOAD MODEL OBJECT AND PARAMETER ENSEMBLE
    # the Copasi model file contains the structure of the model (stoichiometry, species, reactions, rate laws, ODEs); initial concentrations of all species, kinetic parameter values, and the specific optimization problem are set in this script
    model = load_model(model_path)
    # read the data of the parameter estimation sampling (the first column is all zeros, parameter columns start at .iloc[:,1])
    fits_data = pd.read_csv(fits_data_path)

    # SETUP OPTIMIZATION PROBLEM
    # set all species concentrations; ADP is naturally present in ATP stock solution at a rough ratio of 95% ATP/5% ADP but for the model simulation we assume an ideal 100% ATP stock solution (i.e., ADP is not present at t=0): the initial concentrations are from Tuan:Exp36_50mM (planned concentrations)
    set_species('ADP',          initial_concentration=1e-12)
    set_species('AMP',          initial_concentration=1e-12)
    set_species('ATP',          initial_concentration=0.5)
    set_species('ATPP',         initial_concentration=1e-12)
    set_species('E_GLMU',       initial_concentration=0.002032934)
    set_species('E_NAHK',       initial_concentration=0.002506077)
    set_species('E_PPA',        initial_concentration=0.001553358)
    set_species('E_PPK3',       initial_concentration=0.002878526)
    set_species('E_UDK',        initial_concentration=0.00410627)
    set_species('E_UMPK',       initial_concentration=0.004463289)
    set_species('GalNAc',       initial_concentration=50)
    set_species('GalNAc1P',     initial_concentration=1e-12)
    set_species('P',            initial_concentration=1e-12)
    set_species('PolyP',        initial_concentration=256.0)
    set_species('PP',           initial_concentration=1e-12)
    set_species('UDP',          initial_concentration=1e-12)
    set_species('UDP_GalNAc',   initial_concentration=1e-12)
    set_species('UMP',          initial_concentration=1e-12)
    set_species('UTP',          initial_concentration=1e-12)
    set_species('Uri',          initial_concentration=50)
    # create the necessary global quantities and events used for the target function (there can be no spaces in Copasi object names); enzyme concentrations in mM (mmol/l) are multiplied with their molecular weights (in g/mol divided by 1000 to get g/mmol) so that the resulting enzyme load is expressed in g/l
    add_parameter('E_tot_MW', type='assignment', expression='(49190/1000) * [E_GLMU]_0 + (39903/1000) * [E_NAHK]_0 + (19313/1000) * [E_PPA]_0 + (34740/1000) * [E_PPK3]_0 + (24353/1000) * [E_UDK]_0 + (22405/1000) * [E_UMPK]_0')
    # set target function
    set_objective_function(expression='Values[E_tot_MW]', minimize=True)
    # setup optimization variables
    # setup optimization variables
    opt_vars = [{'name': '[E_GLMU]_0',
                'lower': 0.000203,
                'upper': 0.041,
                'start': 0.002032934},
                {'name': '[E_NAHK]_0',
                'lower': 0.000251,
                'upper': 0.05,
                'start': 0.002506077},
                {'name': '[E_PPA]_0',
                'lower': 0.000518,
                'upper': 0.103,
                'start': 0.001553358},
                {'name': '[E_PPK3]_0',
                'lower': 0.000288,
                'upper': 0.029,
                'start': 0.002878526},
                {'name': '[E_UDK]_0',
                'lower': 0.000411,
                'upper': 0.08,
                'start': 0.00410627},
                {'name': '[E_UMPK]_0',
                'lower': 0.000169,
                'upper': 0.034,
                'start': 0.004463289},
                {'name': '[GalNAc]_0',
                'lower': 0.1,
                'upper': 100,
                'start': 50},
                {'name': '[Uri]_0',
                'lower': 0.1,
                'upper': 100,
                'start': 50}]
    set_opt_parameters(opt_vars)
    # create the necessary global quantities used for the constraints (there can be no spaces in Copasi object names)
    add_parameter('UDP_GalNAc_at_24h', type='fixed', initial_value=0)
    add_event('save_UDP-GalNAc_at_24h', 'Time == 24', [('Values[UDP_GalNAc_at_24h]', '[UDP_GalNAc]')])
    add_parameter('Yield', type='assignment', expression='( 2 * Values[UDP_GalNAc_at_24h] ) / ( [GalNAc]_0 + [Uri]_0 )')
    # set constraints; with fixed substrate concentrations an additional yield constraint is not necessary
    set_opt_constraints([
        {'name': 'Values[UDP_GalNAc_at_24h]', 'lower': 40.37, 'upper': 1000}, # the upper limit is just a placeholder
        {'name': 'Values[Yield]', 'lower': 0.8274, 'upper': 1.0}
    ])
    # change time course settings (the optimization task calls the time course task)
    tc_settings = get_task_settings(T.TIME_COURSE)
    tc_settings['problem']['Duration'] = 24
    tc_settings['problem']['StepNumber'] = 100
    tc_settings['problem']['StepSize'] = 24/100
    set_task_settings(T.TIME_COURSE, settings=tc_settings)
    # set subtask and solver
    set_opt_settings(settings={'subtask': T.TIME_COURSE,
                               'method': {'name': PE.GENETIC_ALGORITHM}
                               })

    # SOLVE OPTIMIZATION PROBLEM FOR PARAMETER ENSEMBLE
    # setup for repeated optimizations with all different estimated parameter sets from random parameter sampling result
    opt_results = list()
    # loop over result data frame (every row = one parameter set)
    for index, param_set in tqdm(fits_data.iterrows()):
        # current list of parameter values
        curr_param_set_vals = list(param_set)
        # current list of parameter names; list() to only get content (values or names) of pandas series; param_set.keys() to get names of values in pandas series param_set
        curr_param_set_names = list(param_set.keys())
        # trim both lists (get rid of first entries = name and last entries = obj_val of fit)
        curr_param_set_vals = curr_param_set_vals[1:-1]
        curr_param_set_names = curr_param_set_names[1:-1]
        # set model parameter values to current estimated parameters from rand sampling result
        for i in range(len(list(curr_param_set_vals))): 
            set_reaction_parameters(name = curr_param_set_names[i], value = curr_param_set_vals[i])
        # optimize the model with these updated parameter values
        opt_result = run_optimization() # opt_result is a pandas.core.frame.DataFrame
        opt_stats = get_opt_statistic() # opt_stats is a dictionary
        # get list of optimization variable start values
        opt_vars_start_vals = [opt_vars[0]['start'], 
                               opt_vars[1]['start'], 
                               opt_vars[2]['start'], 
                               opt_vars[3]['start'], 
                               opt_vars[4]['start'], 
                               opt_vars[5]['start'],
                               opt_vars[6]['start'],
                               opt_vars[7]['start']]
        # add start values of optimization variables to result data frame as column 3
        opt_result.insert(2, "start", opt_vars_start_vals, True)
        # store the optimization result data frame and the opt_stats dictionary as a list
        opt_results.append([opt_result, opt_stats])
    # export opt_results as python object
    open_file = open(result_file_name, "wb")
    pickle.dump(opt_results, open_file)
    open_file.close()
    # result structure:
    # outer list opt_results [ inner lists  [ df opt result, dict opt_stats  ]  [...]  [...] ...]
    # indexing: opt_results[0][0]           -> first opt result data frame
    #           opt_results[0][0].iloc[:,0] -> first column in first opt_result data frame
    #           opt_results[0][1]           -> first opt result statistics dictionary
    #           opt_results[0][1]['key']    -> entry of 'key' in first opt result  
    #                                          statistics dictionary
    return opt_results

# ------------------------------------------------------------------------------------------------------ #

# CROSS-VALIDATION AND SCORING

def cross_val_and_score_OP01_7h(model_path, fits_data_path, opt_results_path, result_file_name):
    """ Select the best optimization result of OP01_7h via cross-validation: for each optimization result additional time course simulations are performed with the model parameters set to all of the remaining parameter sets of the ensemble that were not set for the respective optimization. For a parameter ensemble of size n=100 this means that for each of the 100 optimization results exactly 1 parameter set was used to calculate it leaving the other 99 parameter sets which are now used for additional simulations (resulting in a 100 by 100 matrix with different parameters sets per row and different predicted initial concentrations, i.e., optimization results per column which is visualized as a heat map). A scoring function is then applied to each column to identify the best optimization result.

    :param model_path: path to the Copasi model object (.cps file)
    :type model_path: string
    :param fits_data_path: path to the parameter ensemble (.csv file)
    :type fits_data_path: string
    :param opt_results_path: path to the list of optimization results (.pkl file; serialized list; elements: sublists of result pandas data frames and dictionaries with information of the optimization calculation)
    :type opt_results_path: string
    :param result_file_name: path of the calculated optimization result (.pkl file)
    :type fits_data_path: string

    :return df: cross-validation table 
    :type df: pandas dataframe
    :return best_opt_result: information on the best scoring optimization result
    :type best_opt_result: dictionary
    """

    # LOAD DATA
    # load model object
    model = load_model(model_path)
    # create the necessary global quantities and events
    add_parameter('UDP_GalNAc_at_7h', type='fixed', initial_value=0)
    add_event('save_UDP-GalNAc_at_7h', 'Time == 7', [('Values[UDP_GalNAc_at_7h]', '[UDP_GalNAc]')])
    add_parameter('E_tot_MW', type='assignment', expression='(49190/1000) * [E_GLMU]_0 + (39903/1000) * [E_NAHK]_0 + (19313/1000) * [E_PPA]_0 + (34740/1000) * [E_PPK3]_0 + (24353/1000) * [E_UDK]_0 + (22405/1000) * [E_UMPK]_0')
    # set all species concentrations that are not used as optimization variables to their respective constant values - all other species values are set by the optimization task; ADP is naturally present in ATP stock solution at a rough ratio of 95% ATP/5% ADP but for the model simulation we assume an ideal 100% ATP stock solution (i.e., ADP is not present at t=0)
    set_species('ADP',          initial_concentration=1e-12)
    set_species('AMP',          initial_concentration=1e-12)
    set_species('ATP',          initial_concentration=0.5)
    set_species('ATPP',         initial_concentration=1e-12)
    set_species('GalNAc',       initial_concentration=50)
    set_species('GalNAc1P',     initial_concentration=1e-12)
    set_species('P',            initial_concentration=1e-12)
    set_species('PolyP',        initial_concentration=256.0)
    set_species('PP',           initial_concentration=1e-12)
    set_species('UDP',          initial_concentration=1e-12)
    set_species('UDP_GalNAc',   initial_concentration=1e-12)
    set_species('UMP',          initial_concentration=1e-12)
    set_species('UTP',          initial_concentration=1e-12)
    set_species('Uri',          initial_concentration=50)
    # load the data of the parameter ensemble (the first column is all zeros, parameter columns start at .iloc[:,1])
    param_sets = pd.read_csv(fits_data_path)
    # load the data of the repeated optimization; return object: list of list of pandas data frames
    opt_result_sets = pd.read_pickle(opt_results_path)

    # CROSS_VALIDATION LOOP
    # check if the cross-validation has already been calclated; in this case just load the existing results and move to the scoring; otherwise calculate it here
    try:
        with open(f'{result_file_name}.pkl', 'rb') as f:
            all_titers = pickle.load(f)
            print(f'Cross-validation result loaded from {result_file_name}.pkl')
    except:
        # prepare list for overall loop result
        all_titers = []
        # loop over all i parameter sets and j optimization results
        for index, param_set in tqdm(param_sets.iterrows()):
            # prepare list for inner loop results (sublists with titers for all different optimization results for one parameter set)
            titers_row = []
            for opt_result in opt_result_sets:
                # 1) override model kinetic parameters with current parameter set
                # current list of parameter values
                curr_param_set_vals = list(param_set)
                # current list of parameter names; list() to only get content (values or names) of pandas series; param_set .keys() to get names of values in pandas series param_set; trim both lists (get rid of first entries = name and last entries = obj_val of fit)
                curr_param_set_names = list(param_set.keys())
                curr_param_set_vals = curr_param_set_vals[1:-1]
                curr_param_set_names = curr_param_set_names[1:-1]
                # set model parameter values to current estimated parameters from rand sampling result
                for i in range(len(list(curr_param_set_vals))): 
                    set_reaction_parameters(name = curr_param_set_names[i], value = curr_param_set_vals[i])
                # 2) override model initial enzyme concentrations with current optimization result
                curr_opt_set_conc = opt_result[0]['sol']    # pandas.series
                # first set_species needs to be called twice so that the species is actually set ... no idea why that is ...
                set_species('E_GLMU', initial_concentration=curr_opt_set_conc.iloc[0])
                set_species('E_GLMU', initial_concentration=curr_opt_set_conc.iloc[0])
                set_species('E_NAHK', initial_concentration=curr_opt_set_conc.iloc[1])
                set_species('E_PPA',  initial_concentration=curr_opt_set_conc.iloc[2])
                set_species('E_PPK3', initial_concentration=curr_opt_set_conc.iloc[3])
                set_species('E_UDK',  initial_concentration=curr_opt_set_conc.iloc[4])
                set_species('E_UMPK', initial_concentration=curr_opt_set_conc.iloc[5])
                # 3) simulate the model (with automatic step size control)
                sim_result = run_time_course(automatic=True, duration=7)
                # calculate the objective value (here: titer = UDP-GalNAc concentration at 7h); 7h is the last time point of the simulation
                sim_concs_timepoint_7h = sim_result.iloc[-1]
                titer_7h = sim_concs_timepoint_7h['UDP_GalNAc']
                # append the titer to the row of titers (= contains titers for all different optimization results and one parameter set)
                titers_row.append(titer_7h)
            # append row of titers for all different optimization results and one parameter set to overall result list
            all_titers.append(titers_row)
        # export all_titers as python object
        with open(f'{result_file_name}.pkl', 'wb') as f:
            pickle.dump(all_titers, f)
        # result structure:
        # outer list: elements are lists of titers [inner lists: elements are float titers]
        # indexing:  all_titers[0]           -> first list of titers (for all different optimization 
        #                                       results and for one parameter set)
        #            all_titers[0][0]        -> first titer of the inner list

    # SCORING
    # create pandas data frame from all titers list of lists
    # each sublist should be one row of the resulting df
    df = pd.DataFrame(all_titers)
    # column names: optimization result sets O1 ... Oi
    df.columns = ['O{}'.format(i+1) for i in df.columns]
    # row names: parameter sets p1 ... pj
    df.index = ['p{}'.format(i+1) for i in df.index]
    # calculate statistics (count, mean, std, min, 25%, 50%=median, 75%, max) for each column of the p x O matrix (rows: parameter sets p; columns: optimization result sets O)
    df_col_stats = df.describe()
    # calculate sum for each column
    df_col_sums= df.sum(axis=0)
    # select column with the highest median
    df_col_stats_medians = df_col_stats.loc['50%']
    df_max_col_median_idx = df_col_stats_medians.idxmax()
    # select column with the highest minimum
    df_col_stats_mins = df_col_stats.loc['min']
    df_max_col_min_idx = df_col_stats_mins.idxmax()
    # select column with highest sum
    df_max_col_sum_idx = df_col_sums.idxmax()
    # combine different statistical measures to obtain final ranking (idea: for each column calculate 'col_median/overall_max_median + col_min/overall_max_min + col_sum/overall_max_sum = score'; max score = 3)
    df_col_stats_max_median = max(df_col_stats.loc['50%'])     # highest median over all columns
    df_col_stats_max_min    = max(df_col_stats.loc['min'])     # highest min over all columns
    df_col_sums_max_sum     = max(df_col_sums)                 # highest sum over all columns
    # calculate score for each column (= each optimization result)
    scoring_results = []
    for col_name in df_col_stats:
        # content of current column in data frame 'df_col_stats'
        curr_col = df_col_stats.loc[:,col_name]
        # extract median and min
        curr_median = curr_col.loc['50%']
        curr_min = curr_col.loc['min']
        # extract sum from separate data frame 'df_col_sums' with same id (= current col name)
        curr_sum = df_col_sums.loc[col_name]
        # calculate score
        curr_score = (curr_median/df_col_stats_max_median) + (curr_min/df_col_stats_max_min) + (curr_sum/df_col_sums_max_sum)
        # save each term and score in sublist and append to overall result list
        scoring_results.append([curr_median/df_col_stats_max_median, 
                                curr_min/df_col_stats_max_min,
                                curr_sum/df_col_sums_max_sum,
                                curr_score])
    # create data frame from list of sublists 'scoring_result'
    # each sublist should be one row of the resulting df
    df_scores = pd.DataFrame(scoring_results)
    # column names: 'median_score', 'min_score', 'sum_score', 'total_score'
    df_scores.columns = ['median_score', 'min_score', 'sum_score', 'total_score']
    # row names: optimization result sets O1 ... Oi
    df_scores.index = ['O{}'.format(i+1) for i in df_scores.index]
    # identify optimization result with the best total score
    best_total_score_value = max(df_scores.loc[:,'total_score'])
    best_total_score_index = df_scores.loc[:,'total_score'].idxmax()

    # EXPORT
    # export the full scoring data frame
    with open(f'{result_file_name}_Scoring_Result_full.pkl', 'wb') as f:
        pickle.dump(df_scores, f)
    # gather all important information in a result dictionary
    best_opt_result = dict()
    best_opt_result.update({'name': best_total_score_index})
    best_opt_result.update({'index': list(df_scores.index).index(best_total_score_index)})
    best_opt_result.update({'score': df_scores.loc[best_total_score_index]})
    best_opt_result.update({'opt_result': opt_result_sets[list(df_scores.index).index(best_total_score_index)]})
    best_opt_result.update({'titers_data': df.loc[:, best_total_score_index]})
    best_opt_result.update({'titers_stats': df.loc[:, best_total_score_index].describe()})
    # convert predicted init conc. [mM] of best optimization result to [g/l] with molecular weights: (mM value/1000) [mol/l] * molecular weight [g/mol] = lab value [g/l]
    molecular_weights = pd.Series({'E_GLMU': 49190, 'E_NAHK': 39903, 'E_PPA': 19313, 'E_PPK3': 34740, 'E_UDK': 24353, 'E_UMPK': 22405})
    best_opt_result_pred_concs_mM = opt_result_sets[list(df_scores.index).index(best_total_score_index)][0].sol
    # create data frame for vector calculation (apply conversion function from [mM] to [g/l] to all entries); no molecular weights defined for substrate concentrations so NaN placeholders are introduced and those rows are then dropped
    calc_df = pd.concat([molecular_weights, best_opt_result_pred_concs_mM], axis=1).dropna()
    calc_df.columns = ['MW', 'mM_conc']
    calc_df['gram_per_litre_conc'] = (calc_df['mM_conc']/1000) * calc_df['MW']
    best_opt_result.update({'gram_per_litre_concs': calc_df.loc[:, 'gram_per_litre_conc']})
    # export result dictionary as serialized pickle file and write the dictionary to a text file
    file_name = f'{result_file_name}_Scoring_Result_{best_total_score_index}'
    with open(f'{file_name}.pkl', 'wb') as f:
        pickle.dump(best_opt_result, f)
    with open(f'{file_name}.txt', 'w') as f:
        for key, value in best_opt_result.items():
            f.write(f'#{key}#\n{value}\n')

    return df, best_opt_result

def cross_val_and_score_OP01withS_7h(model_path, fits_data_path, opt_results_path, result_file_name):
    """ Select the best optimization result of OP01withS_7h via cross-validation: for each optimization result additional time course simulations are performed with the model parameters set to all of the remaining parameter sets of the ensemble that were not set for the respective optimization. For a parameter ensemble of size n=100 this means that for each of the 100 optimization results exactly 1 parameter set was used to calculate it leaving the other 99 parameter sets which are now used for additional simulations (resulting in a 100 by 100 matrix with different parameters sets per row and different predicted initial concentrations, i.e., optimization results per column which is visualized as a heat map). A scoring function is then applied to each column to identify the best optimization result.

    :param model_path: path to the Copasi model object (.cps file)
    :type model_path: string
    :param fits_data_path: path to the parameter ensemble (.csv file)
    :type fits_data_path: string
    :param opt_results_path: path to the list of optimization results (.pkl file; serialized list; elements: sublists of result pandas data frames and dictionaries with information of the optimization calculation)
    :type opt_results_path: string
    :param result_file_name: path of the calculated optimization result (.pkl file)
    :type fits_data_path: string

    :return df: cross-validation table 
    :type df: pandas dataframe
    :return best_opt_result: information on the best scoring optimization result
    :type best_opt_result: dictionary
    """

    # LOAD DATA
    # load model object
    model = load_model(model_path)
    # create the necessary global quantities and events
    add_parameter('UDP_GalNAc_at_7h', type='fixed', initial_value=0)
    add_event('save_UDP-GalNAc_at_7h', 'Time == 7', [('Values[UDP_GalNAc_at_7h]', '[UDP_GalNAc]')])
    add_parameter('E_tot_MW', type='assignment', expression='(49190/1000) * [E_GLMU]_0 + (39903/1000) * [E_NAHK]_0 + (19313/1000) * [E_PPA]_0 + (34740/1000) * [E_PPK3]_0 + (24353/1000) * [E_UDK]_0 + (22405/1000) * [E_UMPK]_0')
    # set all species concentrations that are not used as optimization variables to their respective constant values - all other species values are set by the optimization task; ADP is naturally present in ATP stock solution at a rough ratio of 95% ATP/5% ADP but for the model simulation we assume an ideal 100% ATP stock solution (i.e., ADP is not present at t=0)
    set_species('ADP',          initial_concentration=1e-12)
    set_species('AMP',          initial_concentration=1e-12)
    set_species('ATP',          initial_concentration=0.5)
    set_species('ATPP',         initial_concentration=1e-12)
    set_species('GalNAc1P',     initial_concentration=1e-12)
    set_species('P',            initial_concentration=1e-12)
    set_species('PolyP',        initial_concentration=256.0)
    set_species('PP',           initial_concentration=1e-12)
    set_species('UDP',          initial_concentration=1e-12)
    set_species('UDP_GalNAc',   initial_concentration=1e-12)
    set_species('UMP',          initial_concentration=1e-12)
    set_species('UTP',          initial_concentration=1e-12)
    # load the data of the parameter ensemble (the first column is all zeros, parameter columns start at .iloc[:,1])
    param_sets = pd.read_csv(fits_data_path)
    # load the data of the repeated optimization; return object: list of list of pandas data frames
    opt_result_sets = pd.read_pickle(opt_results_path)

    # CROSS_VALIDATION LOOP
    # check if the cross-validation has already been calclated; in this case just load the existing results and move to the scoring; otherwise calculate it here
    try:
        with open(f'{result_file_name}.pkl', 'rb') as f:
            all_titers = pickle.load(f)
            print(f'Cross-validation result loaded from {result_file_name}.pkl')
    except:
        # prepare list for overall loop result
        all_titers = []
        # loop over all i parameter sets and j optimization results
        for index, param_set in tqdm(param_sets.iterrows()):
            # prepare list for inner loop results (sublists with titers for all different optimization results for one parameter set)
            titers_row = []
            for opt_result in opt_result_sets:
                # 1) override model kinetic parameters with current parameter set
                # current list of parameter values
                curr_param_set_vals = list(param_set)
                # current list of parameter names; list() to only get content (values or names) of pandas series; param_set .keys() to get names of values in pandas series param_set; trim both lists (get rid of first entries = name and last entries = obj_val of fit)
                curr_param_set_names = list(param_set.keys())
                curr_param_set_vals = curr_param_set_vals[1:-1]
                curr_param_set_names = curr_param_set_names[1:-1]
                # set model parameter values to current estimated parameters from rand sampling result
                for i in range(len(list(curr_param_set_vals))): 
                    set_reaction_parameters(name = curr_param_set_names[i], value = curr_param_set_vals[i])
                # 2) override model initial enzyme concentrations with current optimization result
                curr_opt_set_conc = opt_result[0]['sol']    # pandas.series
                # first set_species needs to be called twice so that the species is actually set ... no idea why that is ...
                set_species('E_GLMU', initial_concentration=curr_opt_set_conc.iloc[0])
                set_species('E_GLMU', initial_concentration=curr_opt_set_conc.iloc[0])
                set_species('E_NAHK', initial_concentration=curr_opt_set_conc.iloc[1])
                set_species('E_PPA',  initial_concentration=curr_opt_set_conc.iloc[2])
                set_species('E_PPK3', initial_concentration=curr_opt_set_conc.iloc[3])
                set_species('E_UDK',  initial_concentration=curr_opt_set_conc.iloc[4])
                set_species('E_UMPK', initial_concentration=curr_opt_set_conc.iloc[5])
                set_species('GalNAc', initial_concentration=curr_opt_set_conc.iloc[6])
                set_species('Uri',    initial_concentration=curr_opt_set_conc.iloc[7])
                # 3) simulate the model (with automatic step size control)
                sim_result = run_time_course(automatic=True, duration=7)
                # calculate the objective value (here: titer = UDP-GalNAc concentration at 7h); 7h is the last time point of the simulation
                sim_concs_timepoint_7h = sim_result.iloc[-1]
                titer_7h = sim_concs_timepoint_7h['UDP_GalNAc']
                # append the titer to the row of titers (= contains titers for all different optimization results and one parameter set)
                titers_row.append(titer_7h)
            # append row of titers for all different optimization results and one parameter set to overall result list
            all_titers.append(titers_row)
        # export all_titers as python object
        with open(f'{result_file_name}.pkl', 'wb') as f:
            pickle.dump(all_titers, f)
        # result structure:
        # outer list: elements are lists of titers [inner lists: elements are float titers]
        # indexing:  all_titers[0]           -> first list of titers (for all different optimization 
        #                                       results and for one parameter set)
        #            all_titers[0][0]        -> first titer of the inner list

    # SCORING
    # create pandas data frame from all titers list of lists
    # each sublist should be one row of the resulting df
    df = pd.DataFrame(all_titers)
    # column names: optimization result sets O1 ... Oi
    df.columns = ['O{}'.format(i+1) for i in df.columns]
    # row names: parameter sets p1 ... pj
    df.index = ['p{}'.format(i+1) for i in df.index]
    # calculate statistics (count, mean, std, min, 25%, 50%=median, 75%, max) for each column of the p x O matrix (rows: parameter sets p; columns: optimization result sets O)
    df_col_stats = df.describe()
    # calculate sum for each column
    df_col_sums= df.sum(axis=0)
    # select column with the highest median
    df_col_stats_medians = df_col_stats.loc['50%']
    df_max_col_median_idx = df_col_stats_medians.idxmax()
    # select column with the highest minimum
    df_col_stats_mins = df_col_stats.loc['min']
    df_max_col_min_idx = df_col_stats_mins.idxmax()
    # select column with highest sum
    df_max_col_sum_idx = df_col_sums.idxmax()
    # combine different statistical measures to obtain final ranking (idea: for each column calculate 'col_median/overall_max_median + col_min/overall_max_min + col_sum/overall_max_sum = score'; max score = 3)
    df_col_stats_max_median = max(df_col_stats.loc['50%'])     # highest median over all columns
    df_col_stats_max_min    = max(df_col_stats.loc['min'])     # highest min over all columns
    df_col_sums_max_sum     = max(df_col_sums)                 # highest sum over all columns
    # calculate score for each column (= each optimization result)
    scoring_results = []
    for col_name in df_col_stats:
        # content of current column in data frame 'df_col_stats'
        curr_col = df_col_stats.loc[:,col_name]
        # extract median and min
        curr_median = curr_col.loc['50%']
        curr_min = curr_col.loc['min']
        # extract sum from separate data frame 'df_col_sums' with same id (= current col name)
        curr_sum = df_col_sums.loc[col_name]
        # calculate score
        curr_score = (curr_median/df_col_stats_max_median) + (curr_min/df_col_stats_max_min) + (curr_sum/df_col_sums_max_sum)
        # save each term and score in sublist and append to overall result list
        scoring_results.append([curr_median/df_col_stats_max_median, 
                                curr_min/df_col_stats_max_min,
                                curr_sum/df_col_sums_max_sum,
                                curr_score])
    # create data frame from list of sublists 'scoring_result'
    # each sublist should be one row of the resulting df
    df_scores = pd.DataFrame(scoring_results)
    # column names: 'median_score', 'min_score', 'sum_score', 'total_score'
    df_scores.columns = ['median_score', 'min_score', 'sum_score', 'total_score']
    # row names: optimization result sets O1 ... Oi
    df_scores.index = ['O{}'.format(i+1) for i in df_scores.index]
    # identify optimization result with the best total score
    best_total_score_value = max(df_scores.loc[:,'total_score'])
    best_total_score_index = df_scores.loc[:,'total_score'].idxmax()

    # EXPORT
    # export the full scoring data frame
    with open(f'{result_file_name}_Scoring_Result_full.pkl', 'wb') as f:
        pickle.dump(df_scores, f)
    # gather all important information in a result dictionary
    best_opt_result = dict()
    best_opt_result.update({'name': best_total_score_index})
    best_opt_result.update({'index': list(df_scores.index).index(best_total_score_index)})
    best_opt_result.update({'score': df_scores.loc[best_total_score_index]})
    best_opt_result.update({'opt_result': opt_result_sets[list(df_scores.index).index(best_total_score_index)]})
    best_opt_result.update({'titers_data': df.loc[:, best_total_score_index]})
    best_opt_result.update({'titers_stats': df.loc[:, best_total_score_index].describe()})
    # convert predicted init conc. [mM] of best optimization result to [g/l] with molecular weights: (mM value/1000) [mol/l] * molecular weight [g/mol] = lab value [g/l]
    molecular_weights = pd.Series({'E_GLMU': 49190, 'E_NAHK': 39903, 'E_PPA': 19313, 'E_PPK3': 34740, 'E_UDK': 24353, 'E_UMPK': 22405})
    best_opt_result_pred_concs_mM = opt_result_sets[list(df_scores.index).index(best_total_score_index)][0].sol
    # create data frame for vector calculation (apply conversion function from [mM] to [g/l] to all entries); no molecular weights defined for substrate concentrations so NaN placeholders are introduced and those rows are then dropped
    calc_df = pd.concat([molecular_weights, best_opt_result_pred_concs_mM], axis=1).dropna()
    calc_df.columns = ['MW', 'mM_conc']
    calc_df['gram_per_litre_conc'] = (calc_df['mM_conc']/1000) * calc_df['MW']
    best_opt_result.update({'gram_per_litre_concs': calc_df.loc[:, 'gram_per_litre_conc']})
    # export result dictionary as serialized pickle file and write the dictionary to a text file
    file_name = f'{result_file_name}_Scoring_Result_{best_total_score_index}'
    with open(f'{file_name}.pkl', 'wb') as f:
        pickle.dump(best_opt_result, f)
    with open(f'{file_name}.txt', 'w') as f:
        for key, value in best_opt_result.items():
            f.write(f'#{key}#\n{value}\n')

    return df, best_opt_result

def cross_val_and_score_OP01_12h(model_path, fits_data_path, opt_results_path, result_file_name):
    """ Select the best optimization result of OP01_12h via cross-validation: for each optimization result additional time course simulations are performed with the model parameters set to all of the remaining parameter sets of the ensemble that were not set for the respective optimization. For a parameter ensemble of size n=100 this means that for each of the 100 optimization results exactly 1 parameter set was used to calculate it leaving the other 99 parameter sets which are now used for additional simulations (resulting in a 100 by 100 matrix with different parameters sets per row and different predicted initial concentrations, i.e., optimization results per column which is visualized as a heat map). A scoring function is then applied to each column to identify the best optimization result.

    :param model_path: path to the Copasi model object (.cps file)
    :type model_path: string
    :param fits_data_path: path to the parameter ensemble (.csv file)
    :type fits_data_path: string
    :param opt_results_path: path to the list of optimization results (.pkl file; serialized list; elements: sublists of result pandas data frames and dictionaries with information of the optimization calculation)
    :type opt_results_path: string
    :param result_file_name: path of the calculated optimization result (.pkl file)
    :type fits_data_path: string

    :return df: cross-validation table 
    :type df: pandas dataframe
    :return best_opt_result: information on the best scoring optimization result
    :type best_opt_result: dictionary
    """

    # LOAD DATA
    # load model object
    model = load_model(model_path)
    # create the necessary global quantities and events
    add_parameter('UDP_GalNAc_at_12h', type='fixed', initial_value=0)
    add_event('save_UDP-GalNAc_at_12h', 'Time == 12', [('Values[UDP_GalNAc_at_12h]', '[UDP_GalNAc]')])
    add_parameter('E_tot_MW', type='assignment', expression='(49190/1000) * [E_GLMU]_0 + (39903/1000) * [E_NAHK]_0 + (19313/1000) * [E_PPA]_0 + (34740/1000) * [E_PPK3]_0 + (24353/1000) * [E_UDK]_0 + (22405/1000) * [E_UMPK]_0')
    # set all species concentrations that are not used as optimization variables to their respective constant values - all other species values are set by the optimization task; ADP is naturally present in ATP stock solution at a rough ratio of 95% ATP/5% ADP but for the model simulation we assume an ideal 100% ATP stock solution (i.e., ADP is not present at t=0)
    set_species('ADP',          initial_concentration=1e-12)
    set_species('AMP',          initial_concentration=1e-12)
    set_species('ATP',          initial_concentration=0.5)
    set_species('ATPP',         initial_concentration=1e-12)
    set_species('GalNAc',       initial_concentration=50)
    set_species('GalNAc1P',     initial_concentration=1e-12)
    set_species('P',            initial_concentration=1e-12)
    set_species('PolyP',        initial_concentration=256.0)
    set_species('PP',           initial_concentration=1e-12)
    set_species('UDP',          initial_concentration=1e-12)
    set_species('UDP_GalNAc',   initial_concentration=1e-12)
    set_species('UMP',          initial_concentration=1e-12)
    set_species('UTP',          initial_concentration=1e-12)
    set_species('Uri',          initial_concentration=50)
    # load the data of the parameter ensemble (the first column is all zeros, parameter columns start at .iloc[:,1])
    param_sets = pd.read_csv(fits_data_path)
    # load the data of the repeated optimization; return object: list of list of pandas data frames
    opt_result_sets = pd.read_pickle(opt_results_path)

    # CROSS_VALIDATION LOOP
    # check if the cross-validation has already been calclated; in this case just load the existing results and move to the scoring; otherwise calculate it here
    try:
        with open(f'{result_file_name}.pkl', 'rb') as f:
            all_titers = pickle.load(f)
            print(f'Cross-validation result loaded from {result_file_name}.pkl')
    except:
        # prepare list for overall loop result
        all_titers = []
        # loop over all i parameter sets and j optimization results
        for index, param_set in tqdm(param_sets.iterrows()):
            # prepare list for inner loop results (sublists with titers for all different optimization results for one parameter set)
            titers_row = []
            for opt_result in opt_result_sets:
                # 1) override model kinetic parameters with current parameter set
                # current list of parameter values
                curr_param_set_vals = list(param_set)
                # current list of parameter names; list() to only get content (values or names) of pandas series; param_set .keys() to get names of values in pandas series param_set; trim both lists (get rid of first entries = name and last entries = obj_val of fit)
                curr_param_set_names = list(param_set.keys())
                curr_param_set_vals = curr_param_set_vals[1:-1]
                curr_param_set_names = curr_param_set_names[1:-1]
                # set model parameter values to current estimated parameters from rand sampling result
                for i in range(len(list(curr_param_set_vals))): 
                    set_reaction_parameters(name = curr_param_set_names[i], value = curr_param_set_vals[i])
                # 2) override model initial enzyme concentrations with current optimization result
                curr_opt_set_conc = opt_result[0]['sol']    # pandas.series
                # first set_species needs to be called twice so that the species is actually set ... no idea why that is ...
                set_species('E_GLMU', initial_concentration=curr_opt_set_conc.iloc[0])
                set_species('E_GLMU', initial_concentration=curr_opt_set_conc.iloc[0])
                set_species('E_NAHK', initial_concentration=curr_opt_set_conc.iloc[1])
                set_species('E_PPA',  initial_concentration=curr_opt_set_conc.iloc[2])
                set_species('E_PPK3', initial_concentration=curr_opt_set_conc.iloc[3])
                set_species('E_UDK',  initial_concentration=curr_opt_set_conc.iloc[4])
                set_species('E_UMPK', initial_concentration=curr_opt_set_conc.iloc[5])
                # 3) simulate the model (with automatic step size control)
                sim_result = run_time_course(automatic=True, duration=12)
                # calculate the objective value (here: titer = UDP-GalNAc concentration at 12h); 12h is the last time point of the simulation
                sim_concs_timepoint_12h = sim_result.iloc[-1]
                titer_12h = sim_concs_timepoint_12h['UDP_GalNAc']
                # append the titer to the row of titers (= contains titers for all different optimization results and one parameter set)
                titers_row.append(titer_12h)
            # append row of titers for all different optimization results and one parameter set to overall result list
            all_titers.append(titers_row)
        # export all_titers as python object
        with open(f'{result_file_name}.pkl', 'wb') as f:
            pickle.dump(all_titers, f)
        # result structure:
        # outer list: elements are lists of titers [inner lists: elements are float titers]
        # indexing:  all_titers[0]           -> first list of titers (for all different optimization 
        #                                       results and for one parameter set)
        #            all_titers[0][0]        -> first titer of the inner list

    # SCORING
    # create pandas data frame from all titers list of lists
    # each sublist should be one row of the resulting df
    df = pd.DataFrame(all_titers)
    # column names: optimization result sets O1 ... Oi
    df.columns = ['O{}'.format(i+1) for i in df.columns]
    # row names: parameter sets p1 ... pj
    df.index = ['p{}'.format(i+1) for i in df.index]
    # calculate statistics (count, mean, std, min, 25%, 50%=median, 75%, max) for each column of the p x O matrix (rows: parameter sets p; columns: optimization result sets O)
    df_col_stats = df.describe()
    # calculate sum for each column
    df_col_sums= df.sum(axis=0)
    # select column with the highest median
    df_col_stats_medians = df_col_stats.loc['50%']
    df_max_col_median_idx = df_col_stats_medians.idxmax()
    # select column with the highest minimum
    df_col_stats_mins = df_col_stats.loc['min']
    df_max_col_min_idx = df_col_stats_mins.idxmax()
    # select column with highest sum
    df_max_col_sum_idx = df_col_sums.idxmax()
    # combine different statistical measures to obtain final ranking (idea: for each column calculate 'col_median/overall_max_median + col_min/overall_max_min + col_sum/overall_max_sum = score'; max score = 3)
    df_col_stats_max_median = max(df_col_stats.loc['50%'])     # highest median over all columns
    df_col_stats_max_min    = max(df_col_stats.loc['min'])     # highest min over all columns
    df_col_sums_max_sum     = max(df_col_sums)                 # highest sum over all columns
    # calculate score for each column (= each optimization result)
    scoring_results = []
    for col_name in df_col_stats:
        # content of current column in data frame 'df_col_stats'
        curr_col = df_col_stats.loc[:,col_name]
        # extract median and min
        curr_median = curr_col.loc['50%']
        curr_min = curr_col.loc['min']
        # extract sum from separate data frame 'df_col_sums' with same id (= current col name)
        curr_sum = df_col_sums.loc[col_name]
        # calculate score
        curr_score = (curr_median/df_col_stats_max_median) + (curr_min/df_col_stats_max_min) + (curr_sum/df_col_sums_max_sum)
        # save each term and score in sublist and append to overall result list
        scoring_results.append([curr_median/df_col_stats_max_median, 
                                curr_min/df_col_stats_max_min,
                                curr_sum/df_col_sums_max_sum,
                                curr_score])
    # create data frame from list of sublists 'scoring_result'
    # each sublist should be one row of the resulting df
    df_scores = pd.DataFrame(scoring_results)
    # column names: 'median_score', 'min_score', 'sum_score', 'total_score'
    df_scores.columns = ['median_score', 'min_score', 'sum_score', 'total_score']
    # row names: optimization result sets O1 ... Oi
    df_scores.index = ['O{}'.format(i+1) for i in df_scores.index]
    # identify optimization result with the best total score
    best_total_score_value = max(df_scores.loc[:,'total_score'])
    best_total_score_index = df_scores.loc[:,'total_score'].idxmax()

    # EXPORT
    # export the full scoring data frame
    with open(f'{result_file_name}_Scoring_Result_full.pkl', 'wb') as f:
        pickle.dump(df_scores, f)
    # gather all important information in a result dictionary
    best_opt_result = dict()
    best_opt_result.update({'name': best_total_score_index})
    best_opt_result.update({'index': list(df_scores.index).index(best_total_score_index)})
    best_opt_result.update({'score': df_scores.loc[best_total_score_index]})
    best_opt_result.update({'opt_result': opt_result_sets[list(df_scores.index).index(best_total_score_index)]})
    best_opt_result.update({'titers_data': df.loc[:, best_total_score_index]})
    best_opt_result.update({'titers_stats': df.loc[:, best_total_score_index].describe()})
    # convert predicted init conc. [mM] of best optimization result to [g/l] with molecular weights: (mM value/1000) [mol/l] * molecular weight [g/mol] = lab value [g/l]
    molecular_weights = pd.Series({'E_GLMU': 49190, 'E_NAHK': 39903, 'E_PPA': 19313, 'E_PPK3': 34740, 'E_UDK': 24353, 'E_UMPK': 22405})
    best_opt_result_pred_concs_mM = opt_result_sets[list(df_scores.index).index(best_total_score_index)][0].sol
    # create data frame for vector calculation (apply conversion function from [mM] to [g/l] to all entries); no molecular weights defined for substrate concentrations so NaN placeholders are introduced and those rows are then dropped
    calc_df = pd.concat([molecular_weights, best_opt_result_pred_concs_mM], axis=1).dropna()
    calc_df.columns = ['MW', 'mM_conc']
    calc_df['gram_per_litre_conc'] = (calc_df['mM_conc']/1000) * calc_df['MW']
    best_opt_result.update({'gram_per_litre_concs': calc_df.loc[:, 'gram_per_litre_conc']})
    # export result dictionary as serialized pickle file and write the dictionary to a text file
    file_name = f'{result_file_name}_Scoring_Result_{best_total_score_index}'
    with open(f'{file_name}.pkl', 'wb') as f:
        pickle.dump(best_opt_result, f)
    with open(f'{file_name}.txt', 'w') as f:
        for key, value in best_opt_result.items():
            f.write(f'#{key}#\n{value}\n')

    return df, best_opt_result

def cross_val_and_score_OP01withS_12h(model_path, fits_data_path, opt_results_path, result_file_name):
    """ Select the best optimization result of OP01withS_12h via cross-validation: for each optimization result additional time course simulations are performed with the model parameters set to all of the remaining parameter sets of the ensemble that were not set for the respective optimization. For a parameter ensemble of size n=100 this means that for each of the 100 optimization results exactly 1 parameter set was used to calculate it leaving the other 99 parameter sets which are now used for additional simulations (resulting in a 100 by 100 matrix with different parameters sets per row and different predicted initial concentrations, i.e., optimization results per column which is visualized as a heat map). A scoring function is then applied to each column to identify the best optimization result.

    :param model_path: path to the Copasi model object (.cps file)
    :type model_path: string
    :param fits_data_path: path to the parameter ensemble (.csv file)
    :type fits_data_path: string
    :param opt_results_path: path to the list of optimization results (.pkl file; serialized list; elements: sublists of result pandas data frames and dictionaries with information of the optimization calculation)
    :type opt_results_path: string
    :param result_file_name: path of the calculated optimization result (.pkl file)
    :type fits_data_path: string

    :return df: cross-validation table 
    :type df: pandas dataframe
    :return best_opt_result: information on the best scoring optimization result
    :type best_opt_result: dictionary
    """

    # LOAD DATA
    # load model object
    model = load_model(model_path)
    # create the necessary global quantities and events
    add_parameter('UDP_GalNAc_at_12h', type='fixed', initial_value=0)
    add_event('save_UDP-GalNAc_at_12h', 'Time == 12', [('Values[UDP_GalNAc_at_12h]', '[UDP_GalNAc]')])
    add_parameter('E_tot_MW', type='assignment', expression='(49190/1000) * [E_GLMU]_0 + (39903/1000) * [E_NAHK]_0 + (19313/1000) * [E_PPA]_0 + (34740/1000) * [E_PPK3]_0 + (24353/1000) * [E_UDK]_0 + (22405/1000) * [E_UMPK]_0')
    # set all species concentrations that are not used as optimization variables to their respective constant values - all other species values are set by the optimization task; ADP is naturally present in ATP stock solution at a rough ratio of 95% ATP/5% ADP but for the model simulation we assume an ideal 100% ATP stock solution (i.e., ADP is not present at t=0)
    set_species('ADP',          initial_concentration=1e-12)
    set_species('AMP',          initial_concentration=1e-12)
    set_species('ATP',          initial_concentration=0.5)
    set_species('ATPP',         initial_concentration=1e-12)
    set_species('GalNAc1P',     initial_concentration=1e-12)
    set_species('P',            initial_concentration=1e-12)
    set_species('PolyP',        initial_concentration=256.0)
    set_species('PP',           initial_concentration=1e-12)
    set_species('UDP',          initial_concentration=1e-12)
    set_species('UDP_GalNAc',   initial_concentration=1e-12)
    set_species('UMP',          initial_concentration=1e-12)
    set_species('UTP',          initial_concentration=1e-12)
    # load the data of the parameter ensemble (the first column is all zeros, parameter columns start at .iloc[:,1])
    param_sets = pd.read_csv(fits_data_path)
    # load the data of the repeated optimization; return object: list of list of pandas data frames
    opt_result_sets = pd.read_pickle(opt_results_path)

    # CROSS_VALIDATION LOOP
    # check if the cross-validation has already been calclated; in this case just load the existing results and move to the scoring; otherwise calculate it here
    try:
        with open(f'{result_file_name}.pkl', 'rb') as f:
            all_titers = pickle.load(f)
            print(f'Cross-validation result loaded from {result_file_name}.pkl')
    except:
        # prepare list for overall loop result
        all_titers = []
        # loop over all i parameter sets and j optimization results
        for index, param_set in tqdm(param_sets.iterrows()):
            # prepare list for inner loop results (sublists with titers for all different optimization results for one parameter set)
            titers_row = []
            for opt_result in opt_result_sets:
                # 1) override model kinetic parameters with current parameter set
                # current list of parameter values
                curr_param_set_vals = list(param_set)
                # current list of parameter names; list() to only get content (values or names) of pandas series; param_set .keys() to get names of values in pandas series param_set; trim both lists (get rid of first entries = name and last entries = obj_val of fit)
                curr_param_set_names = list(param_set.keys())
                curr_param_set_vals = curr_param_set_vals[1:-1]
                curr_param_set_names = curr_param_set_names[1:-1]
                # set model parameter values to current estimated parameters from rand sampling result
                for i in range(len(list(curr_param_set_vals))): 
                    set_reaction_parameters(name = curr_param_set_names[i], value = curr_param_set_vals[i])
                # 2) override model initial enzyme concentrations with current optimization result
                curr_opt_set_conc = opt_result[0]['sol']    # pandas.series
                # first set_species needs to be called twice so that the species is actually set ... no idea why that is ...
                set_species('E_GLMU', initial_concentration=curr_opt_set_conc.iloc[0])
                set_species('E_GLMU', initial_concentration=curr_opt_set_conc.iloc[0])
                set_species('E_NAHK', initial_concentration=curr_opt_set_conc.iloc[1])
                set_species('E_PPA',  initial_concentration=curr_opt_set_conc.iloc[2])
                set_species('E_PPK3', initial_concentration=curr_opt_set_conc.iloc[3])
                set_species('E_UDK',  initial_concentration=curr_opt_set_conc.iloc[4])
                set_species('E_UMPK', initial_concentration=curr_opt_set_conc.iloc[5])
                set_species('GalNAc', initial_concentration=curr_opt_set_conc.iloc[6])
                set_species('Uri',    initial_concentration=curr_opt_set_conc.iloc[7])
                # 3) simulate the model (with automatic step size control)
                sim_result = run_time_course(automatic=True, duration=12)
                # calculate the objective value (here: titer = UDP-GalNAc concentration at 12h); 12h is the last time point of the simulation
                sim_concs_timepoint_12h = sim_result.iloc[-1]
                titer_12h = sim_concs_timepoint_12h['UDP_GalNAc']
                # append the titer to the row of titers (= contains titers for all different optimization results and one parameter set)
                titers_row.append(titer_12h)
            # append row of titers for all different optimization results and one parameter set to overall result list
            all_titers.append(titers_row)
        # export all_titers as python object
        with open(f'{result_file_name}.pkl', 'wb') as f:
            pickle.dump(all_titers, f)
        # result structure:
        # outer list: elements are lists of titers [inner lists: elements are float titers]
        # indexing:  all_titers[0]           -> first list of titers (for all different optimization 
        #                                       results and for one parameter set)
        #            all_titers[0][0]        -> first titer of the inner list

    # SCORING
    # create pandas data frame from all titers list of lists
    # each sublist should be one row of the resulting df
    df = pd.DataFrame(all_titers)
    # column names: optimization result sets O1 ... Oi
    df.columns = ['O{}'.format(i+1) for i in df.columns]
    # row names: parameter sets p1 ... pj
    df.index = ['p{}'.format(i+1) for i in df.index]
    # calculate statistics (count, mean, std, min, 25%, 50%=median, 75%, max) for each column of the p x O matrix (rows: parameter sets p; columns: optimization result sets O)
    df_col_stats = df.describe()
    # calculate sum for each column
    df_col_sums= df.sum(axis=0)
    # select column with the highest median
    df_col_stats_medians = df_col_stats.loc['50%']
    df_max_col_median_idx = df_col_stats_medians.idxmax()
    # select column with the highest minimum
    df_col_stats_mins = df_col_stats.loc['min']
    df_max_col_min_idx = df_col_stats_mins.idxmax()
    # select column with highest sum
    df_max_col_sum_idx = df_col_sums.idxmax()
    # combine different statistical measures to obtain final ranking (idea: for each column calculate 'col_median/overall_max_median + col_min/overall_max_min + col_sum/overall_max_sum = score'; max score = 3)
    df_col_stats_max_median = max(df_col_stats.loc['50%'])     # highest median over all columns
    df_col_stats_max_min    = max(df_col_stats.loc['min'])     # highest min over all columns
    df_col_sums_max_sum     = max(df_col_sums)                 # highest sum over all columns
    # calculate score for each column (= each optimization result)
    scoring_results = []
    for col_name in df_col_stats:
        # content of current column in data frame 'df_col_stats'
        curr_col = df_col_stats.loc[:,col_name]
        # extract median and min
        curr_median = curr_col.loc['50%']
        curr_min = curr_col.loc['min']
        # extract sum from separate data frame 'df_col_sums' with same id (= current col name)
        curr_sum = df_col_sums.loc[col_name]
        # calculate score
        curr_score = (curr_median/df_col_stats_max_median) + (curr_min/df_col_stats_max_min) + (curr_sum/df_col_sums_max_sum)
        # save each term and score in sublist and append to overall result list
        scoring_results.append([curr_median/df_col_stats_max_median, 
                                curr_min/df_col_stats_max_min,
                                curr_sum/df_col_sums_max_sum,
                                curr_score])
    # create data frame from list of sublists 'scoring_result'
    # each sublist should be one row of the resulting df
    df_scores = pd.DataFrame(scoring_results)
    # column names: 'median_score', 'min_score', 'sum_score', 'total_score'
    df_scores.columns = ['median_score', 'min_score', 'sum_score', 'total_score']
    # row names: optimization result sets O1 ... Oi
    df_scores.index = ['O{}'.format(i+1) for i in df_scores.index]
    # identify optimization result with the best total score
    best_total_score_value = max(df_scores.loc[:,'total_score'])
    best_total_score_index = df_scores.loc[:,'total_score'].idxmax()

    # EXPORT
    # export the full scoring data frame
    with open(f'{result_file_name}_Scoring_Result_full.pkl', 'wb') as f:
        pickle.dump(df_scores, f)
    # gather all important information in a result dictionary
    best_opt_result = dict()
    best_opt_result.update({'name': best_total_score_index})
    best_opt_result.update({'index': list(df_scores.index).index(best_total_score_index)})
    best_opt_result.update({'score': df_scores.loc[best_total_score_index]})
    best_opt_result.update({'opt_result': opt_result_sets[list(df_scores.index).index(best_total_score_index)]})
    best_opt_result.update({'titers_data': df.loc[:, best_total_score_index]})
    best_opt_result.update({'titers_stats': df.loc[:, best_total_score_index].describe()})
    # convert predicted init conc. [mM] of best optimization result to [g/l] with molecular weights: (mM value/1000) [mol/l] * molecular weight [g/mol] = lab value [g/l]
    molecular_weights = pd.Series({'E_GLMU': 49190, 'E_NAHK': 39903, 'E_PPA': 19313, 'E_PPK3': 34740, 'E_UDK': 24353, 'E_UMPK': 22405})
    best_opt_result_pred_concs_mM = opt_result_sets[list(df_scores.index).index(best_total_score_index)][0].sol
    # create data frame for vector calculation (apply conversion function from [mM] to [g/l] to all entries); no molecular weights defined for substrate concentrations so NaN placeholders are introduced and those rows are then dropped
    calc_df = pd.concat([molecular_weights, best_opt_result_pred_concs_mM], axis=1).dropna()
    calc_df.columns = ['MW', 'mM_conc']
    calc_df['gram_per_litre_conc'] = (calc_df['mM_conc']/1000) * calc_df['MW']
    best_opt_result.update({'gram_per_litre_concs': calc_df.loc[:, 'gram_per_litre_conc']})
    # export result dictionary as serialized pickle file and write the dictionary to a text file
    file_name = f'{result_file_name}_Scoring_Result_{best_total_score_index}'
    with open(f'{file_name}.pkl', 'wb') as f:
        pickle.dump(best_opt_result, f)
    with open(f'{file_name}.txt', 'w') as f:
        for key, value in best_opt_result.items():
            f.write(f'#{key}#\n{value}\n')

    return df, best_opt_result

def cross_val_and_score_OP05(model_path, fits_data_path, opt_results_path, result_file_name):
    """ Select the best optimization result of OP05 via cross-validation: for each optimization result additional time course simulations are performed with the model parameters set to all of the remaining parameter sets of the ensemble that were not set for the respective optimization. For a parameter ensemble of size n=100 this means that for each of the 100 optimization results exactly 1 parameter set was used to calculate it leaving the other 99 parameter sets which are now used for additional simulations (resulting in a 100 by 100 matrix with different parameters sets per row and different predicted initial concentrations, i.e., optimization results per column which is visualized as a heat map). A scoring function is then applied to each column to identify the best optimization result.

    :param model_path: path to the Copasi model object (.cps file)
    :type model_path: string
    :param fits_data_path: path to the parameter ensemble (.csv file)
    :type fits_data_path: string
    :param opt_results_path: path to the list of optimization results (.pkl file; serialized list; elements: sublists of result pandas data frames and dictionaries with information of the optimization calculation)
    :type opt_results_path: string
    :param result_file_name: path of the calculated optimization result (.pkl file)
    :type fits_data_path: string

    :return df_E_tot: cross-validation table of enzyme load values
    :type df_E_tot: pandas dataframe
    :return df_titers: cross-validation table of titer values 
    :type df_titers: pandas dataframe
    :return best_opt_result: information on the best scoring optimization result
    :type best_opt_result: dictionary
    """

    # LOAD DATA
    # load model object
    model = load_model(model_path)
    # create necessary global quantities and events
    add_parameter('E_tot_MW', type='assignment', expression='(49190/1000) * [E_GLMU]_0 + (39903/1000) * [E_NAHK]_0 + (19313/1000) * [E_PPA]_0 + (34740/1000) * [E_PPK3]_0 + (24353/1000) * [E_UDK]_0 + (22405/1000) * [E_UMPK]_0')
    add_parameter('UDP_GalNAc_at_24h', type='fixed', initial_value=0)
    add_event('save_UDP-GalNAc_at_24h', 'Time == 24', [('Values[UDP_GalNAc_at_24h]', '[UDP_GalNAc]')])
    add_parameter('Uri_Yield', type='assignment', expression='Values[UDP_GalNAc_at_24h] / [Uri]_0')
    # set all species concentrations that are not used as optimization variables to their respective constant values - all other species values are set by the optimization task; ADP is naturally present in ATP stock solution at a rough ratio of 95% ATP/5% ADP but for the model simulation we assume an ideal 100% ATP stock solution (i.e., ADP is not present at t=0)
    set_species('ADP',          initial_concentration=1e-12)
    set_species('AMP',          initial_concentration=1e-12)
    set_species('ATP',          initial_concentration=0.5)
    set_species('ATPP',         initial_concentration=1e-12)
    set_species('GalNAc',       initial_concentration=50)
    set_species('GalNAc1P',     initial_concentration=1e-12)
    set_species('P',            initial_concentration=1e-12)
    set_species('PolyP',        initial_concentration=256.0)
    set_species('PP',           initial_concentration=1e-12)
    set_species('UDP',          initial_concentration=1e-12)
    set_species('UDP_GalNAc',   initial_concentration=1e-12)
    set_species('UMP',          initial_concentration=1e-12)
    set_species('UTP',          initial_concentration=1e-12)
    set_species('Uri',          initial_concentration=50)
    # load the data of the parameter ensemble (the first column is all zeros, parameter columns start at .iloc[:,1])
    param_sets = pd.read_csv(fits_data_path)
    # load the data of the repeated optimization; return object: list of list of pandas data frames
    opt_result_sets = pd.read_pickle(opt_results_path)

    # CROSS-VALIDATION LOOP
    # check if the cross-validation has already been calclated; in this case just load the existing results and move to the scoring; otherwise calculate it here
    try:
        with open(f'{result_file_name}_allTit.pkl', 'rb') as f:
            all_titers = pickle.load(f)
        with open(f'{result_file_name}_allEtot.pkl', 'rb') as f:
            all_E_tot = pickle.load(f)
        print(f'Cross-validation results loaded from {result_file_name}_allTit.pkl and {result_file_name}_allEtot.pkl')
    except:
        # prepare lists for overall loop results
        all_E_tot = []
        all_titers = []
        # loop over all i parameter sets and j optimization results
        for index, param_set in tqdm(param_sets.iterrows()):
            # prepare lists for inner loop results (sublists with enzyme loads and titers for all different optimization results for one parameter set)
            E_tot_row = []
            titer_row = []
            for opt_result in opt_result_sets:
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
                # 2) override model initial enzyme concentrations with current optimization result
                curr_opt_set_conc = opt_result[0]['sol']    # pandas.series
                # first set_species needs to be called twice so that the species is actually set ... no idea why that is ...
                set_species('E_GLMU', initial_concentration=curr_opt_set_conc.iloc[0])
                set_species('E_GLMU', initial_concentration=curr_opt_set_conc.iloc[0])
                set_species('E_NAHK', initial_concentration=curr_opt_set_conc.iloc[1])
                set_species('E_PPA',  initial_concentration=curr_opt_set_conc.iloc[2])
                set_species('E_PPK3', initial_concentration=curr_opt_set_conc.iloc[3])
                set_species('E_UDK',  initial_concentration=curr_opt_set_conc.iloc[4])
                set_species('E_UMPK', initial_concentration=curr_opt_set_conc.iloc[5])
                # 3) simulate the model
                sim_result = run_time_course(automatic=True, duration=24)
                # calculate the objective value (here: enzyme load = sum of all initial enzyme concentrations in mM; E_tot is constant across all time points)
                sim_vals_timepoint_24h = sim_result.iloc[-1]
                E_tot_24h = sim_vals_timepoint_24h['Values[E_tot_MW]']
                # titer at 24h is part of a constraint in this case so it is also stored
                titer_24h = sim_vals_timepoint_24h['UDP_GalNAc']
                # append the enzyme load to the row of enzyme loads (= contains enzyme loads for all different optimization results and one parameter set)
                E_tot_row.append(E_tot_24h)
                # append the titer to the row of titers (= contains titers for all different optimization results and one parameter set)
                titer_row.append(titer_24h)
            # append row of enzyme loads and titers for all different optimization results and one parameter set to overall result list
            all_E_tot.append(E_tot_row)
            all_titers.append(titer_row)
        # export all_E_tot as python object
        open_file = open(result_file_name+'_allEtot.pkl', "wb")
        pickle.dump(all_E_tot, open_file)
        open_file.close()
        # export all_E_tot as python object
        with open(f'{result_file_name}_allEtot.pkl', 'wb') as f:
            pickle.dump(all_E_tot, f)
        # export all_titers as python object
        with open(f'{result_file_name}_allTit.pkl', 'wb') as f:
            pickle.dump(all_titers, f)
    # result structure (example for all enzyme loads object):
    # outer list: elements are lists of enzyme loads [inner lists: elements are enzyme loads]
    # indexing:  all_E_tot[0]           -> first list of enzyme loads (for all different optimization 
    #                                       results and for one parameter set)
    #            all_E_tot[0][0]        -> first enzyme load of the inner list

    # SCORING
    # create pandas data frame from all enzyme loads list of lists
    # each sublist should be one row of the resulting df
    df_E_tot = pd.DataFrame(all_E_tot)
    # column names: optimization result sets O1 ... Oi
    df_E_tot.columns = ['O{}'.format(i+1) for i in df_E_tot.columns]
    # row names: parameter sets p1 ... pj
    df_E_tot.index = ['p{}'.format(i+1) for i in df_E_tot.index]
    # create pandas data frame from all titers list of lists
    # each sublist should be one row of the resulting df
    df_titers = pd.DataFrame(all_titers)
    # column names: optimization result sets O1 ... Oi
    df_titers.columns = ['O{}'.format(i+1) for i in df_titers.columns]
    # row names: parameter sets p1 ... pj
    df_titers.index = ['p{}'.format(i+1) for i in df_titers.index]
    # calculate statistics (count, mean, std, min, 25%, 50%=median, 75%, max) for each column of the p x O matrix (rows: parameter sets p; columns: optimization result sets O)
    df_E_tot_col_stats = df_E_tot.describe()
    df_titers_col_stats = df_titers.describe()
    # calculate sum for each column
    df_E_tot_col_sums= df_E_tot.sum(axis=0)
    df_titers_col_sums= df_titers.sum(axis=0)
    # combine different statistical measures to obtain final ranking
    # idea: (col_enyzme_load_value-0.53)/(overall_enzyme_load_min-0.53) + min(1,(titers_curr_median/40.37)) because constraints in optimization problem were titer at 24h >= 40.37 mM and E_tot <= 0.53 g/l; all terms are normalized between [0,1]: closer to 1 means better (lower) E_tot and better (higher) titer; titer term is capped at 1 since we only care that the constraint is fulfilled (even higher titers are not relevant for the scoring of this particular optimization problem where the focus is on low E_tot values -> we don't want to select an optimization result where bad E_tot values are compensated for by high titers)
    df_E_tot_col_stats_min_min  = min(df_E_tot_col_stats.loc['min'])   # lowest min enzyme load
    # calculate score for each column (= each optimization result)
    scoring_results = []
    for E_tot_col_name, titers_col_name in zip(df_titers_col_stats, df_E_tot_col_stats):
        # content of current column in data frame 'df_titers_col_stats'
        titers_curr_col = df_titers_col_stats.loc[:,titers_col_name]
        # extract median
        titers_curr_median = titers_curr_col.loc['50%']
        # content of current column in data frame 'df_E_tot_col_stats'
        E_tot_curr_col = df_E_tot_col_stats.loc[:,E_tot_col_name]
        E_tot_curr_min = E_tot_curr_col.loc['min']
        # calculate score
        curr_score = (E_tot_curr_min-0.53)/(df_E_tot_col_stats_min_min-0.53) + min(1,(titers_curr_median/40.37))
        # save each term and score in sublist and append to overall result list
        scoring_results.append([(E_tot_curr_min-0.53)/(df_E_tot_col_stats_min_min-0.53),
                                min(1,titers_curr_median/40.37),
                                curr_score])
    # create data frame from list of sublists 'scoring_result'
    # each sublist should be one row of the resulting df
    df_scores = pd.DataFrame(scoring_results)
    # column names: 'E_tot_min_score', 'titers_min_constraint_score', 'total_score'
    df_scores.columns = ['E_tot_min_score', 'titers_median_constraint_score', 'total_score']
    # row names: optimization result sets O1 ... Oi
    df_scores.index = ['O{}'.format(i+1) for i in df_scores.index]
    # identify optimization result with the best total score
    best_total_score_value = max(df_scores.loc[:,'total_score'])
    best_total_score_index = df_scores.loc[:,'total_score'].idxmax()

    # EXPORT
    # export the full scoring data frame
    with open(f'{result_file_name}_Scoring_Result_full.pkl', 'wb') as f:
        pickle.dump(df_scores, f)
    # gather all important information in a result dictionary
    best_opt_result = dict()
    best_opt_result.update({'name': best_total_score_index})
    best_opt_result.update({'index': list(df_scores.index).index(best_total_score_index)})
    best_opt_result.update({'score': df_scores.loc[best_total_score_index]})
    best_opt_result.update({'opt_result': opt_result_sets[list(df_scores.index).index(best_total_score_index)]})
    best_opt_result.update({'Etot_data': df_E_tot.loc[:, best_total_score_index]})
    best_opt_result.update({'Etot_stats': df_E_tot.loc[:, best_total_score_index].describe()})
    best_opt_result.update({'titers_data': df_titers.loc[:, best_total_score_index]})
    best_opt_result.update({'titers_stats': df_titers.loc[:, best_total_score_index].describe()})
    # convert predicted init conc. [mM] of best optimization result to [g/l] with molecular weights: (mM value/1000) [mol/l] * molecular weight [g/mol] = lab value [g/l]
    molecular_weights = pd.Series({'E_GLMU': 49190, 'E_NAHK': 39903, 'E_PPA': 19313, 'E_PPK3': 34740, 'E_UDK': 24353, 'E_UMPK': 22405})
    best_opt_result_pred_concs_mM = opt_result_sets[list(df_scores.index).index(best_total_score_index)][0].sol
    # create data frame for vector calculation (apply conversion function from [mM] to [g/l] to all entries); no molecular weights defined for substrate concentrations so NaN placeholders are introduced and those rows are then dropped
    calc_df = pd.concat([molecular_weights, best_opt_result_pred_concs_mM], axis=1).dropna()
    calc_df.columns = ['MW', 'mM_conc']
    calc_df['gram_per_litre_conc'] = (calc_df['mM_conc']/1000) * calc_df['MW']
    best_opt_result.update({'gram_per_litre_concs': calc_df.loc[:, 'gram_per_litre_conc']})
    # export result dictionary as serialized pickle file and write the dictionary to a text file
    file_name = f'{result_file_name}_Scoring_Result_{best_total_score_index}'
    with open(f'{file_name}.pkl', 'wb') as f:
        pickle.dump(best_opt_result, f)
    with open(f'{file_name}.txt', 'w') as f:
        for key, value in best_opt_result.items():
            f.write(f'#{key}#\n{value}\n')

    return df_E_tot, df_titers, best_opt_result

def cross_val_and_score_OP05withS(model_path, fits_data_path, opt_results_path, result_file_name):
    """ Select the best optimization result of OP05withS via cross-validation: for each optimization result additional time course simulations are performed with the model parameters set to all of the remaining parameter sets of the ensemble that were not set for the respective optimization. For a parameter ensemble of size n=100 this means that for each of the 100 optimization results exactly 1 parameter set was used to calculate it leaving the other 99 parameter sets which are now used for additional simulations (resulting in a 100 by 100 matrix with different parameters sets per row and different predicted initial concentrations, i.e., optimization results per column which is visualized as a heat map). A scoring function is then applied to each column to identify the best optimization result.

    :param model_path: path to the Copasi model object (.cps file)
    :type model_path: string
    :param fits_data_path: path to the parameter ensemble (.csv file)
    :type fits_data_path: string
    :param opt_results_path: path to the list of optimization results (.pkl file; serialized list; elements: sublists of result pandas data frames and dictionaries with information of the optimization calculation)
    :type opt_results_path: string
    :param result_file_name: path of the calculated optimization result (.pkl file)
    :type fits_data_path: string

    :return df_E_tot: cross-validation table of enzyme load values
    :type df_E_tot: pandas dataframe
    :return df_titers: cross-validation table of titer values 
    :type df_titers: pandas dataframe
    :return df_yields: cross-validation table of yield values
    :type df_yields: pandas dataframe
    :return best_opt_result: information on the best scoring optimization result
    :type best_opt_result: dictionary
    """

    # LOAD DATA
    # load model object
    model = load_model(model_path)
    # create necessary global quantities and events
    add_parameter('E_tot_MW', type='assignment', expression='(49190/1000) * [E_GLMU]_0 + (39903/1000) * [E_NAHK]_0 + (19313/1000) * [E_PPA]_0 + (34740/1000) * [E_PPK3]_0 + (24353/1000) * [E_UDK]_0 + (22405/1000) * [E_UMPK]_0')
    add_parameter('UDP_GalNAc_at_24h', type='fixed', initial_value=0)
    add_event('save_UDP-GalNAc_at_24h', 'Time == 24', [('Values[UDP_GalNAc_at_24h]', '[UDP_GalNAc]')])
    add_parameter('Yield', type='assignment', expression='( 2 * Values[UDP_GalNAc_at_24h] ) / ( [GalNAc]_0 + [Uri]_0 )')
    # set  all species concentrations that are not used as optimization variables to their respective constant values - all other species values are set by the optimization task; ADP is naturally present in ATP stock solution at a rough ratio of 95% ATP/5% ADP but for the model simulation we assume an ideal 100% ATP stock solution (i.e., ADP is not present at t=0)
    set_species('ADP',          initial_concentration=1e-12)
    set_species('AMP',          initial_concentration=1e-12)
    set_species('ATP',          initial_concentration=0.5)
    set_species('ATPP',         initial_concentration=1e-12)
    set_species('GalNAc1P',     initial_concentration=1e-12)
    set_species('P',            initial_concentration=1e-12)
    set_species('PolyP',        initial_concentration=256.0)
    set_species('PP',           initial_concentration=1e-12)
    set_species('UDP',          initial_concentration=1e-12)
    set_species('UDP_GalNAc',   initial_concentration=1e-12)
    set_species('UMP',          initial_concentration=1e-12)
    set_species('UTP',          initial_concentration=1e-12)
    # load the data of the parameter ensemble (the first column is all zeros, parameter columns start at .iloc[:,1])
    param_sets = pd.read_csv(fits_data_path)
    # load the data of the repeated optimization; return object: list of list of pandas data frames
    opt_result_sets = pd.read_pickle(opt_results_path)

    # CROSS-VALIDATION LOOP
    # check if the cross-validation has already been calclated; in this case just load the existing results and move to the scoring; otherwise calculate it here
    try:
        with open(f'{result_file_name}_allEtot.pkl', 'rb') as f:
            all_E_tot = pickle.load(f)
        with open(f'{result_file_name}_allTit.pkl', 'rb') as f:
            all_titers = pickle.load(f)
        with open(f'{result_file_name}_allYields.pkl', 'rb') as f:
            all_yields = pickle.load(f)
        print(f'Cross-validation results loaded from {result_file_name}_allEtot.pkl, {result_file_name}_allTit.pkl, and {result_file_name}_allYields.pkl')
    except:
        # prepare lists for overall loop results
        all_E_tot = []
        all_titers = []
        all_yields = []
        # loop over all i parameter sets and j optimization results
        for index, param_set in tqdm(param_sets.iterrows()):
            # prepare lists for inner loop results (sublists with enzyme loads and titers for all different optimization results for one parameter set)
            E_tot_row = []
            titer_row = []
            yields_row = []
            for opt_result in opt_result_sets:
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
                # 2) override model initial enzyme concentrations with current optimization result
                curr_opt_set_conc = opt_result[0]['sol']    # pandas.series
                # first set_species needs to be called twice so that the species is actually set ... no idea why that is ...
                set_species('E_GLMU', initial_concentration=curr_opt_set_conc.iloc[0])
                set_species('E_GLMU', initial_concentration=curr_opt_set_conc.iloc[0])
                set_species('E_NAHK', initial_concentration=curr_opt_set_conc.iloc[1])
                set_species('E_PPA',  initial_concentration=curr_opt_set_conc.iloc[2])
                set_species('E_PPK3', initial_concentration=curr_opt_set_conc.iloc[3])
                set_species('E_UDK',  initial_concentration=curr_opt_set_conc.iloc[4])
                set_species('E_UMPK', initial_concentration=curr_opt_set_conc.iloc[5])
                set_species('GalNAc', initial_concentration=curr_opt_set_conc.iloc[6])
                set_species('Uri',    initial_concentration=curr_opt_set_conc.iloc[7])
                # 3) simulate the model
                sim_result = run_time_course(automatic=True, duration=24)
                # calculate the objective value (here: enzyme load = sum of all initial enzyme concentrations in mM; E_tot is constant across all time points)
                sim_vals_timepoint_24h = sim_result.iloc[-1]
                E_tot_24h = sim_vals_timepoint_24h['Values[E_tot_MW]']
                # titer and yield at 24h are part of the constraints so they are also stored
                titer_24h = sim_vals_timepoint_24h['UDP_GalNAc']
                yield_24h = sim_vals_timepoint_24h['Values[Yield]']
                # append the enzyme load to the row of enzyme loads (= contains enzyme loads for all different optimization results and one parameter set)
                E_tot_row.append(E_tot_24h)
                # append the titer to the row of titers (= contains titers for all different optimization results and one parameter set)
                titer_row.append(titer_24h)
                # append the yield to the row of yields (= contains yields for all different optimization results and one parameter set)
                yields_row.append(yield_24h)
            # append row of enzyme loads and titers for all different optimization results and one parameter set to overall result list
            all_E_tot.append(E_tot_row)
            all_titers.append(titer_row)
            all_yields.append(yields_row)
        # export all_E_tot as python object
        with open(f'{result_file_name}_allEtot.pkl', 'wb') as f:
             pickle.dump(all_E_tot, f)
        # export all_titers as python object
        with open(f'{result_file_name}_allTit.pkl', 'wb') as f:
             pickle.dump(all_titers, f)
        # export all_yields as python object
        with open(f'{result_file_name}_allYields.pkl', 'wb') as f:
             pickle.dump(all_yields, f)
        # result structure (example for all enzyme loads object):
        # outer list: elements are lists of enzyme loads [inner lists: elements are enzyme loads]
        # indexing:  all_E_tot[0]           -> first list of enzyme loads (for all different optimization 
        #                                       results and for one parameter set)
        #            all_E_tot[0][0]        -> first enzyme load of the inner list

    # SCORING
    # create pandas data frame from all enzyme loads list of lists
    # each sublist should be one row of the resulting df
    df_E_tot = pd.DataFrame(all_E_tot)
    # column names: optimization result sets O1 ... Oi
    df_E_tot.columns = ['O{}'.format(i+1) for i in df_E_tot.columns]
    # row names: parameter sets p1 ... pj
    df_E_tot.index = ['p{}'.format(i+1) for i in df_E_tot.index]
    # create pandas data frame from all titers list of lists
    # each sublist should be one row of the resulting df
    df_titers = pd.DataFrame(all_titers)
    # column names: optimization result sets O1 ... Oi
    df_titers.columns = ['O{}'.format(i+1) for i in df_titers.columns]
    # row names: parameter sets p1 ... pj
    df_titers.index = ['p{}'.format(i+1) for i in df_titers.index]
    # create pandas data frame from all yields list of lists
    # each sublist should be one row of the resulting df
    df_yields = pd.DataFrame(all_yields)
    # column names: optimization result sets O1 ... Oi
    df_yields.columns = ['O{}'.format(i+1) for i in df_yields.columns]
    # row names: parameter sets p1 ... pj
    df_yields.index = ['p{}'.format(i+1) for i in df_yields.index]
    # calculate statistics (count, mean, std, min, 25%, 50%=median, 75%, max) for each column of the p x O matrix (rows: parameter sets p; columns: optimization result sets O)
    df_E_tot_col_stats = df_E_tot.describe()
    df_titers_col_stats = df_titers.describe()
    df_yields_col_stats = df_yields.describe()
    # calculate sum for each column
    df_E_tot_col_sums= df_E_tot.sum(axis=0)
    df_titers_col_sums= df_titers.sum(axis=0)
    df_yields_col_sums= df_yields.sum(axis=0)
    # combine different statistical measures to obtain final ranking
    # idea: (col_enyzme_load_value-0.53)/(overall_enzyme_load_min-0.53) + min(1,(titers_curr_median/40.37)) + min(1,(yields_curr_median/0.8274)) because constraints in optimization problem were titer at 24h >= 40.37 mM, yield at 24h >= 82.74% and E_tot <= 0.53 g/l; all three terms are normalized between [0,1]: closer to 1 means better (lower) E_tot, better (higher) titer and better (higher) yield; titer and yield terms are capped at 1 since we only care that the constraints are fulfilled (even higher titers and yields are not relevant for the scoring of this particular optimization problem where the focus is on low E_tot values -> we don't want to select an optimization result where bad E_tot values are compensated for by high titers and yields)
    df_E_tot_col_stats_min_min  = min(df_E_tot_col_stats.loc['min'])   # lowest min enzyme load
    # calculate score for each column (= each optimization result)
    scoring_results = []
    for E_tot_col_name, titers_col_name, yields_col_name in zip(df_E_tot_col_stats, df_titers_col_stats, df_yields_col_stats):
        # content of current column in data frame 'df_E_tot_col_stats'
        E_tot_curr_col = df_E_tot_col_stats.loc[:,E_tot_col_name]
        E_tot_curr_min = E_tot_curr_col.loc['min']
        # content of current column in data frame 'df_titers_col_stats'
        titers_curr_col = df_titers_col_stats.loc[:,titers_col_name]
        # extract median
        titers_curr_median = titers_curr_col.loc['50%']
        # content of current column in data frame 'df_yields_col_stats'
        yields_curr_col = df_yields_col_stats.loc[:,yields_col_name]
        # extract median
        yields_curr_median = yields_curr_col.loc['50%']
        # calculate score
        curr_score = (E_tot_curr_min-0.53)/(df_E_tot_col_stats_min_min-0.53) + min(1,(titers_curr_median/40.37)) + min(1,(yields_curr_median/0.8274))
        # save each term and score in sublist and append to overall result list
        scoring_results.append([(E_tot_curr_min-0.53)/(df_E_tot_col_stats_min_min-0.53),
                                min(1,titers_curr_median/40.37),
                                min(1,yields_curr_median/0.8274),
                                curr_score])
    # create data frame from list of sublists 'scoring_result'
    # each sublist should be one row of the resulting df
    df_scores = pd.DataFrame(scoring_results)
    # column names: 'E_tot_min_score', 'titers_min_constraint_score', 'total_score'
    df_scores.columns = ['E_tot_min_score', 'titers_median_constraint_score', 'yields_median_constraint_score', 'total_score']
    # row names: optimization result sets O1 ... Oi
    df_scores.index = ['O{}'.format(i+1) for i in df_scores.index]
    # identify optimization result with the best total score
    best_total_score_value = max(df_scores.loc[:,'total_score'])
    best_total_score_index = df_scores.loc[:,'total_score'].idxmax()

    # EXPORT
    # export the full scoring data frame
    with open(f'{result_file_name}_Scoring_Result_full.pkl', 'wb') as f:
        pickle.dump(df_scores, f)
    # gather all important information in a result dictionary
    best_opt_result = dict()
    best_opt_result.update({'name': best_total_score_index})
    best_opt_result.update({'index': list(df_scores.index).index(best_total_score_index)})
    best_opt_result.update({'score': df_scores.loc[best_total_score_index]})
    best_opt_result.update({'opt_result': opt_result_sets[list(df_scores.index).index(best_total_score_index)]})
    best_opt_result.update({'Etot_data': df_E_tot.loc[:, best_total_score_index]})
    best_opt_result.update({'Etot_stats': df_E_tot.loc[:, best_total_score_index].describe()})
    best_opt_result.update({'titers_data': df_titers.loc[:, best_total_score_index]})
    best_opt_result.update({'titers_stats': df_titers.loc[:, best_total_score_index].describe()})
    best_opt_result.update({'yields_data': df_yields.loc[:, best_total_score_index]})
    best_opt_result.update({'yields_stats': df_yields.loc[:, best_total_score_index].describe()})
    # convert predicted init conc. [mM] of best optimization result to [g/l] with molecular weights: (mM value/1000) [mol/l] * molecular weight [g/mol] = lab value [g/l]
    molecular_weights = pd.Series({'E_GLMU': 49190, 'E_NAHK': 39903, 'E_PPA': 19313, 'E_PPK3': 34740, 'E_UDK': 24353, 'E_UMPK': 22405})
    best_opt_result_pred_concs_mM = opt_result_sets[list(df_scores.index).index(best_total_score_index)][0].sol
    # create data frame for vector calculation (apply conversion function from [mM] to [g/l] to all entries); no molecular weights defined for substrate concentrations so NaN placeholders are introduced and those rows are then dropped
    calc_df = pd.concat([molecular_weights, best_opt_result_pred_concs_mM], axis=1).dropna()
    calc_df.columns = ['MW', 'mM_conc']
    calc_df['gram_per_litre_conc'] = (calc_df['mM_conc']/1000) * calc_df['MW']
    best_opt_result.update({'gram_per_litre_concs': calc_df.loc[:, 'gram_per_litre_conc']})
    # export result dictionary as serialized pickle file and write the dictionary to a text file
    file_name = f'{result_file_name}_Scoring_Result_{best_total_score_index}'
    with open(f'{file_name}.pkl', 'wb') as f:
        pickle.dump(best_opt_result, f)
    with open(f'{file_name}.txt', 'w') as f:
        for key, value in best_opt_result.items():
            f.write(f'#{key}#\n{value}\n')

    return df_E_tot, df_titers, df_yields, best_opt_result

