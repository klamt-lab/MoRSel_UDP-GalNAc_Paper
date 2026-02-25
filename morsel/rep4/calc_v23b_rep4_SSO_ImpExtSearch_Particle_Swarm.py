#!/usr/bin/env python
# -*- coding: utf-8 -*-

import time
import os.path
import sys

sys.path.append("..")
from func_lib import *


# ------------------------------------------------------------------------------------------------------ #

# NEGATIVE CONTROL VARIANT PARAMETER ENSEMBLE
# -------------------------------------------
# Idea: as all methods in this script file are fundamentally based on testing the inclusion of various model modifications (e.g., regulation terms) it is necessary to first establish a baseline as a negative control: it only includes reactions that are part of the cascade design (and implements them with basic convenience kinetics rate laws without any extra modifications and/or regulations)
# define negative control baseline variant
# set model name
model_variant_name = 'v23brep4NegControl_Exp36.37.39'
n_runs = 50
# calculate a negative control ensemble only if it doesn't already exist
if os.path.isfile(f"{model_variant_name}_{n_runs}runs_evaluated_vars_log") == False:
    #               base term1
    neg_ctrl_var = [0,   0,   # ADP_Decay_v1
                    0,   0,   # ADP_Decay_v2
                    1,  10,   # GLMU
                    1,  10,   # NAHK_ATP
                    0,  10,   # NAHK_ADP
                    0,  10,   # NAHK_ATPP
                    1,  10,   # PPA
                    1,  10,   # PPK3_A
                    1,  10,   # PPK3_U
                    0,  10,   # PPK3_tetra
                    1,  10,   # UDK_ATP
                    0,  10,   # UDK_ADP
                    0,  10,   # UDK_ATPP
                    1,  10,   # UMPK_ATP
                    0,  10,   # UMPK_ADP
                    0,  10]   # UMPK_ATPP
    neg_ctrl_var = pd.DataFrame(np.array(neg_ctrl_var).reshape(16, 2),
                                columns=['base', 'term1'],
                                index=['ADP_Decay_v1', 'ADP_Decay_v2', 'GLMU', 'NAHK_ATP', 'NAHK_ADP', 'NAHK_ATPP', 'PPA', 'PPK3_A', 'PPK3_U', 'PPK3_tetra', 'UDK_ATP', 'UDK_ADP', 'UDK_ATPP', 'UMPK_ATP', 'UMPK_ADP', 'UMPK_ATPP'])
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
    # repeat the evaluation of the same model variant n times (n replicates) and store the results
    evaluated_vars_log = create_parameter_ensemble(neg_ctrl_var, n_runs, exp_data_file_names, exp_data_dataframes, model_variant_name)
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
#            base term1 term2 term3 term4
start_var = [0,   0,    0,    0,    0,    # ADP_Decay_v1
             0,   0,    0,    0,    0,    # ADP_Decay_v2
             1,  10,    0,    0,    0,    # GLMU
             1,  10,    0,    0,    0,    # NAHK_ATP
             0,  10,    0,    0,    0,    # NAHK_ADP
             0,  10,    0,    0,    0,    # NAHK_ATPP
             1,  10,    0,    0,    0,    # PPA
             1,  10,    0,    0,    0,    # PPK3_A
             1,  10,    0,    0,    0,    # PPK3_U
             0,  10,    0,    0,    0,    # PPK3_tetra
             1,  10,    0,    0,    0,    # UDK_ATP
             0,  10,    0,    0,    0,    # UDK_ADP
             0,  10,    0,    0,    0,    # UDK_ATPP
             1,  10,    0,    0,    0,    # UMPK_ATP
             0,  10,    0,    0,    0,    # UMPK_ADP
             0,  10,    0,    0,    0]    # UMPK_ATPP
start_var = pd.DataFrame(np.array(start_var).reshape(16, 5),
                         columns=['base', 'term1', 'term2', 'term3', 'term4'],
                         index=['ADP_Decay_v1', 'ADP_Decay_v2', 'GLMU', 'NAHK_ATP', 'NAHK_ADP', 'NAHK_ATPP', 'PPA', 'PPK3_A', 'PPK3_U', 'PPK3_tetra', 'UDK_ATP', 'UDK_ADP', 'UDK_ATPP', 'UMPK_ATP', 'UMPK_ADP', 'UMPK_ATPP'])
# define variable terms (those that will be successively added) where the first two elements are the row and column coordinates pointing to elements in the start variant data frame and the third element is the index of the base or regulation term; tuples are placed in inner lists because combinations of tuples are also valid elements of the vari_terms list
vari_terms = [[(0,0,1)],                                    # ADP_Decay_v1.on
              [(1,0,1)],                                    # ADP_Decay_v2.on

              [(2,2,2)],                                    # GLMU.ADP_Inhib
              [(2,3,16)],                                   # GLMU.UDP_GalNAc_Inhib

              [(3,0,2), (4,0,1), (5,0,0)],                  # change base kin. of NAHK_ATP to Conv. Kin. + comp. inhib. (#2) and activate NAHK_ADP with Conv. Kin. + comp. inhib. (#1)
              [(3,0,3), (4,0,0), (5,0,1)],                  # change base kin. of NAHK_ATP to Conv. Kin. + comp. inhib. (#3) and activate NAHK_ATPP with Conv. Kin. + comp. inhib. (#1)
              [(3,0,4), (4,0,2), (5,0,2)],                  # change base kin. of NAHK_ATP to Conv. Kin. + comp. inhib. (#4) and activate NAHK_ADP with Conv. Kin. + comp. inhib. (#2) and NAHK_ATPP with Conv. Kin. + comp. inhib. (#2)
              [(3,2,2), (4,2,2), (5,2,2)],                  # NAHK.ADP_Inhib (for all associated reactions: NAHK_ATP, NAHK_ADP, NAHK_ATPP)

              [(6,2,2)],                                    # PPA.ADP_Inhib
              [(6,3,5)],                                    # PPA.ATP_Act
              [(6,4,11)],                                   # PPA.PP_Act

              [(7,0,2), (8,0,2), (9,0,1)],                  # change base kin. of PPK3_A to Conv. Kin. + comp. inhib. (#2) and activate PPK3_U with Conv. Kin. + comp. inhib. (#2) and PPK3_tetra with Conv. Kin + comp. inhib (#1)
              [(7,2,2), (8,2,2), (9,2,2)],                  # PPK3.ADP_Inhib (for all associated reactions: PPK3_A, PPK3_U, and PPK3_tetra)
              [(7,3,12), (8,3,12), (9,3,12)],               # PPK3.PP_Inhib (for all associated reactions: PPK3_A, PPK3_U, and PPK3_tetra)

              [(10,0,2), (11,0,1), (12,0,0)],               # change base kin. of UDK_ATP to Conv. Kin. + comp. inhib. (#2) and activate UDK_ADP with Conv. Kin. + comp. inhib. (#1)
              [(10,0,3), (11,0,0), (12,0,1)],               # change base kin. of UDK_ATP to Conv. Kin. + comp. inhib. (#3) and activate UDK_ATPP with Conv. Kin. + comp. inhib. (#1)
              [(10,0,4), (11,0,2), (12,0,2)],               # change base kin. of UDK_ATP to Conv. Kin. + comp. inhib. (#4) and activate UDK_ADP with Conv. Kin. + comp. inhib. (#2) and UDK_ATPP with Conv. Kin. + comp. inhib. (#2)
              [(10,2,2), (11,2,2), (12,2,2)],               # UDK.ADP_Inhib (for all associated reactions: UDK_ATP, UDK_ADP, UDK_ATPP)
              [(10,3,22), (11,3,22), (12,3,22)],            # UDK.UTP_Inhib (for all associated reactions: UDK_ATP, UDK_ADP, UDK_ATPP)

              [(13,0,2), (14,0,1), (15,0,0)],               # change base kin. of UMPK_ATP to Conv. Kin. + comp. inhib. (#2) and activate UMPK_ADP with Conv. Kin. + comp. inhib. (#1)
              [(13,0,3), (14,0,0), (15,0,1)],               # change base kin. of UMPK_ATP to Conv. Kin. + comp. inhib. (#3) and activate UMPK_ATPP with Conv. Kin. + comp. inhib. (#1)
              [(13,0,4), (14,0,2), (15,0,2)],               # change base kin. of UMPK_ATP to Conv. Kin. + comp. inhib. (#4) and activate UMPK_ADP with Conv. Kin. + comp. inhib. (#2) and UMPK_ATPP with Conv. Kin. + comp. inhib. (#2)
              [(13,2,2) , (14,2,2), (15,2,2)]]              # UMPK.ADP_Inhib (for all associated reactions: UMPK_ATP, UMPK_ADP, UMPK_ATPP)
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
# define global settings for all parameter estimations
PE_algorithm_info={'name': 'Particle Swarm', 
                   'settings': {'method': {'Iteration Limit': 800,
                                           'Swarm Size': 50,
                                           'Std. Deviation': 1e-06}}}
# define which model selection criteria are to be used
selected_MSC = ['AIC', 'AICc', 'BIC', 'CIC1', 'CIC2', 'CIC3']
# use the existing parameter ensemble of the initial core model (negative control; here the start variant) instead of evaluating the start variant inside the improved_extension_search function
reading_file = open('v23brep4NegControl_Exp36.37.39_50runs_evaluated_vars_log', 'rb')
neg_ctrl_log = pickle.load(reading_file)
reading_file.close()
start_var_ensemble = neg_ctrl_log

# run improved extension search and store result (output dictionary that contains list of constructed 
# model variants and analysis result of the final model selection step)
start = time.perf_counter()
output_dict = improved_extension_search(start_var, vari_terms,
                                        exp_data_file_names, exp_data_dataframes,
                                        PE_algorithm_info=PE_algorithm_info,
                                        selected_MSC=selected_MSC,
                                        start_var_ensemble=start_var_ensemble)
end = time.perf_counter()
elapsed_hours = (end-start)/3600 # convert seconds to hours
print(f'Time taken: {elapsed_hours:.2f} hours')

# store result as pickle file
storage_file = open('SSO_ImpExtSearch_v23b_rep4_output_dict', 'wb')
pickle.dump(output_dict, storage_file)
storage_file.close()

# remove 'dump' text files that are created by the 'add_experiment' function that is called by the 
# eval_struct_var function; windows and unix systems use different file path separators for name in 
# exp_data_names
clean_up_exp_dump_files(exp_data_file_names)

# ------------------------------------------------------------------------------------------------------ #