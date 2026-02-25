#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

sys.path.append("..")
from batch_opt_func_lib import *

""" Analysis of the influence of different sets of kinetic model parameters (taken from an ensemble of parameter set) on the results of optimization problems (i.e., the predicted optimal initial enzyme and substrate concentrations). The optimization problems are defined in the 'opt_prob_defs.pdf' document (definitions are based on templates stored in the 'Collection_of_Optimization_Problems' directory).

    Approach: Define the optimization problems and solve them for each of the estimated parameter sets from the parameter ensemble (100 times). Plot the optimization results (= predicted initial enzyme and substrate concentrations) as combination of box and strip plot.

    Main package: basiCO - simplified Copasi Python API
                  <https://github.com/copasi/basico
"""

# SETUP
# load model structure and parameter ensemble
model_path = 'ImpExtSearch_v23b_rep2_MV5.cps'
fits_data_path = 'sampling_output_Particle_Swarm50runs_ImpExtSearch_v23b_rep2_ModelVar5.csv'

# OP05withS
""" Description
    -----------
    Target function: min E_tot
    Opt. variables: all initial enzyme and substrate concentrations
    Constraints: [UDP-GalNAc at 24h] >= 40.37 mM
                 Uri_Yield >= 82.74%
    Method: Genetic Algorithm (default settings)
"""
# solve optimization problem
result_file_name = "ImpExtSearch_v23b_rep2_MV5_PartSwarm50x_OP05withS.pkl"
opt_results = solve_OP05withS(model_path, fits_data_path, result_file_name)
