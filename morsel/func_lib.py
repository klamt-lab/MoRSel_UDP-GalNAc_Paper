#!/usr/bin/env python
# -*- coding: utf-8 -*-

from basico import *
from collections import OrderedDict
import itertools
from joblib import Parallel, delayed 
import math
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path
import pickle
import re
import seaborn as sns
from tqdm import tqdm


""" Collection of functions used by the super structure optimization calculation, analysis and visualization 
scripts as well as within other functions of this library.
"""


# ------------------------------------------------------------------------------------------------------ #

# MODEL SPECIFIC FUNCTIONS

# 1) <model_name>_data function: define reaction schemes (network stoichiometry), literature values for kinetic parameters of the initial core model and initial species concentrations (as placeholder - the parameter estimations take the initial concentrations from the first rows of the provided experimental data files at t=0)

# 2) <model_name>_term_libs function: define the model-specific 'building blocks' (variable terms) for the super structure optimization (SSO): each rate law can be a combination of a base kinetics term (multiple reaction-specific options) and one or more regulation terms (reaction-independent)


# UDP-GalNAc Model

def UDP_GalNAc_model_data():
    """ Define (1) string representations (schemes) of all reactions of the UDP-GalNAc model according to the syntax used by COPASI, (2) a dictionary of parameter values, and (3) a dictionary of initial concentrations for all species. Model units are hours, mmol, and liters.
    
    :return reaction_scheme_dict:
    :type reaction_scheme_dict:
    :return param_dict:
    :type param_dict:
    :return init_conc_dict:
    :type init_conc_dict:
    """

    # define reaction schemes; all of them are set even if a reaction is disabled in the structural variant matrix (in this case the underlying rate law is set to 0); the schemes are used as default and are modified according to the variable terms set in the structural variant matrix; modifiers that are not used in the selected rate law are inactive
    ## v23: add P to products of ADP_Decay_v1
    reaction_scheme_dict = {'ADP_Decay_v1': 'ADP -> AMP + P',
                            'ADP_Decay_v2': 'ADP + ADP -> ATP + AMP',
                            'GLMU': 'GalNAc1P + UTP = UDP_GalNAc + PP;  E_GLMU',
                            'NAHK_ATP': 'GalNAc + ATP = GalNAc1P + ADP;  E_NAHK ATPP AMP',
                            'NAHK_ADP': 'GalNAc + ADP = GalNAc1P + AMP;  E_NAHK ATPP ATP',
                            'NAHK_ATPP': 'GalNAc + ATPP = GalNAc1P + ATP;  E_NAHK ADP AMP',
                            'PPA': 'PP -> 2 * P;  E_PPA',
                            'PPK3_A': 'ADP + PolyP = ATP;  E_PPK3 ATPP UTP UDP',
                            'PPK3_U': 'UDP + PolyP = UTP;  E_PPK3 ATPP ATP ADP',
                            'PPK3_tetra': 'ATP + PolyP = ATPP;  E_PPK3 ADP UTP UDP',
                            'UDK_ATP': 'Uri + ATP = UMP + ADP;  E_UDK ATPP AMP',
                            'UDK_ADP': 'Uri + ADP = UMP + AMP;  E_UDK ATPP ATP',
                            'UDK_ATPP': 'Uri + ATPP = UMP + ATP;  E_UDK ADP AMP',
                            'UMPK_ATP': 'UMP + ATP = UDP + ADP;  E_UMPK ATPP AMP',
                            'UMPK_ADP': 'UMP + ADP = UDP + AMP;  E_UMPK ATPP ATP',
                            'UMPK_ATPP': 'UMP + ATPP = UDP + ATP;  E_UMPK ADP AMP'}

    # collect known literature values of kinetic parameters; not all parameters have to appear here - some may not be used depending on the kinetics that are selected according to the structural variant matrix; kcat [1/h]; Km, ki, ka [mM]; K_eq [-]
    ## v20c: remove all lit. data except known Keq values of the main _ATP variants
    ## v23b: add NAHK literature information (as no experimental data points are available for NAHK substrates and products so we need some information for the fit)
    param_dict = {'ADP_Decay_v1': {},
                  'ADP_Decay_v2': {},
                  'GLMU': {'K_eq_GLMU': 210},               # [eQuilibrator; pH 9; pMg=-log10(0.120)=0.9]
                  'NAHK_ATP': {'K_eq_NAHK': 60,             # [eQuilibrator; pH 9; pMg=-log10(0.120)=0.9]
                               'Km_GalNAc': 0.065,          # [Nishimoto2007]
                               'Km_ATP': 0.172,             # [Nishimoto2007]
                               'kcat_F': 2707.2},           # [Nishimoto2007]
                  'NAHK_ADP': {'K_eq_NAHK': 60,             # [eQuilibrator; pH 9; pMg=-log10(0.120)=0.9]
                               'Km_GalNAc': 0.065,          # [Nishimoto2007]
                               'kcat_F': 2707.2},           # [Nishimoto2007],
                  'NAHK_ATPP': {'K_eq_NAHK': 60,            # [eQuilibrator; pH 9; pMg=-log10(0.120)=0.9]
                                'Km_GalNAc': 0.065,         # [Nishimoto2007]
                                'kcat_F': 2707.2},          # [Nishimoto2007],
                  'PPA': {},
                  'PPK3_A': {},
                  'PPK3_U': {},
                  'PPK3_tetra': {},
                  'UDK_ATP': {'K_eq_UDK': 12000},           # [eQuilibrator; pH 9; pMg=-log10(0.120)=0.9]
                  'UDK_ADP': {},
                  'UDK_ATPP': {},
                  'UMPK_ATP': {'K_eq_UMPK': 5},             # [eQuilibrator; pH 9; pMg=-log10(0.120)=0.9]
                  'UMPK_ADP': {},
                  'UMPK_ATPP': {}
    }

    # define all non-zero initial species concentrations (the initial concentrations for the remaining species will be set to a very low value (1e-12) - setting them to zero can make problems with certain rate laws); initial concentrations come from Tuan:Exp36_5 and are just used as a placeholder when creating the COPASI model object - during the estimation loop basico changes the initial species concentrations to those found in the provided experimental data files
    init_conc_dict = {'GalNAc': 5.26361053933333,
                      'ATP': 0.509636666666667, # ATP+ATPP+ATPPP at t0
                      'ADP': 0.007756666666667,
                      'Uri': 5.26361053933333,
                      'E_PPA': 0.001553358,
                      'E_PPK3': 0.002878526,
                      'E_UDK': 0.00410627,
                      'E_NAHK': 0.002506077,
                      'E_GLMU': 0.002032934,
                      'E_UMPK': 0.004463289,
                      'PolyP': 256} # 32 mM PolyP * 8 available P bonds per average length chain (n=14)

    return reaction_scheme_dict, param_dict, init_conc_dict


def UDP_GalNAc_model_term_libs():
    """ Define the libraries of model specific base and regulation terms. Base terms are stored in a dictionary where the keys are reaction names and the values are lists of dictionaries; each of them contains a kinetic rate law equation and the associated mapping of all variables and parameters. Regulation terms are stored in a list of dictionaries which contain names, equations, and associated mappings.
    
    :return base_kinetics: dictionary of base kinetics (equations and mappings) per reaction
    :type base_kinetics: dictionary
    :return regulatory_terms: list of regulation terms (each entry is a dictionary with the name, equation, and mapping of the regulation term)
    :type regulatory_terms: list
    """

# create repositories of reaction specific base kinetics and general regulatory terms that can be indexed according to the integers contained in the structural variant vector (0: no base kinetics <-> reaction turned off, 1: first base kinetics in sublist, and so on); all kinetics come with dictionaries mapping variable names to their usage, i.e., 'parameter', 'substrate', 'product', or 'modifier'
# regarding the base kinetics: some reactions are 'paired', i.e., they are catalyzed by the same enzyme (and some substrates/products are the same); there are extra base kinetics with built-in competitive inhibition terms for combinations of these paired reactions (competitive inhibition terms are equations that scale the Km values of all substrates and products that compete for binding; the implementation is based on the approach described in section 3.1  of the supplementary material of Schäuble2013, <https://doi.org/10.1016/j.febslet.2013.06.025>); possible combinations of base kinetics for ([NAHK/UDK/UMPK]_ATP, [NAHK/UDK/UMPK]_ADP, [NAHK/UDK/UMPK]_ATPP): either (1,0,0), (2,1,0), (3,0,1), or (4,2,2); for (PPK3_A,PPK3_U,PPK3_tetra): either (1,1,0) or (2,2,2)
    base_kinetics = {
        'ADP_Decay_v1': [{'equation': '0', 'mapping': {}}, #0
                         {'equation': '(k_MA*ADP)', 'mapping': {'k_MA': 'parameter', 'ADP': 'substrate'}}], #1
        'ADP_Decay_v2': [{'equation': '0', 'mapping': {}}, #0
                         {'equation': '(k_MA*(ADP)^2)', 'mapping': {'k_MA': 'parameter', 'ADP': 'substrate'}}], #1
        'GLMU': [{'equation': '0', 'mapping': {}}, #0
                 {'equation': '(E_GLMU*(((kcat_F*(GalNAc1P/Km_GalNAc1P)*(UTP/Km_UTP))-(kcat_F*(1/K_eq_GLMU)*((UDP_GalNAc*PP)/(Km_GalNAc1P*Km_UTP))))/(((1+(GalNAc1P/Km_GalNAc1P))*(1+(UTP/Km_UTP)))+((1+(UDP_GalNAc/Km_UDP_GalNAc))*(1+(PP/Km_PP)))-1)))', 'mapping': {'E_GLMU': 'modifier', 'kcat_F': 'parameter', 'UTP': 'substrate', 'Km_UTP': 'parameter', 'PP': 'product', 'Km_PP': 'parameter', 'GalNAc1P': 'substrate', 'Km_GalNAc1P': 'parameter', 'UDP_GalNAc': 'product', 'Km_UDP_GalNAc': 'parameter', 'K_eq_GLMU': 'parameter'}}], #1
        'NAHK_ATP': [{'equation': '0', 'mapping': {}}, #0
                     {'equation': '(E_NAHK*(((kcat_F*(GalNAc/Km_GalNAc)*(ATP/Km_ATP))-(kcat_F*(1/K_eq_NAHK)*((GalNAc1P*ADP)/(Km_GalNAc*Km_ATP))))/(((1+(GalNAc/Km_GalNAc))*(1+(ATP/Km_ATP)))+(((1+(GalNAc1P/Km_GalNAc1P))*(1+(ADP/Km_ADP)))-1))))', 'mapping': {'E_NAHK': 'modifier', 'kcat_F': 'parameter', 'ATP': 'substrate', 'Km_ATP': 'parameter', 'ADP': 'product', 'Km_ADP': 'parameter', 'GalNAc': 'substrate', 'Km_GalNAc': 'parameter', 'GalNAc1P': 'product', 'Km_GalNAc1P': 'parameter', 'K_eq_NAHK': 'parameter'}}, #1 Convenience kinetics
                     {'equation': '(E_NAHK*(((kcat_F*(GalNAc/Km_GalNAc)*(ATP/(Km_ATP*(1+(ADP/Km_ADP)+(AMP/Km_AMP)))))-(kcat_F*(1/K_eq_NAHK)*((GalNAc1P*ADP)/(Km_GalNAc*(Km_ATP*(1+(ADP/Km_ADP)+(AMP/Km_AMP)))))))/(((1+(GalNAc/Km_GalNAc))*(1+(ATP/(Km_ATP*(1+(ADP/Km_ADP)+(AMP/Km_AMP))))))+(((1+(GalNAc1P/Km_GalNAc1P))*(1+(ADP/(Km_ADP*(1+(ADP/Km_ADP)+(AMP/Km_AMP))))))-1))))', 'mapping': {'E_NAHK': 'modifier', 'kcat_F': 'parameter', 'ATP': 'substrate', 'Km_ATP': 'parameter', 'ADP': 'product', 'Km_ADP': 'parameter', 'GalNAc': 'substrate', 'Km_GalNAc': 'parameter', 'GalNAc1P': 'product', 'Km_GalNAc1P': 'parameter', 'K_eq_NAHK': 'parameter', 'AMP': 'modifier', 'Km_AMP': 'parameter'}}, #2 Convenience kinetics with comp. inhib. for NAHK_ADP
                     {'equation': '(E_NAHK*(((kcat_F*(GalNAc/Km_GalNAc)*(ATP/(Km_ATP*(1+(ATPP/Km_ATPP)+(ATP/Km_ATP)))))-(kcat_F*(1/K_eq_NAHK)*((GalNAc1P*ADP)/(Km_GalNAc*(Km_ATP*(1+(ATPP/Km_ATPP)+(ATP/Km_ATP)))))))/(((1+(GalNAc/Km_GalNAc))*(1+(ATP/(Km_ATP*(1+(ATPP/Km_ATPP)+(ATP/Km_ATP))))))+(((1+(GalNAc1P/Km_GalNAc1P))*(1+(ADP/(Km_ADP*(1+(ATPP/Km_ATPP)+(ATP/Km_ATP))))))-1))))', 'mapping': {'E_NAHK': 'modifier', 'kcat_F': 'parameter', 'ATP': 'substrate', 'Km_ATP': 'parameter', 'ADP': 'product', 'Km_ADP': 'parameter', 'GalNAc': 'substrate', 'Km_GalNAc': 'parameter', 'GalNAc1P': 'product', 'Km_GalNAc1P': 'parameter', 'K_eq_NAHK': 'parameter', 'ATPP': 'modifier', 'Km_ATPP': 'parameter'}}, #3 Convenience kinetics with comp. inhib. for NAHK_ATPP
                     {'equation': '(E_NAHK*(((kcat_F*(GalNAc/Km_GalNAc)*(ATP/(Km_ATP*(1+(ADP/Km_ADP)+(ATPP/Km_ATPP)+(AMP/Km_AMP)+(ATP/Km_ATP)))))-(kcat_F*(1/K_eq_NAHK)*((GalNAc1P*ADP)/(Km_GalNAc*(Km_ATP*(1+(ADP/Km_ADP)+(ATPP/Km_ATPP)+(AMP/Km_AMP)+(ATP/Km_ATP)))))))/(((1+(GalNAc/Km_GalNAc))*(1+(ATP/(Km_ATP*(1+(ADP/Km_ADP)+(ATPP/Km_ATPP)+(AMP/Km_AMP)+(ATP/Km_ATP))))))+(((1+(GalNAc1P/Km_GalNAc1P))*(1+(ADP/(Km_ADP*(1+(ADP/Km_ADP)+(ATPP/Km_ATPP)+(AMP/Km_AMP)+(ATP/Km_ATP))))))-1))))', 'mapping': {'E_NAHK': 'modifier', 'kcat_F': 'parameter', 'ATP': 'substrate', 'Km_ATP': 'parameter', 'ADP': 'product', 'Km_ADP': 'parameter', 'GalNAc': 'substrate', 'Km_GalNAc': 'parameter', 'GalNAc1P': 'product', 'Km_GalNAc1P': 'parameter', 'K_eq_NAHK': 'parameter', 'ATPP': 'modifier', 'Km_ATPP': 'parameter', 'AMP': 'modifier', 'Km_AMP': 'parameter'}}], #4 Convenience kinetics with comp. inhib. for NAHK_ADP and NAHK_ATPP
        'NAHK_ADP': [{'equation': '0', 'mapping': {}}, #0
                     {'equation': '(E_NAHK*(((kcat_F*(GalNAc/Km_GalNAc)*(ADP/(Km_ADP*(1+(ATP/Km_ATP)+(ADP/Km_ADP)))))-(kcat_F*(1/K_eq_NAHK)*((GalNAc1P*AMP)/(Km_GalNAc*(Km_ADP*(1+(ATP/Km_ATP)+(ADP/Km_ADP)))))))/(((1+(GalNAc/Km_GalNAc))*(1+(ADP/(Km_ADP*(1+(ATP/Km_ATP)+(ADP/Km_ADP))))))+(((1+(GalNAc1P/Km_GalNAc1P))*(1+(AMP/(Km_AMP*(1+(ATP/Km_ATP)+(ADP/Km_ADP))))))-1))))', 'mapping': {'E_NAHK': 'modifier', 'kcat_F': 'parameter', 'ADP': 'substrate', 'Km_ADP': 'parameter', 'AMP': 'product', 'Km_AMP': 'parameter', 'GalNAc': 'substrate', 'Km_GalNAc': 'parameter', 'GalNAc1P': 'product', 'Km_GalNAc1P': 'parameter', 'K_eq_NAHK': 'parameter', 'ATP': 'modifier', 'Km_ATP': 'parameter'}}, #1 Convenience kinetics with comp. inhib. for NAHK_ATP
                     {'equation': '(E_NAHK*(((kcat_F*(GalNAc/Km_GalNAc)*(ADP/(Km_ADP*(1+2*(ATP/Km_ATP)+(ATPP/Km_ATPP)+(ADP/Km_ADP)))))-(kcat_F*(1/K_eq_NAHK)*((GalNAc1P*AMP)/(Km_GalNAc*(Km_ADP*(1+2*(ATP/Km_ATP)+(ATPP/Km_ATPP)+(ADP/Km_ADP)))))))/(((1+(GalNAc/Km_GalNAc))*(1+(ADP/(Km_ADP*(1+2*(ATP/Km_ATP)+(ATPP/Km_ATPP)+(ADP/Km_ADP))))))+(((1+(GalNAc1P/Km_GalNAc1P))*(1+(AMP/(Km_AMP*(1+2*(ATP/Km_ATP)+(ATPP/Km_ATPP)+(ADP/Km_ADP))))))-1))))', 'mapping': {'E_NAHK': 'modifier', 'kcat_F': 'parameter', 'ADP': 'substrate', 'Km_ADP': 'parameter', 'AMP': 'product', 'Km_AMP': 'parameter', 'GalNAc': 'substrate', 'Km_GalNAc': 'parameter', 'GalNAc1P': 'product', 'Km_GalNAc1P': 'parameter', 'K_eq_NAHK': 'parameter', 'ATP': 'modifier', 'Km_ATP': 'parameter', 'ATPP': 'modifier', 'Km_ATPP': 'parameter'}}], #2 Convenience kinetics with comp. inhib. for NAHK_ATP and NAHK_ATPP
        'NAHK_ATPP': [{'equation': '0', 'mapping': {}}, #0
                      {'equation': '(E_NAHK*(((kcat_F*(GalNAc/Km_GalNAc)*(ATPP/(Km_ATPP*(1+(ATP/Km_ATP)+(ADP/Km_ADP)))))-(kcat_F*(1/K_eq_NAHK)*((GalNAc1P*ATP)/(Km_GalNAc*(Km_ATPP*(1+(ATP/Km_ATP)+(ADP/Km_ADP)))))))/(((1+(GalNAc/Km_GalNAc))*(1+(ATPP/(Km_ATPP*(1+(ATP/Km_ATP)+(ADP/Km_ADP))))))+(((1+(GalNAc1P/Km_GalNAc1P))*(1+(ATP/(Km_ATP*(1+(ATP/Km_ATP)+(ADP/Km_ADP))))))-1))))', 'mapping': {'E_NAHK': 'modifier', 'kcat_F': 'parameter', 'ATPP': 'substrate', 'Km_ATPP': 'parameter', 'ATP': 'product', 'Km_ATP': 'parameter', 'GalNAc': 'substrate', 'Km_GalNAc': 'parameter', 'GalNAc1P': 'product', 'Km_GalNAc1P': 'parameter', 'K_eq_NAHK': 'parameter', 'ADP': 'modifier', 'Km_ADP': 'parameter'}}, #1 Convenience kinetics with comp. inhibition for NAHK_ATP
                      {'equation': '(E_NAHK*(((kcat_F*(GalNAc/Km_GalNAc)*(ATPP/(Km_ATPP*(1+(ATP/Km_ATP)+2*(ADP/Km_ADP)+(AMP/Km_AMP)))))-(kcat_F*(1/K_eq_NAHK)*((GalNAc1P*ATP)/(Km_GalNAc*(Km_ATPP*(1+(ATP/Km_ATP)+2*(ADP/Km_ADP)+(AMP/Km_AMP)))))))/(((1+(GalNAc/Km_GalNAc))*(1+(ATPP/(Km_ATPP*(1+(ATP/Km_ATP)+2*(ADP/Km_ADP)+(AMP/Km_AMP))))))+(((1+(GalNAc1P/Km_GalNAc1P))*(1+(ATP/(Km_ATP*(1+(ATP/Km_ATP)+2*(ADP/Km_ADP)+(AMP/Km_AMP))))))-1))))', 'mapping': {'E_NAHK': 'modifier', 'kcat_F': 'parameter', 'ATPP': 'substrate', 'Km_ATPP': 'parameter', 'ATP': 'product', 'Km_ATP': 'parameter', 'GalNAc': 'substrate', 'Km_GalNAc': 'parameter', 'GalNAc1P': 'product', 'Km_GalNAc1P': 'parameter', 'K_eq_NAHK': 'parameter', 'ADP': 'modifier', 'Km_ADP': 'parameter', 'AMP': 'modifier', 'Km_AMP': 'parameter'}}], #2 Convenience kinetics with comp. inhib. for NAHK_ATP and NAHK_ADP
        'PPA': [{'equation': '0', 'mapping': {}}, #0
                {'equation': '(E_PPA*kcat_F*(1/(1+(Km_PP/PP))))', 'mapping': {'E_PPA': 'modifier', 'kcat_F': 'parameter', 'PP': 'substrate', 'Km_PP': 'parameter'}}], #1
        'PPK3_A': [{'equation': '0', 'mapping': {}}, #0
                   {'equation': '(E_PPK3*(((kcat_F*(ADP/(Km_ADP*(1+(UDP/Km_UDP)+(UTP/Km_UTP))))*(PolyP/Km_PolyP))-(kcat_F*(1/K_eq_PPK3_A)*((ATP)/((Km_ADP*(1+(UDP/Km_UDP)+(UTP/Km_UTP)))*Km_PolyP))))/(((1+(ADP/(Km_ADP*(1+(UDP/Km_UDP)+(UTP/Km_UTP)))))*(1+(PolyP/Km_PolyP)))+(1+(ATP/(Km_ATP*(1+(UDP/Km_UDP)+(UTP/Km_UTP)))))-1)))', 'mapping': {'E_PPK3': 'modifier', 'kcat_F': 'parameter', 'ADP': 'substrate', 'Km_ADP': 'parameter', 'ATP': 'product', 'Km_ATP': 'parameter', 'K_eq_PPK3_A': 'parameter', 'UDP': 'modifier', 'Km_UDP': 'parameter', 'UTP': 'modifier', 'Km_UTP': 'parameter', 'PolyP': 'substrate', 'Km_PolyP': 'parameter'}}, #1 Convenience kinetics with comp. inhib. for PPK3_U
                   {'equation': '(E_PPK3*(((kcat_F*(ADP/(Km_ADP*(1+(UDP/Km_UDP)+(ATP/Km_ATP)+(UTP/Km_UTP)+(ATPP/Km_ATPP))))*(PolyP/Km_PolyP))-(kcat_F*(1/K_eq_PPK3_A)*((ATP)/((Km_ADP*(1+(UDP/Km_UDP)+(ATP/Km_ATP)+(UTP/Km_UTP)+(ATPP/Km_ATPP)))*Km_PolyP))))/(((1+(ADP/(Km_ADP*(1+(UDP/Km_UDP)+(ATP/Km_ATP)+(UTP/Km_UTP)+(ATPP/Km_ATPP)))))*(1+(PolyP/Km_PolyP)))+(1+(ATP/(Km_ATP*(1+(UDP/Km_UDP)+(ATP/Km_ATP)+(UTP/Km_UTP)+(ATPP/Km_ATPP)))))-1)))', 'mapping': {'E_PPK3': 'modifier', 'kcat_F': 'parameter', 'ADP': 'substrate', 'Km_ADP': 'parameter', 'ATP': 'product', 'Km_ATP': 'parameter', 'K_eq_PPK3_A': 'parameter', 'UDP': 'modifier', 'Km_UDP': 'parameter', 'UTP': 'modifier', 'Km_UTP': 'parameter', 'PolyP': 'substrate', 'Km_PolyP': 'parameter', 'ATPP': 'modifier', 'Km_ATPP': 'parameter'}}], #2 Convenience kinetics with comp. inhib. for PPK3_U and PPK3_tetra
        'PPK3_U': [{'equation': '0', 'mapping': {}}, #0
                   {'equation': '(E_PPK3*(((kcat_F*(UDP/(Km_UDP*(1+(ADP/Km_ADP)+(ATP/Km_ATP))))*(PolyP/Km_PolyP))-(kcat_F*(1/K_eq_PPK3_U)*((UTP)/((Km_UDP*(1+(ADP/Km_ADP)+(ATP/Km_ATP)))*Km_PolyP))))/(((1+(UDP/(Km_UDP*(1+(ADP/Km_ADP)+(ATP/Km_ATP)))))*(1+(PolyP/Km_PolyP)))+(1+(UTP/(Km_UTP*(1+(ADP/Km_ADP)+(ATP/Km_ATP)))))-1)))', 'mapping': {'E_PPK3': 'modifier', 'kcat_F': 'parameter', 'UDP': 'substrate', 'Km_UDP': 'parameter', 'UTP': 'product', 'Km_UTP': 'parameter', 'K_eq_PPK3_U': 'parameter', 'ADP': 'modifier', 'Km_ADP': 'parameter', 'ATP': 'modifier', 'Km_ATP': 'parameter', 'PolyP': 'substrate', 'Km_PolyP': 'parameter'}}, #1 Convenience kinetics with comp. inhib. for PPK3_A
                   {'equation': '(E_PPK3*(((kcat_F*(UDP/(Km_UDP*(1+(ADP/Km_ADP)+2*(ATP/Km_ATP)+(ATPP/Km_ATPP))))*(PolyP/Km_PolyP))-(kcat_F*(1/K_eq_PPK3_U)*((UTP)/((Km_UDP*(1+(ADP/Km_ADP)+2*(ATP/Km_ATP)+(ATPP/Km_ATPP)))*Km_PolyP))))/(((1+(UDP/(Km_UDP*(1+(ADP/Km_ADP)+2*(ATP/Km_ATP)+(ATPP/Km_ATPP)))))*(1+(PolyP/Km_PolyP)))+(1+(UTP/(Km_UTP*(1+(ADP/Km_ADP)+2*(ATP/Km_ATP)+(ATPP/Km_ATPP)))))-1)))', 'mapping': {'E_PPK3': 'modifier', 'kcat_F': 'parameter', 'UDP': 'substrate', 'Km_UDP': 'parameter', 'UTP': 'product', 'Km_UTP': 'parameter', 'K_eq_PPK3_U': 'parameter', 'ADP': 'modifier', 'Km_ADP': 'parameter', 'ATP': 'modifier', 'Km_ATP': 'parameter', 'PolyP': 'substrate', 'Km_PolyP': 'parameter', 'ATPP': 'modifier', 'Km_ATPP': 'parameter'}}], #2 Convenience kinetics with comp. inhib. for PPK3_U and PPK3_tetra
        'PPK3_tetra': [{'equation': '0', 'mapping': {}}, #0
                       {'equation': '(E_PPK3*(((kcat_F*(ATP/(Km_ATP*(1+(ADP/Km_ADP)+(UDP/Km_UDP)+(ATP/Km_ATP)+(UTP/Km_UTP))))*(PolyP/Km_PolyP))-(kcat_F*(1/K_eq_PPK3_tetra)*((ATPP)/((Km_ATP*(1+(ADP/Km_ADP)+(UDP/Km_UDP)+(ATP/Km_ATP)+(UTP/Km_UTP)))*Km_PolyP))))/(((1+(ATP/(Km_ATP*(1+(ADP/Km_ADP)+(UDP/Km_UDP)+(ATP/Km_ATP)+(UTP/Km_UTP)))))*(1+(PolyP/Km_PolyP)))+(1+(ATPP/(Km_ATPP*(1+(ADP/Km_ADP)+(UDP/Km_UDP)+(ATP/Km_ATP)+(UTP/Km_UTP)))))-1)))', 'mapping': {'E_PPK3': 'modifier', 'kcat_F': 'parameter', 'ATP': 'substrate', 'Km_ATP': 'parameter', 'ATPP': 'product', 'Km_ATPP': 'parameter', 'K_eq_PPK3_tetra': 'parameter', 'UDP': 'modifier', 'Km_UDP': 'parameter', 'UTP': 'modifier', 'Km_UTP': 'parameter', 'PolyP': 'substrate', 'Km_PolyP': 'parameter', 'ADP': 'modifier', 'Km_ADP': 'parameter'}}], #1
        'UDK_ATP': [{'equation': '0', 'mapping': {}}, #0
                    {'equation': '(E_UDK*(((kcat_F*(Uri/Km_Uri)*(ATP/Km_ATP))-(kcat_F*(1/K_eq_UDK)*((UMP*ADP)/(Km_Uri*Km_ATP))))/(((1+(Uri/Km_Uri))*(1+(ATP/Km_ATP)))+((1+(UMP/Km_UMP))*(1+(ADP/Km_ADP)))-1)))', 'mapping': {'E_UDK': 'modifier', 'kcat_F': 'parameter', 'ATP': 'substrate', 'Km_ATP': 'parameter', 'ADP': 'product', 'Km_ADP': 'parameter', 'Uri': 'substrate', 'Km_Uri': 'parameter', 'UMP': 'product', 'Km_UMP': 'parameter', 'K_eq_UDK': 'parameter'}}, #1 Convenience kinetics
                    {'equation': '(E_UDK*(((kcat_F*(Uri/Km_Uri)*(ATP/(Km_ATP*(1+(ADP/Km_ADP)+(AMP/Km_AMP)))))-(kcat_F*(1/K_eq_UDK)*((UMP*ADP)/(Km_Uri*(Km_ATP*(1+(ADP/Km_ADP)+(AMP/Km_AMP)))))))/(((1+(Uri/Km_Uri))*(1+(ATP/(Km_ATP*(1+(ADP/Km_ADP)+(AMP/Km_AMP))))))+((1+(UMP/Km_UMP))*(1+(ADP/(Km_ADP*(1+(ADP/Km_ADP)+(AMP/Km_AMP))))))-1)))', 'mapping': {'E_UDK': 'modifier', 'kcat_F': 'parameter', 'ATP': 'substrate', 'Km_ATP': 'parameter', 'ADP': 'product', 'Km_ADP': 'parameter', 'Uri': 'substrate', 'Km_Uri': 'parameter', 'UMP': 'product', 'Km_UMP': 'parameter', 'K_eq_UDK': 'parameter', 'AMP': 'modifier', 'Km_AMP': 'parameter'}}, #2 Convenience kinetics with comp. inhib. for UDK_ADP
                    {'equation': '(E_UDK*(((kcat_F*(Uri/Km_Uri)*(ATP/(Km_ATP*(1+(ATPP/Km_ATPP)+(ATP/Km_ATP)))))-(kcat_F*(1/K_eq_UDK)*((UMP*ADP)/(Km_Uri*(Km_ATP*(1+(ATPP/Km_ATPP)+(ATP/Km_ATP)))))))/(((1+(Uri/Km_Uri))*(1+(ATP/(Km_ATP*(1+(ATPP/Km_ATPP)+(ATP/Km_ATP))))))+((1+(UMP/Km_UMP))*(1+(ADP/(Km_ADP*(1+(ATPP/Km_ATPP)+(ATP/Km_ATP))))))-1)))', 'mapping': {'E_UDK': 'modifier', 'kcat_F': 'parameter', 'ATP': 'substrate', 'Km_ATP': 'parameter', 'ADP': 'product', 'Km_ADP': 'parameter', 'Uri': 'substrate', 'Km_Uri': 'parameter', 'UMP': 'product', 'Km_UMP': 'parameter', 'K_eq_UDK': 'parameter', 'ATPP': 'modifier', 'Km_ATPP': 'parameter'}}, #3 Convenience kinetics with comp. inhib. for UDK_ATPP
                    {'equation': '(E_UDK*(((kcat_F*(Uri/Km_Uri)*(ATP/(Km_ATP*(1+(ADP/Km_ADP)+(ATPP/Km_ATPP)+(AMP/Km_AMP)+(ATP/Km_ATP)))))-(kcat_F*(1/K_eq_UDK)*((UMP*ADP)/(Km_Uri*(Km_ATP*(1+(ADP/Km_ADP)+(ATPP/Km_ATPP)+(AMP/Km_AMP)+(ATP/Km_ATP)))))))/(((1+(Uri/Km_Uri))*(1+(ATP/(Km_ATP*(1+(ADP/Km_ADP)+(ATPP/Km_ATPP)+(AMP/Km_AMP)+(ATP/Km_ATP))))))+((1+(UMP/Km_UMP))*(1+(ADP/(Km_ADP*(1+(ADP/Km_ADP)+(ATPP/Km_ATPP)+(AMP/Km_AMP)+(ATP/Km_ATP))))))-1)))', 'mapping': {'E_UDK': 'modifier', 'kcat_F': 'parameter', 'ATP': 'substrate', 'Km_ATP': 'parameter', 'ADP': 'product', 'Km_ADP': 'parameter', 'Uri': 'substrate', 'Km_Uri': 'parameter', 'UMP': 'product', 'Km_UMP': 'parameter', 'K_eq_UDK': 'parameter', 'ATPP': 'modifier', 'Km_ATPP': 'parameter', 'AMP': 'modifier', 'Km_AMP': 'parameter' }}], #4 Convenience kinetics with comp. inhib. for UDK_ADP and UDK_ATPP
        'UDK_ADP': [{'equation': '0', 'mapping': {}}, #0
                    {'equation': '(E_UDK*(((kcat_F*(Uri/Km_Uri)*(ADP/(Km_ADP*(1+(ATP/Km_ATP)+(ADP/Km_ADP)))))-(kcat_F*(1/K_eq_UDK)*((UMP*AMP)/(Km_Uri*(Km_ADP*(1+(ATP/Km_ATP)+(ADP/Km_ADP)))))))/(((1+(Uri/Km_Uri))*(1+(ADP/(Km_ADP*(1+(ATP/Km_ATP)+(ADP/Km_ADP))))))+((1+(UMP/Km_UMP))*(1+(AMP/(Km_AMP*(1+(ATP/Km_ATP)+(ADP/Km_ADP))))))-1)))', 'mapping': {'E_UDK': 'modifier', 'kcat_F': 'parameter', 'ADP': 'substrate', 'Km_ADP': 'parameter', 'AMP': 'product', 'Km_AMP': 'parameter', 'Uri': 'substrate', 'Km_Uri': 'parameter', 'UMP': 'product', 'Km_UMP': 'parameter', 'K_eq_UDK': 'parameter', 'ATP': 'modifier', 'Km_ATP': 'parameter'}}, #1 Convenience kinetics with comp. inhib for UDK_ATP
                    {'equation': '(E_UDK*(((kcat_F*(Uri/Km_Uri)*(ADP/(Km_ADP*(1+2*(ATP/Km_ATP)+(ATPP/Km_ATPP)+(ADP/Km_ADP)))))-(kcat_F*(1/K_eq_UDK)*((UMP*AMP)/(Km_Uri*(Km_ADP*(1+2*(ATP/Km_ATP)+(ATPP/Km_ATPP)+(ADP/Km_ADP)))))))/(((1+(Uri/Km_Uri))*(1+(ADP/(Km_ADP*(1+2*(ATP/Km_ATP)+(ATPP/Km_ATPP)+(ADP/Km_ADP))))))+((1+(UMP/Km_UMP))*(1+(AMP/(Km_AMP*(1+2*(ATP/Km_ATP)+(ATPP/Km_ATPP)+(ADP/Km_ADP))))))-1)))', 'mapping': {'E_UDK': 'modifier', 'kcat_F': 'parameter', 'ADP': 'substrate', 'Km_ADP': 'parameter', 'AMP': 'product', 'Km_AMP': 'parameter', 'Uri': 'substrate', 'Km_Uri': 'parameter', 'UMP': 'product', 'Km_UMP': 'parameter', 'K_eq_UDK': 'parameter', 'ATP': 'modifier', 'Km_ATP': 'parameter', 'ATPP': 'modifier', 'Km_ATPP': 'parameter'}}], #2 Convenience kinetics with comp. inhib. for UDK_ATP and UDK_ATPP
        'UDK_ATPP': [{'equation': '0', 'mapping': {}}, #0
                     {'equation': '(E_UDK*(((kcat_F*(Uri/Km_Uri)*(ATPP/(Km_ATPP*(1+(ATP/Km_ATP)+(ADP/Km_ADP)))))-(kcat_F*(1/K_eq_UDK)*((UMP*ATP)/(Km_Uri*(Km_ATPP*(1+(ATP/Km_ATP)+(ADP/Km_ADP)))))))/(((1+(Uri/Km_Uri))*(1+(ATPP/(Km_ATPP*(1+(ATP/Km_ATP)+(ADP/Km_ADP))))))+((1+(UMP/Km_UMP))*(1+(ATP/(Km_ATP*(1+(ATP/Km_ATP)+(ADP/Km_ADP))))))-1)))', 'mapping': {'E_UDK': 'modifier', 'kcat_F': 'parameter', 'ATPP': 'substrate', 'Km_ATPP': 'parameter', 'ATP': 'product', 'Km_ATP': 'parameter', 'Uri': 'substrate', 'Km_Uri': 'parameter', 'UMP': 'product', 'Km_UMP': 'parameter', 'K_eq_UDK': 'parameter', 'ADP': 'modifier', 'Km_ADP': 'parameter'}}, #1 Convenience kinetics with comp. inhib. for UDK_ATP
                     {'equation': '(E_UDK*(((kcat_F*(Uri/Km_Uri)*(ATPP/(Km_ATPP*(1+(ATP/Km_ATP)+2*(ADP/Km_ADP)+(AMP/Km_AMP)))))-(kcat_F*(1/K_eq_UDK)*((UMP*ATP)/(Km_Uri*(Km_ATPP*(1+(ATP/Km_ATP)+2*(ADP/Km_ADP)+(AMP/Km_AMP)))))))/(((1+(Uri/Km_Uri))*(1+(ATPP/(Km_ATPP*(1+(ATP/Km_ATP)+2*(ADP/Km_ADP)+(AMP/Km_AMP))))))+((1+(UMP/Km_UMP))*(1+(ATP/(Km_ATP*(1+(ATP/Km_ATP)+2*(ADP/Km_ADP)+(AMP/Km_AMP))))))-1)))', 'mapping': {'E_UDK': 'modifier', 'kcat_F': 'parameter', 'ATPP': 'substrate', 'Km_ATPP': 'parameter', 'ATP': 'product', 'Km_ATP': 'parameter', 'Uri': 'substrate', 'Km_Uri': 'parameter', 'UMP': 'product', 'Km_UMP': 'parameter', 'K_eq_UDK': 'parameter', 'ADP': 'modifier', 'Km_ADP': 'parameter', 'AMP': 'modifier', 'Km_AMP': 'parameter'}}], #2 Convenience kinetics with comp. inhib. for UDK_ATP and UDK_ADP
        'UMPK_ATP': [{'equation': '0', 'mapping': {}}, #0
                     {'equation': '(E_UMPK*(((kcat_F*(UMP/Km_UMP)*(ATP/Km_ATP))-(kcat_F*(1/K_eq_UMPK)*((UDP*ADP)/(Km_UMP*Km_ATP))))/(((1+(UMP/Km_UMP))*(1+(ATP/Km_ATP)))+((1+(UDP/Km_UDP))*(1+(ADP/Km_ADP)))-1)))', 'mapping': {'E_UMPK': 'modifier', 'kcat_F': 'parameter', 'UMP': 'substrate', 'Km_UMP': 'parameter', 'ATP': 'substrate', 'Km_ATP': 'parameter', 'UDP': 'product', 'Km_UDP': 'parameter', 'ADP': 'product', 'Km_ADP': 'parameter', 'K_eq_UMPK': 'parameter'}}, #1 Convenience kinetics
                     {'equation': '(E_UMPK*(((kcat_F*(UMP/Km_UMP)*(ATP/(Km_ATP*(1+(ADP/Km_ADP)+(AMP/Km_AMP)))))-(kcat_F*(1/K_eq_UMPK)*((UDP*ADP)/(Km_UMP*(Km_ATP*(1+(ADP/Km_ADP)+(AMP/Km_AMP)))))))/(((1+(UMP/Km_UMP))*(1+(ATP/(Km_ATP*(1+(ADP/Km_ADP)+(AMP/Km_AMP))))))+((1+(UDP/Km_UDP))*(1+(ADP/(Km_ADP*(1+(ADP/Km_ADP)+(AMP/Km_AMP))))))-1)))', 'mapping': {'E_UMPK': 'modifier', 'kcat_F': 'parameter', 'UMP': 'substrate', 'Km_UMP': 'parameter', 'ATP': 'substrate', 'Km_ATP': 'parameter', 'UDP': 'product', 'Km_UDP': 'parameter', 'ADP': 'product', 'Km_ADP': 'parameter', 'K_eq_UMPK': 'parameter', 'AMP': 'modifier', 'Km_AMP': 'parameter'}}, #2 Convenience kinetics with comp. inhib. for UMPK_ADP
                     {'equation': '(E_UMPK*(((kcat_F*(UMP/Km_UMP)*(ATP/(Km_ATP*(1+(ATPP/Km_ATPP)+(ATP/Km_ATP)))))-(kcat_F*(1/K_eq_UMPK)*((UDP*ADP)/(Km_UMP*(Km_ATP*(1+(ATPP/Km_ATPP)+(ATP/Km_ATP)))))))/(((1+(UMP/Km_UMP))*(1+(ATP/(Km_ATP*(1+(ATPP/Km_ATPP)+(ATP/Km_ATP))))))+((1+(UDP/Km_UDP))*(1+(ADP/(Km_ADP*(1+(ATPP/Km_ATPP)+(ATP/Km_ATP))))))-1)))', 'mapping': {'E_UMPK': 'modifier', 'kcat_F': 'parameter', 'UMP': 'substrate', 'Km_UMP': 'parameter', 'ATP': 'substrate', 'Km_ATP': 'parameter', 'UDP': 'product', 'Km_UDP': 'parameter', 'ADP': 'product', 'Km_ADP': 'parameter', 'K_eq_UMPK': 'parameter', 'ATPP': 'modifier', 'Km_ATPP': 'parameter'}}, #3 Convenience kinetics with comp. inhib. for UMPK_ATPP
                     {'equation': '(E_UMPK*(((kcat_F*(UMP/Km_UMP)*(ATP/(Km_ATP*(1+(ADP/Km_ADP)+(ATPP/Km_ATPP)+(AMP/Km_AMP)+(ATP/Km_ATP)))))-(kcat_F*(1/K_eq_UMPK)*((UDP*ADP)/(Km_UMP*(Km_ATP*(1+(ADP/Km_ADP)+(ATPP/Km_ATPP)+(AMP/Km_AMP)+(ATP/Km_ATP)))))))/(((1+(UMP/Km_UMP))*(1+(ATP/(Km_ATP*(1+(ADP/Km_ADP)+(ATPP/Km_ATPP)+(AMP/Km_AMP)+(ATP/Km_ATP))))))+((1+(UDP/Km_UDP))*(1+(ADP/(Km_ADP*(1+(ADP/Km_ADP)+(ATPP/Km_ATPP)+(AMP/Km_AMP)+(ATP/Km_ATP))))))-1)))', 'mapping': {'E_UMPK': 'modifier', 'kcat_F': 'parameter', 'UMP': 'substrate', 'Km_UMP': 'parameter', 'ATP': 'substrate', 'Km_ATP': 'parameter', 'UDP': 'product', 'Km_UDP': 'parameter', 'ADP': 'product', 'Km_ADP': 'parameter', 'K_eq_UMPK': 'parameter', 'ATPP': 'modifier', 'Km_ATPP': 'parameter', 'AMP': 'modifier', 'Km_AMP': 'parameter'}}], #4 Convenience kinetics with comp. inhib. for UMPK_ADP and UMPK_ATPP
        'UMPK_ADP': [{'equation': '0', 'mapping': {}}, #0
                     {'equation': '(E_UMPK*(((kcat_F*(UMP/Km_UMP)*(ADP/(Km_ADP*(1+(ATP/Km_ATP)+(ADP/Km_ADP)))))-(kcat_F*(1/K_eq_UMPK)*((UDP*AMP)/(Km_UMP*(Km_ADP*(1+(ATP/Km_ATP)+(ADP/Km_ADP)))))))/(((1+(UMP/Km_UMP))*(1+(ADP/(Km_ADP*(1+(ATP/Km_ATP)+(ADP/Km_ADP))))))+((1+(UDP/Km_UDP))*(1+(AMP/(Km_AMP*(1+(ATP/Km_ATP)+(ADP/Km_ADP))))))-1)))', 'mapping': {'E_UMPK': 'modifier', 'kcat_F': 'parameter', 'UMP': 'substrate', 'Km_UMP': 'parameter', 'ADP': 'substrate', 'Km_ADP': 'parameter', 'UDP': 'product', 'Km_UDP': 'parameter', 'AMP': 'product', 'Km_AMP': 'parameter', 'K_eq_UMPK': 'parameter', 'ATP': 'modifier', 'Km_ATP': 'parameter'}}, #1 Convenience kinetics with comp. inhib. for UMPK_ATP
                     {'equation': '(E_UMPK*(((kcat_F*(UMP/Km_UMP)*(ADP/(Km_ADP*(1+2*(ATP/Km_ATP)+(ATPP/Km_ATPP)+(ADP/Km_ADP)))))-(kcat_F*(1/K_eq_UMPK)*((UDP*AMP)/(Km_UMP*(Km_ADP*(1+2*(ATP/Km_ATP)+(ATPP/Km_ATPP)+(ADP/Km_ADP)))))))/(((1+(UMP/Km_UMP))*(1+(ADP/(Km_ADP*(1+2*(ATP/Km_ATP)+(ATPP/Km_ATPP)+(ADP/Km_ADP))))))+((1+(UDP/Km_UDP))*(1+(AMP/(Km_AMP*(1+2*(ATP/Km_ATP)+(ATPP/Km_ATPP)+(ADP/Km_ADP))))))-1)))', 'mapping': {'E_UMPK': 'modifier', 'kcat_F': 'parameter', 'UMP': 'substrate', 'Km_UMP': 'parameter', 'ADP': 'substrate', 'Km_ADP': 'parameter', 'UDP': 'product', 'Km_UDP': 'parameter', 'AMP': 'product', 'Km_AMP': 'parameter', 'K_eq_UMPK': 'parameter', 'ATP': 'modifier', 'Km_ATP': 'parameter', 'ATPP': 'modifier', 'Km_ATPP': 'parameter'}}], #2 Convenience kinetics with comp. inhib. for UMPK_ATP and UMPK_ATPP
        'UMPK_ATPP': [{'equation': '0', 'mapping': {}}, #0
                      {'equation': '(E_UMPK*(((kcat_F*(UMP/Km_UMP)*(ATPP/(Km_ATPP*(1+(ATP/Km_ATP)+(ADP/Km_ADP)))))-(kcat_F*(1/K_eq_UMPK)*((UDP*ATP)/(Km_UMP*(Km_ATPP*(1+(ATP/Km_ATP)+(ADP/Km_ADP)))))))/(((1+(UMP/Km_UMP))*(1+(ATPP/(Km_ATPP*(1+(ATP/Km_ATP)+(ADP/Km_ADP))))))+((1+(UDP/Km_UDP))*(1+(ATP/(Km_ATP*(1+(ATP/Km_ATP)+(ADP/Km_ADP))))))-1)))', 'mapping': {'E_UMPK': 'modifier', 'kcat_F': 'parameter', 'UMP': 'substrate', 'Km_UMP': 'parameter', 'ATPP': 'substrate', 'Km_ATPP': 'parameter', 'UDP': 'product', 'Km_UDP': 'parameter', 'ATP': 'product', 'Km_ATP': 'parameter', 'K_eq_UMPK': 'parameter', 'ADP': 'modifier', 'Km_ADP': 'parameter'}}, #1 Convenience kinetics with comp. inhib. for UMPK_ATP
                      {'equation': '(E_UMPK*(((kcat_F*(UMP/Km_UMP)*(ATPP/(Km_ATPP*(1+(ATP/Km_ATP)+2*(ADP/Km_ADP)+(AMP/Km_AMP)))))-(kcat_F*(1/K_eq_UMPK)*((UDP*ATP)/(Km_UMP*(Km_ATPP*(1+(ATP/Km_ATP)+2*(ADP/Km_ADP)+(AMP/Km_AMP)))))))/(((1+(UMP/Km_UMP))*(1+(ATPP/(Km_ATPP*(1+(ATP/Km_ATP)+2*(ADP/Km_ADP)+(AMP/Km_AMP))))))+((1+(UDP/Km_UDP))*(1+(ATP/(Km_ATP*(1+(ATP/Km_ATP)+2*(ADP/Km_ADP)+(AMP/Km_AMP))))))-1)))', 'mapping': {'E_UMPK': 'modifier', 'kcat_F': 'parameter', 'UMP': 'substrate', 'Km_UMP': 'parameter', 'ATPP': 'substrate', 'Km_ATPP': 'parameter', 'UDP': 'product', 'Km_UDP': 'parameter', 'ATP': 'product', 'Km_ATP': 'parameter', 'K_eq_UMPK': 'parameter', 'ADP': 'modifier', 'Km_ADP': 'parameter', 'AMP': 'modifier', 'Km_AMP': 'parameter'}}] #2 Convenience kinetics with comp. inhib. for UMPK_ATP and UMPK_ADP
    }
# regarding the regulation terms: all model species are represented as possible activators and possible inhibitors; the inhibition equation is based on Liebermeister and Klipp, 2006 <https://pmc.ncbi.nlm.nih.gov/articles/PMC1781438/> while the activation equation is a custom approach with a redefined ka parameter (custom ka = 1/(normal ka); this way the activation can be completely turned off by setting ka to 0)
##v23+: add new regulation terms to the library: Hill inhibition via ATP and ADP
    regulatory_terms = [
        {'name': 'no_reg',  # 0
         'equation': '1',
         'mapping': {}},
        {'name': 'ADP_Act',  # 1
         'equation': '(1+(ka_ADP*ADP))',
         'mapping': {'ka_ADP': 'parameter', 'ADP': 'modifier'}},
        {'name': 'ADP_Inhib',  # 2
         'equation': '((ki_ADP)/(ki_ADP+ADP))',
         'mapping': {'ki_ADP': 'parameter', 'ADP': 'modifier'}},
        {'name': 'AMP_Act',  # 3
         'equation': '(1+(ka_AMP*AMP))',
         'mapping': {'ka_AMP': 'parameter', 'AMP': 'modifier'}},
        {'name': 'AMP_Inhib',  # 4
         'equation': '((ki_AMP)/(ki_AMP+AMP))',
         'mapping': {'ki_AMP': 'parameter', 'AMP': 'modifier'}},
        {'name': 'ATP_Act',  # 5
         'equation': '(1+(ka_ATP*ATP))',
         'mapping': {'ka_ATP': 'parameter', 'ATP': 'modifier'}},
        {'name': 'ATP_Inhib',  # 6
         'equation': '((ki_ATP)/(ki_ATP+ATP))',
         'mapping': {'ki_ATP': 'parameter', 'ATP': 'modifier'}},
        {'name': 'GalNAc_Act',  # 7
         'equation': '(1+(ka_GalNAc*GalNAc))',
         'mapping': {'ka_GalNAc': 'parameter', 'GalNAc': 'modifier'}},
        {'name': 'GalNAc_Inhib',  # 8
         'equation': '((ki_GalNAc)/(ki_GalNAc+GalNAc))',
         'mapping': {'ki_GalNAc': 'parameter', 'GalNAc': 'modifier'}},
        {'name': 'P_Act',  # 9
         'equation': '(1+(ka_P*P))',
         'mapping': {'ka_P': 'parameter', 'P': 'modifier'}},
        {'name': 'P_Inhib',  # 10
         'equation': '((ki_P)/(ki_P+P))',
         'mapping': {'ki_P': 'parameter', 'P': 'modifier'}},
        {'name': 'PP_Act',  # 11
         'equation': '(1+(ka_PP*PP))',
         'mapping': {'ka_PP': 'parameter', 'PP': 'modifier'}},
        {'name': 'PP_Inhib',  # 12
         'equation': '((ki_PP)/(ki_PP+PP))',
         'mapping': {'ki_PP': 'parameter', 'PP': 'modifier'}},
        {'name': 'UDP_Act',  # 13
         'equation': '(1+(ka_UDP*UDP))',
         'mapping': {'ka_UDP': 'parameter', 'UDP': 'modifier'}},
        {'name': 'UDP_Inhib',  # 14
         'equation': '((ki_UDP)/(ki_UDP+UDP))',
         'mapping': {'ki_UDP': 'parameter', 'UDP': 'modifier'}},
        {'name': 'UDP_GalNAc_Act',  # 15
         'equation': '(1+(ka_UDP_GalNAc*UDP_GalNAc))',
         'mapping': {'ka_UDP_GalNAc': 'parameter', 'UDP_GalNAc': 'modifier'}},
        {'name': 'UDP_GalNAc_Inhib',  # 16
         'equation': '((ki_UDP_GalNAc)/(ki_UDP_GalNAc+UDP_GalNAc))',
         'mapping': {'ki_UDP_GalNAc': 'parameter', 'UDP_GalNAc': 'modifier'}},
        {'name': 'UMP_Act',  # 17
         'equation': '(1+(ka_UMP*UMP))',
         'mapping': {'ka_UMP': 'parameter', 'UMP': 'modifier'}},
        {'name': 'UMP_Inhib',  # 18
         'equation': '((ki_UMP)/(ki_UMP+UMP))',
         'mapping': {'ki_UMP': 'parameter', 'UMP': 'modifier'}},
        {'name': 'Uri_Act',  # 19
         'equation': '(1+(ka_Uri*Uri))',
         'mapping': {'ka_Uri': 'parameter', 'Uri': 'modifier'}},
        {'name': 'Uri_Inhib',  # 20
         'equation': '((ki_Uri)/(ki_Uri+Uri))',
         'mapping': {'ki_Uri': 'parameter', 'Uri': 'modifier'}},
        {'name': 'UTP_Act',  # 21
         'equation': '(1+(ka_UTP*UTP))',
         'mapping': {'ka_UTP': 'parameter', 'UTP': 'modifier'}},
        {'name': 'UTP_Inhib',  # 22
         'equation': '((ki_UTP)/(ki_UTP+UTP))',
         'mapping': {'ki_UTP': 'parameter', 'UTP': 'modifier'}},
        {'name': 'PolyP_Act',  # 23
         'equation': '(1+(ka_PolyP*PolyP))',
         'mapping': {'ka_PolyP': 'parameter', 'PolyP': 'modifier'}},
        {'name': 'PolyP_Inhib',  # 24
         'equation': '((ki_PolyP)/(ki_PolyP+PolyP))',
         'mapping': {'ki_PolyP': 'parameter', 'PolyP': 'modifier'}},
        {'name': 'ATPP_Act',  # 25
         'equation': '(1+(ka_ATPP*ATPP))',
         'mapping': {'ka_ATPP': 'parameter', 'ATPP': 'modifier'}},
        {'name': 'ATPP_Inhib',  # 26
         'equation': '((ki_ATPP)/(ki_ATPP+ATPP))',
         'mapping': {'ki_ATPP': 'parameter', 'ATPP': 'modifier'}},
        {'name': 'ATP_Inhib_Hill',  # 27
         'equation': '((ki_ATP^n_Hill)/((ki_ATP^n_Hill)+(ATP^n_Hill)))',
         'mapping': {'ki_ATP': 'parameter', 'n_Hill': 'parameter', 'ATP': 'modifier'}},
        {'name': 'ADP_Inhib_Hill',  # 28
         'equation': '((ki_ADP^n_Hill)/((ki_ADP^n_Hill)+(ADP^n_Hill)))',
         'mapping': {'ki_ADP': 'parameter', 'n_Hill': 'parameter', 'ADP': 'modifier'}}
    ]

    return base_kinetics, regulatory_terms


# Simple toymodel (for the rediscovery test)

def toymodel_data():
    """ Define (1) string representations (schemes) of all reactions of the toymodel according to the syntax used by COPASI, (2) a dictionary of parameter values, and (3) a dictionary of initial concentrations for all species. Model units are hours, mmol, and liters.
    
    :return reaction_scheme_dict:
    :type reaction_scheme_dict:
    :return param_dict:
    :type param_dict:
    :return init_conc_dict:
    :type init_conc_dict:
    """

    # define reaction schemes; all of them are set even if a reaction is disabled in the structural variant matrix (in this case the underlying rate law is set to 0); the schemes are used as default and are modified according to the variable terms set in the structural variant matrix; modifiers that are not used in the selected rate law are inactive
    reaction_scheme_dict = {'r1': 'S = P',
                            'r2': 'S = X',
                            'r3': 'P = X',
                            'r4': 'S + P = X'}

    # collect known literature values of kinetic parameters; not all parameters have to appear here - some may not be used depending on the kinetics that are selected according to the structural variant matrix; kcat [1/h]; Km, ki, ka [mM]; K_eq [-]
    param_dict = {'r1': {},
                  'r2': {},
                  'r3': {},
                  'r4': {}}

    # define all non-zero initial species concentrations (the initial concentrations for the remaining species will be set to a very low value (1e-12) - setting them to zero can make problems with certain rate laws); initial concentrations are just used as a placeholder when creating the COPASI model object - during the estimation loop basico changes the initial species concentrations to those found in the provided experimental data files
    init_conc_dict = {'S': 10}

    return reaction_scheme_dict, param_dict, init_conc_dict


def toymodel_term_libs():
    """ Define the libraries of model specific base and regulation terms. Base terms are stored in a dictionary where the keys are reaction names and the values are lists of dictionaries; each of them contains a kinetic rate law equation and the associated mapping of all variables and parameters. Regulation terms are stored in a list of dictionaries which contain names, equations, and associated mappings.
    
    :return base_kinetics: dictionary of base kinetics (equations and mappings) per reaction
    :type base_kinetics: dictionary
    :return regulatory_terms: list of regulation terms (each entry is a dictionary with the name, equation, and mapping of the regulation term)
    :type regulatory_terms: list
    """

# create repositories of reaction specific base kinetics and general regulatory terms that can be indexed according to the integers contained in the structural variant vector (0: no base kinetics <-> reaction turned off, 1: first base kinetics in sublist, and so on); all kinetics come with dictionaries mapping variable names to their usage, i.e., 'parameter', 'substrate', 'product', or 'modifier'
    base_kinetics = {
        'r1': [{'equation': '0', 'mapping': {}}, #0
                         {'equation': '(k_MAforward*S-k_MAreverse*P)', 'mapping': {'k_MAforward': 'parameter', 'k_MAreverse': 'parameter', 'S': 'substrate', 'P': 'product'}}], #1
        'r2': [{'equation': '0', 'mapping': {}}, #0
                         {'equation': '(k_MAforward*S-k_MAreverse*X)', 'mapping': {'k_MAforward': 'parameter', 'k_MAreverse': 'parameter', 'S': 'substrate', 'X': 'product'}}], #1
        'r3': [{'equation': '0', 'mapping': {}}, #0
                         {'equation': '(k_MAforward*P-k_MAreverse*X)', 'mapping': {'k_MAforward': 'parameter', 'k_MAreverse': 'parameter', 'P': 'substrate', 'X': 'product'}}], #1
        'r4': [{'equation': '0', 'mapping': {}}, #0
                         {'equation': '(k_MAforward*S*P-k_MAreverse*X)', 'mapping': {'k_MAforward': 'parameter', 'k_MAreverse': 'parameter', 'S': 'substrate', 'P': 'substrate', 'X': 'product'}}], #1
    }
# regarding the regulation terms: all model species are represented as possible activators and possible inhibitors; the inhibition equation is based on Liebermeister and Klipp, 2006 <https://pmc.ncbi.nlm.nih.gov/articles/PMC1781438/> while the activation equation is a custom approach with a redefined ka parameter (custom ka = 1/(normal ka); this way the activation can be completely turned off by setting ka to 0)
    regulatory_terms = [
        {'name': 'no_reg',  # 0
         'equation': '1',
         'mapping': {}},
        {'name': 'S_Act',  # 1
         'equation': '(1+(ka_S*S))',
         'mapping': {'ka_S': 'parameter', 'S': 'modifier'}},
        {'name': 'S_Inhib',  # 2
         'equation': '((ki_S)/(ki_S+S))',
         'mapping': {'ki_S': 'parameter', 'S': 'modifier'}},
        {'name': 'P_Act',  # 3
         'equation': '(1+(ka_P*P))',
         'mapping': {'ka_P': 'parameter', 'P': 'modifier'}},
        {'name': 'P_Inhib',  # 4
         'equation': '((ki_P)/(ki_P+P))',
         'mapping': {'ki_P': 'parameter', 'P': 'modifier'}},
        {'name': 'X_Act',  # 5
         'equation': '(1+(ka_X*X))',
         'mapping': {'ka_X': 'parameter', 'X': 'modifier'}},
        {'name': 'X_Inhib',  # 6
         'equation': '((ki_X)/(ki_X+X))',
         'mapping': {'ki_X': 'parameter', 'X': 'modifier'}},
    ]

    return base_kinetics, regulatory_terms

# ------------------------------------------------------------------------------------------------------ #

# prepare experimental data for any task that involves parameter estimations (UDP-GalNAc model specific!)
def load_and_process_exp_data(exp_data_file_names):
    """ Load experimental data then change the following things for each data frame: (a) add up the measurements of adenosine tetra- and penta- phosphates (ATPP and ATPPP); (b) scale PolyP to get the amount of available phosphate bonds; (c) replace all zeros with 1e-12 to avoid incompatibilities with certain rate laws like the reversible Hill rate equation.

    :param exp_data_file_names: list of experimental data file name strings (with file extensions)
    :type exp_data_file_names: list

    :return exp_data_dataframes: list of pandas data frames that contain the modified (pre-processed) experimental data
    :type exp_data_dataframes: list
    """

    # create a pandas data frame from each data file
    exp_data_dataframes = list()
    for file_name in exp_data_file_names:
        # load the raw data
        exp_data_df = pd.read_csv(file_name, sep='\t')
        # change the data (pre-processing)
        # (a) add up the measurements of adenosine tetra- and penta- phosphates (ATPP, and ATPPP) in the ATPP column
        exp_data_df['temp_ATPP_total'] = exp_data_df[['[ATPP]', '[ATPPP]']].sum(axis=1)
        # remove the old individual columns that were summed up
        exp_data_df = exp_data_df.drop(columns=['[ATPP]', '[ATPPP]'])
        # change the temporary name to the proper name that Copasi can recognize as dependent species
        exp_data_df.rename(columns={'temp_ATPP_total': '[ATPP]'}, inplace=True)
        # reorder the columns so that the new ATPP column is inserted after the ATP column
        exp_data_df = exp_data_df[['Time', '[Uri]', '[UMP]', '[UDP_GalNAc]', '[UDP]', '[UTP]', '[AMP]', '[ADP]', '[ATP]', '[ATPP]', '[GalNAc]_0', '[ATP]_0', '[ADP]_0', '[Uri]_0', '[E_PPA]_0', '[E_PPK3]_0', '[E_UDK]_0', '[E_NAHK]_0', '[E_GLMU]_0', '[E_UMPK]_0', '[PolyP]_0']]
        # (b) scale PolyP by *8 to get the amount of available phosphate bonds
        exp_data_df.at[0, '[PolyP]_0'] = exp_data_df.loc[0, '[PolyP]_0']*8
        # (c) replace all zeros with 1e-12 to avoid incompatibilities with certain rate laws like the reversible Hill rate equation
        for col in exp_data_df.columns:
            if col != 'Time':
                exp_data_df.replace({col: 0}, 1e-12, inplace=True)
        # store the modified data frame
        exp_data_dataframes.append(exp_data_df)

    return exp_data_dataframes

# ------------------------------------------------------------------------------------------------------ #

# MODEL SELECTION - SUPER STRUCTURE OPTIMIZATION (SSO) FUNCTIONS

# core function to create a COPASI model object from a model variant data frame
def create_model_obj_from_struct_var(struct_var, model_data=UDP_GalNAc_model_data, term_libs=UDP_GalNAc_model_term_libs):
    """ Create a COPASI model object with kinetics derived from the provided structural variant matrix (pandas data frame). Integers are used as indices to look up base kinetics and regulation terms in the respective libraries. The model object and the new reaction schemes are returned.

    :param struct_var: a structural variant matrix
    :type struct_var: pandas.core.frame.DataFrame
    :param model_data: model-specific function that contains information on the reaction schemes, known literature values for reaction parameters, and placeholder values for all non-zero initial concentrations (default: UDP_GalNAc model)
    :type model_data: function
    :param term_libs: model-specific function that contains all potential base rate law equations as well as all potential regulation terms (default: UDP_GalNAc model)
    :type term_libs: function

    :return model_object: a Copasi model object
    :type model_object: COPASI.CDataModel
    :return new_reaction_schemes_dict: a dictionary of the new reaction schemes
    :type new_reaction_schemes_dict: dictionary
    """

    # -------------------------------------------------------------------
    # LOAD MODEL SPECIFIC LIBRARIES OF BASE KINETICS AND REGULATORY TERMS
    # -------------------------------------------------------------------
    base_kinetics, regulatory_terms = term_libs()

    # -----------------------------------------------------
    # DECODE STRUCTURAL VARIANT VECTOR AND CREATE RATE LAWS
    # -----------------------------------------------------
    # convert all entries in the structural variant data frame to integers (they should already be integers but pygmo likes to convert them to floats even if integers are supplied, so we have to make sure that we work with integers here because they can are used as indices)
    for i in struct_var.columns:
        try:
            struct_var[[i]] = struct_var[[i]].astype(int)
        except (Exception,):
            pass
    # iterate over rows of the structural variant matrix: for each row vector take the first entry as base kinetics integer and then check if any other entries exist - if yes then columns with regulation terms exist and the regulations terms that they specify need to be concatenated to the rate law equation as multipliers (the associated mappings also need to be updated)
    combined_eq_list = list()
    func_to_usage_mapping_dict = dict()
    # decode structural variant integer list by combining base kinetics with selected regulatory terms 
    for row_name, row_vector in struct_var.iterrows():
        base_term = base_kinetics[row_name][row_vector.iloc[0]]['equation']
        # only build rate laws if there's a non-zero base term
        if base_term != '0':
            # if there are more entries in the row vector beyond the one for the base term: look them up and add them to the equation as multiplicative factors
            combined_eq = None
            func_to_usage_mapping = dict()
            if len(row_vector) > 1:
                reg_terms = list()
                for int_val in row_vector[1::]:
                    reg_terms.append(regulatory_terms[int_val]['equation'])
                    func_to_usage_mapping.update(regulatory_terms[int_val]['mapping'])
                combined_eq = (base_term + '*' + '*'.join(reg_terms))
            else:
                combined_eq = base_term
            combined_eq_list.append(combined_eq)
            # merge [rate law variable name to usage] mapping dictionaries of base kinetics and selected regulatory terms (or empty dict) for all reactions; append base kinetics dictionary after regulatory terms dictionaries have been added so that species which already appear in the regulatory terms as 'modifiers' are properly set as 'substrate' or 'product' if they appear as such in the base kinetics
            func_to_usage_mapping.update(base_kinetics[row_name][row_vector.iloc[0]]['mapping'])
            func_to_usage_mapping_dict[row_name] = func_to_usage_mapping
        # if the base term is zero then the entire rate law must be zero (since any regulation terms are just multiplied with the base term so if it's zero they don't matter)
        elif base_term == '0':
            combined_eq = '0'
            combined_eq_list.append(combined_eq)
            func_to_usage_mapping = {}
            func_to_usage_mapping_dict[row_name] = func_to_usage_mapping

    # --------------------------
    # CREATE COPASI MODEL OBJECT
    # --------------------------
    # remove any variant functions from the global function list if they exist (there might be remnants from previous iterations that need to be removed so that we can start from a 'clean' slate)
    for name in sorted(list(get_functions().index)):
        if '_Variant' in name:
            remove_function(name)
    # add rate law functions created from newly combined equation strings and updated [function to usage] mapping dictionaries to global function list; get global function list with get_functions()
    for (name, eq, mapping) in zip(list(struct_var.index),
                                   combined_eq_list,
                                   list(func_to_usage_mapping_dict.values())):
        reaction_type = str()
        if name == 'PPA':
            reaction_type = 'irreversible'
        else:
            reaction_type = 'reversible'
        add_function(name + '_Variant', eq, type=reaction_type, mapping=mapping)
    # load model specific data
    default_reaction_schemes_dict, param_dict, init_conc_dict = model_data()
    # select only those default reaction schemes for reactions that are present in the structural variant
    default_reaction_schemes_dict = {key: value for key, value in default_reaction_schemes_dict.items() if
                                     key in list(struct_var.index)}
    # update default reaction scheme strings with selected regulators: first get scheme without any mods via partition(), then build list of selected regulators from updated mapping dictionary, then remove those regulator strings that are already part of the core scheme (i.e., only keep the 'unique' modifiers) and then add the rest after the ';' build the new scheme
    new_reaction_schemes_dict = dict()
    for (reaction_name, reaction_scheme, func_to_usage_mapping) in zip(list(struct_var.index),
                                                                       list(default_reaction_schemes_dict.values()),
                                                                       list(func_to_usage_mapping_dict.values())):
        reaction_scheme_core, _, _ = reaction_scheme.partition(';')
        selected_regulators = [k for k, v in func_to_usage_mapping.items() if v == 'modifier']
        if len(selected_regulators) != 0:
            new_scheme = reaction_scheme_core + '; ' + ' '.join(selected_regulators)
        else:
            # if no modifiers (neither enzymes nor regulators) were added to the reaction then the new scheme is the same as the core
            new_scheme = reaction_scheme_core
        new_reaction_schemes_dict[reaction_name] = new_scheme
    # create [function to reaction] mapping dictionaries for all reactions; the keys are the function variables, the values are either the parameter values (for parameters) or the names used in the reaction scheme (for species); the parameter values that are set here are mainly placeholders and will be overwritten by eval_struct_var in case the model is used for a parameter estimation
    func_to_reac_mapping_dict = dict()
    for reaction_name in list(struct_var.index):
        func_to_reac_mapping = dict()
        for var_name, var_type in func_to_usage_mapping_dict[reaction_name].items():
            if var_type == 'parameter':
                # set all ki and ka parameters to values that correspond to low regulation strengths (regulation term close or equal to 1):
                # - inhibition terms: (ki/(ki+[I])) -> very large ki = term close to 1 => very low inhibition strength (close to 0%)
                # - activation terms: (1+ka*[A]) -> ka set to 0 = term equal to 1 => activation strength equal to 0% (no activation at all)
                if 'ki_'in var_name:
                    func_to_reac_mapping[var_name] = 1e04 # (1e04/(1e04+[I])) => very close to 1
                elif 'ka_' in var_name:
                    func_to_reac_mapping[var_name] = 0 # (1+(0*[A])) => exactly at 1
                # for all other kinds of parameters: use either the available literature value or 0.1 as 
                # fall back value
                else:
                    if var_name in param_dict[reaction_name]:
                        func_to_reac_mapping[var_name] = param_dict[reaction_name][var_name]
                    else:
                        func_to_reac_mapping[var_name] = 0.1  # COPASI default parameter value
            else:
                # all non-parameter variables are species and their names are the same in the functions and the schemes
                func_to_reac_mapping[var_name] = var_name
        func_to_reac_mapping_dict[reaction_name] = func_to_reac_mapping
    # initialize an empty model object
    model_object = new_model()
    set_model_unit(time_unit='h', substance_unit='mmol', volume_unit='l', model=model_object)
    # add reactions to model object
    for name, scheme in new_reaction_schemes_dict.items():
        add_reaction(name=name,
                     scheme=scheme,
                     function=name + '_Variant',
                     model=model_object)
        # set reaction mapping: species variables (substrates, products, modifiers) are the same in the function and the reaction scheme; for all parameter variables their float values are assigned; due to a bug in basico the mapping needs to be done separately from the add_reaction call AND one by one (instead of passing a complete dictionary with all mapping information) => only then will the mapping be correct
        for key, val in func_to_reac_mapping_dict[name].items():
            set_reaction_mapping(reaction=name,
                                 mapping={key: val},
                                 model=model_object)

    #--------------------
    # MODEL-SPECIFIC CODE
    #--------------------
    # only run these code blocks if the right model-specific functions were loaded
    if model_data.__name__ == 'UDP_GalNAc_model_data' and term_libs.__name__ == 'UDP_GalNAc_model_term_libs':
        # create global model objects (global quantities) that can be used to share selected parameters between specific reaction groups (in these cases the same enzyme catalyzes the reversible conversion of multiple competing substrates and products so their Km values are used to modify the Km values of the substrates and products that are affected by the competition in the corresponding reactions)
        # (1) PPK3 catalyzes the conversion of ADP to ATP, UDP to UTP, and ATP to ATPP (adenosine tetraphosphate); allowed combinations of base kinetics: (1,1,0)[= PPK3_A and PPK3_U are both active and parameters are shared for competitive inhibition, PPK3_tetra is inactive] and (2,2,1)[= PPK3_A, PPK3_U, and PPK3_tetra are all active and parameters are shared for competitive inhibition]
        if (struct_var.loc['PPK3_A']['base'] == 1 and struct_var.loc['PPK3_U']['base'] == 1 and struct_var.loc['PPK3_tetra']['base'] == 0):
            # first the global quantity needs to be initialized with a placeholder value, then it can be changed to the proper form (here, 'parameter' refers to a global quantity; sometimes also called 'model value')
            add_parameter('(PPK3_U).Km_UDP_Link', 
                          value=get_reaction_parameters().loc['(PPK3_U).Km_UDP', 'value'])
            set_parameters('(PPK3_U).Km_UDP_Link',
                           type='assignment',
                           expression="(PPK3_U).Km_UDP",
                           model=model_object)
            add_parameter('(PPK3_U).Km_UTP_Link', 
                          value=get_reaction_parameters().loc['(PPK3_U).Km_UTP', 'value'])
            set_parameters('(PPK3_U).Km_UTP_Link',
                           type='assignment',
                           expression="(PPK3_U).Km_UTP",
                           model=model_object)
            add_parameter('(PPK3_A).Km_ADP_Link',
                          value=get_reaction_parameters().loc['(PPK3_A).Km_ADP', 'value'])
            set_parameters('(PPK3_A).Km_ADP_Link',
                           type='assignment',
                           expression="(PPK3_A).Km_ADP",
                           model=model_object)
            add_parameter('(PPK3_A).Km_ATP_Link',
                          value=get_reaction_parameters().loc['(PPK3_A).Km_ATP', 'value'])
            set_parameters('(PPK3_A).Km_ATP_Link',
                           type='assignment',
                           expression="(PPK3_A).Km_ATP",
                           model=model_object)
            # update mapping of selected reactions with shared parameters via global quantities
            set_reaction_mapping(reaction='PPK3_U',
                         mapping={'Km_ADP': '(PPK3_A).Km_ADP_Link', 'Km_ATP': '(PPK3_A).Km_ATP_Link'},
                         model=model_object)
            set_reaction_mapping(reaction='PPK3_A',
                         mapping={'Km_UDP': '(PPK3_U).Km_UDP_Link', 'Km_UTP': '(PPK3_U).Km_UTP_Link'},
                         model=model_object)
        if (struct_var.loc['PPK3_A']['base'] == 2 and struct_var.loc['PPK3_U']['base'] == 2 and struct_var.loc['PPK3_tetra']['base'] == 1):
            # first the global quantity needs to be initialized with a placeholder value, then it can be changed to the proper form (here, 'parameter' refers to a global quantity; sometimes also called 'model value')
            add_parameter('(PPK3_U).Km_UDP_Link', 
                          value=get_reaction_parameters().loc['(PPK3_U).Km_UDP', 'value'])
            set_parameters('(PPK3_U).Km_UDP_Link',
                           type='assignment',
                           expression="(PPK3_U).Km_UDP",
                           model=model_object)
            add_parameter('(PPK3_U).Km_UTP_Link', 
                          value=get_reaction_parameters().loc['(PPK3_U).Km_UTP', 'value'])
            set_parameters('(PPK3_U).Km_UTP_Link',
                           type='assignment',
                           expression="(PPK3_U).Km_UTP",
                           model=model_object)
            add_parameter('(PPK3_A).Km_ADP_Link',
                          value=get_reaction_parameters().loc['(PPK3_A).Km_ADP', 'value'])
            set_parameters('(PPK3_A).Km_ADP_Link',
                           type='assignment',
                           expression="(PPK3_A).Km_ADP",
                           model=model_object)
            add_parameter('(PPK3_A).Km_ATP_Link',
                          value=get_reaction_parameters().loc['(PPK3_A).Km_ATP', 'value'])
            set_parameters('(PPK3_A).Km_ATP_Link',
                           type='assignment',
                           expression="(PPK3_A).Km_ATP",
                           model=model_object)
            add_parameter('(PPK3_tetra).Km_ATPP_Link',
                          value=get_reaction_parameters().loc['(PPK3_tetra).Km_ATPP', 'value'])
            set_parameters('(PPK3_tetra).Km_ATPP_Link',
                           type='assignment',
                           expression="(PPK3_tetra).Km_ATPP",
                           model=model_object)
            # update mapping of selected reactions with shared parameters via global quantities
            set_reaction_mapping(reaction='PPK3_U',
                         mapping={'Km_ADP': '(PPK3_A).Km_ADP_Link', 'Km_ATP': '(PPK3_A).Km_ATP_Link', 'Km_ATPP': '(PPK3_tetra).Km_ATPP_Link'},
                         model=model_object)
            set_reaction_mapping(reaction='PPK3_A',
                         mapping={'Km_UDP': '(PPK3_U).Km_UDP_Link', 'Km_UTP': '(PPK3_U).Km_UTP_Link', 'Km_ATPP': '(PPK3_tetra).Km_ATPP_Link'},
                         model=model_object)
            set_reaction_mapping(reaction='PPK3_tetra',
                         mapping={'Km_ADP': '(PPK3_A).Km_ADP_Link', 'Km_ATP': '(PPK3_A).Km_ATP_Link', 'Km_UDP': '(PPK3_U).Km_UDP_Link', 'Km_UTP': '(PPK3_U).Km_UTP_Link'},
                         model=model_object)
        # (2) NAHK, UDK, and UMPK can utilize ATP, ADP, and ATPP as competing phosphate donors (and since the reactions are reversible ADP, AMP, and ATP compete as phosphate acceptors); allowed combinations of base kinetics: (1,0,0)[= _ATP version is active, _ADP and _ATPP versions are inactive - no parameter sharing necessary], (2,1,0)[= _ATP and _ADP versions are active and parameters are shared for competitive inhibition], (3,0,1)[= _ATP and _ATPP versions are active and parameters are shared for competitive inhibition], (4,2,2)[= _ATP, _ADP and _ATPP versions are active and parameters are shared for competitive inhibition]
        if (struct_var.loc['NAHK_ATP']['base'] == 2 and struct_var.loc['NAHK_ADP']['base'] == 1 and struct_var.loc['NAHK_ATPP']['base'] == 0):
            add_parameter('(NAHK_ATP).Km_ATP_Link', 
                          value=get_reaction_parameters().loc['(NAHK_ATP).Km_ATP', 'value'])
            set_parameters('(NAHK_ATP).Km_ATP_Link',
                           type='assignment',
                           expression="(NAHK_ATP).Km_ATP",
                           model=model_object)
            add_parameter('(NAHK_ATP).Km_ADP_Link', 
                          value=get_reaction_parameters().loc['(NAHK_ATP).Km_ADP', 'value'])
            set_parameters('(NAHK_ATP).Km_ADP_Link',
                           type='assignment',
                           expression="(NAHK_ATP).Km_ADP",
                           model=model_object)
            add_parameter('(NAHK_ADP).Km_AMP_Link', 
                          value=get_reaction_parameters().loc['(NAHK_ADP).Km_AMP', 'value'])
            set_parameters('(NAHK_ADP).Km_AMP_Link',
                           type='assignment',
                           expression="(NAHK_ADP).Km_AMP",
                           model=model_object)
            # update mapping of selected reactions with shared parameters via global quantities
            set_reaction_mapping(reaction='NAHK_ATP',
                         mapping={'Km_AMP': '(NAHK_ADP).Km_AMP_Link'},
                         model=model_object)
            set_reaction_mapping(reaction='NAHK_ADP',
                         mapping={'Km_ATP': '(NAHK_ATP).Km_ATP_Link', 'Km_ADP': '(NAHK_ATP).Km_ADP_Link'},
                         model=model_object)
        if (struct_var.loc['NAHK_ATP']['base'] == 3 and struct_var.loc['NAHK_ADP']['base'] == 0 and struct_var.loc['NAHK_ATPP']['base'] == 1):
            add_parameter('(NAHK_ATP).Km_ATP_Link', 
                          value=get_reaction_parameters().loc['(NAHK_ATP).Km_ATP', 'value'])
            set_parameters('(NAHK_ATP).Km_ATP_Link',
                           type='assignment',
                           expression="(NAHK_ATP).Km_ATP",
                           model=model_object)
            add_parameter('(NAHK_ATP).Km_ADP_Link', 
                          value=get_reaction_parameters().loc['(NAHK_ATP).Km_ADP', 'value'])
            set_parameters('(NAHK_ATP).Km_ADP_Link',
                           type='assignment',
                           expression="(NAHK_ATP).Km_ADP",
                           model=model_object)
            add_parameter('(NAHK_ATPP).Km_ATPP_Link', 
                          value=get_reaction_parameters().loc['(NAHK_ATPP).Km_ATPP', 'value'])
            set_parameters('(NAHK_ATPP).Km_ATPP_Link',
                           type='assignment',
                           expression="(NAHK_ATPP).Km_ATPP",
                           model=model_object)
            # update mapping of selected reactions with shared parameters via global quantities
            set_reaction_mapping(reaction='NAHK_ATP',
                         mapping={'Km_ATPP': '(NAHK_ATPP).Km_ATPP_Link'},
                         model=model_object)
            set_reaction_mapping(reaction='NAHK_ATPP',
                         mapping={'Km_ATP': '(NAHK_ATP).Km_ATP_Link' ,'Km_ADP': '(NAHK_ATP).Km_ADP_Link'},
                         model=model_object)
        if (struct_var.loc['NAHK_ATP']['base'] == 4 and struct_var.loc['NAHK_ADP']['base'] == 2 and struct_var.loc['NAHK_ATPP']['base'] == 2):
            add_parameter('(NAHK_ATP).Km_ATP_Link', 
                          value=get_reaction_parameters().loc['(NAHK_ATP).Km_ATP', 'value'])
            set_parameters('(NAHK_ATP).Km_ATP_Link',
                           type='assignment',
                           expression="(NAHK_ATP).Km_ATP",
                           model=model_object)
            add_parameter('(NAHK_ATP).Km_ADP_Link', 
                          value=get_reaction_parameters().loc['(NAHK_ATP).Km_ADP', 'value'])
            set_parameters('(NAHK_ATP).Km_ADP_Link',
                           type='assignment',
                           expression="(NAHK_ATP).Km_ADP",
                           model=model_object)
            add_parameter('(NAHK_ADP).Km_AMP_Link', 
                          value=get_reaction_parameters().loc['(NAHK_ADP).Km_AMP', 'value'])
            set_parameters('(NAHK_ADP).Km_AMP_Link',
                           type='assignment',
                           expression="(NAHK_ADP).Km_AMP",
                           model=model_object)
            add_parameter('(NAHK_ATPP).Km_ATPP_Link', 
                          value=get_reaction_parameters().loc['(NAHK_ATPP).Km_ATPP', 'value'])
            set_parameters('(NAHK_ATPP).Km_ATPP_Link',
                           type='assignment',
                           expression="(NAHK_ATPP).Km_ATPP",
                           model=model_object)
            # update mapping of selected reactions with shared parameters via global quantities
            set_reaction_mapping(reaction='NAHK_ATP',
                         mapping={'Km_AMP': '(NAHK_ADP).Km_AMP_Link', 'Km_ATPP': '(NAHK_ATPP).Km_ATPP_Link'},
                         model=model_object)
            set_reaction_mapping(reaction='NAHK_ADP',
                         mapping={'Km_ATP': '(NAHK_ATP).Km_ATP_Link', 'Km_ADP': '(NAHK_ATP).Km_ADP_Link', 'Km_ATPP': '(NAHK_ATPP).Km_ATPP_Link'},
                         model=model_object)
            set_reaction_mapping(reaction='NAHK_ATPP',
                         mapping={'Km_ATP': '(NAHK_ATP).Km_ATP_Link', 'Km_ADP': '(NAHK_ATP).Km_ADP_Link', 'Km_AMP': '(NAHK_ADP).Km_AMP_Link'},
                         model=model_object)
        if (struct_var.loc['UDK_ATP']['base'] == 2 and struct_var.loc['UDK_ADP']['base'] == 1 and struct_var.loc['UDK_ATPP']['base'] == 0):
            add_parameter('(UDK_ATP).Km_ATP_Link', 
                          value=get_reaction_parameters().loc['(UDK_ATP).Km_ATP', 'value'])
            set_parameters('(UDK_ATP).Km_ATP_Link',
                           type='assignment',
                           expression="(UDK_ATP).Km_ATP",
                           model=model_object)
            add_parameter('(UDK_ATP).Km_ADP_Link', 
                          value=get_reaction_parameters().loc['(UDK_ATP).Km_ADP', 'value'])
            set_parameters('(UDK_ATP).Km_ADP_Link',
                           type='assignment',
                           expression="(UDK_ATP).Km_ADP",
                           model=model_object)
            add_parameter('(UDK_ADP).Km_AMP_Link', 
                          value=get_reaction_parameters().loc['(UDK_ADP).Km_AMP', 'value'])
            set_parameters('(UDK_ADP).Km_AMP_Link',
                           type='assignment',
                           expression="(UDK_ADP).Km_AMP",
                           model=model_object)
            # update mapping of selected reactions with shared parameters via global quantities
            set_reaction_mapping(reaction='UDK_ATP',
                         mapping={'Km_AMP': '(UDK_ADP).Km_AMP_Link'},
                         model=model_object)
            set_reaction_mapping(reaction='UDK_ADP',
                         mapping={'Km_ATP': '(UDK_ATP).Km_ATP_Link', 'Km_ADP': '(UDK_ATP).Km_ADP_Link'},
                         model=model_object)
        if (struct_var.loc['UDK_ATP']['base'] == 3 and struct_var.loc['UDK_ADP']['base'] == 0 and struct_var.loc['UDK_ATPP']['base'] == 1):
            add_parameter('(UDK_ATP).Km_ATP_Link', 
                          value=get_reaction_parameters().loc['(UDK_ATP).Km_ATP', 'value'])
            set_parameters('(UDK_ATP).Km_ATP_Link',
                           type='assignment',
                           expression="(UDK_ATP).Km_ATP",
                           model=model_object)
            add_parameter('(UDK_ATP).Km_ADP_Link', 
                          value=get_reaction_parameters().loc['(UDK_ATP).Km_ADP', 'value'])
            set_parameters('(UDK_ATP).Km_ADP_Link',
                           type='assignment',
                           expression="(UDK_ATP).Km_ADP",
                           model=model_object)
            add_parameter('(UDK_ATPP).Km_ATPP_Link', 
                          value=get_reaction_parameters().loc['(UDK_ATPP).Km_ATPP', 'value'])
            set_parameters('(UDK_ATPP).Km_ATPP_Link',
                           type='assignment',
                           expression="(UDK_ATPP).Km_ATPP",
                           model=model_object)
            # update mapping of selected reactions with shared parameters via global quantities
            set_reaction_mapping(reaction='UDK_ATP',
                         mapping={'Km_ATPP': '(UDK_ATPP).Km_ATPP_Link'},
                         model=model_object)
            set_reaction_mapping(reaction='UDK_ATPP',
                         mapping={'Km_ATP': '(UDK_ATP).Km_ATP_Link' ,'Km_ADP': '(UDK_ATP).Km_ADP_Link'},
                         model=model_object)
        if (struct_var.loc['UDK_ATP']['base'] == 4 and struct_var.loc['UDK_ADP']['base'] == 2 and struct_var.loc['UDK_ATPP']['base'] == 2):
            add_parameter('(UDK_ATP).Km_ATP_Link', 
                          value=get_reaction_parameters().loc['(UDK_ATP).Km_ATP', 'value'])
            set_parameters('(UDK_ATP).Km_ATP_Link',
                           type='assignment',
                           expression="(UDK_ATP).Km_ATP",
                           model=model_object)
            add_parameter('(UDK_ATP).Km_ADP_Link', 
                          value=get_reaction_parameters().loc['(UDK_ATP).Km_ADP', 'value'])
            set_parameters('(UDK_ATP).Km_ADP_Link',
                           type='assignment',
                           expression="(UDK_ATP).Km_ADP",
                           model=model_object)
            add_parameter('(UDK_ADP).Km_AMP_Link', 
                          value=get_reaction_parameters().loc['(UDK_ADP).Km_AMP', 'value'])
            set_parameters('(UDK_ADP).Km_AMP_Link',
                           type='assignment',
                           expression="(UDK_ADP).Km_AMP",
                           model=model_object)
            add_parameter('(UDK_ATPP).Km_ATPP_Link', 
                          value=get_reaction_parameters().loc['(UDK_ATPP).Km_ATPP', 'value'])
            set_parameters('(UDK_ATPP).Km_ATPP_Link',
                           type='assignment',
                           expression="(UDK_ATPP).Km_ATPP",
                           model=model_object)
            # update mapping of selected reactions with shared parameters via global quantities
            set_reaction_mapping(reaction='UDK_ATP',
                         mapping={'Km_AMP': '(UDK_ADP).Km_AMP_Link', 'Km_ATPP': '(UDK_ATPP).Km_ATPP_Link'},
                         model=model_object)
            set_reaction_mapping(reaction='UDK_ADP',
                         mapping={'Km_ATP': '(UDK_ATP).Km_ATP_Link', 'Km_ADP': '(UDK_ATP).Km_ADP_Link', 'Km_ATPP': '(UDK_ATPP).Km_ATPP_Link'},
                         model=model_object)
            set_reaction_mapping(reaction='UDK_ATPP',
                         mapping={'Km_ATP': '(UDK_ATP).Km_ATP_Link', 'Km_ADP': '(UDK_ATP).Km_ADP_Link', 'Km_AMP': '(UDK_ADP).Km_AMP_Link'},
                         model=model_object)
        if (struct_var.loc['UMPK_ATP']['base'] == 2 and struct_var.loc['UMPK_ADP']['base'] == 1 and struct_var.loc['UMPK_ATPP']['base'] == 0):
            add_parameter('(UMPK_ATP).Km_ATP_Link', 
                          value=get_reaction_parameters().loc['(UMPK_ATP).Km_ATP', 'value'])
            set_parameters('(UMPK_ATP).Km_ATP_Link',
                           type='assignment',
                           expression="(UMPK_ATP).Km_ATP",
                           model=model_object)
            add_parameter('(UMPK_ATP).Km_ADP_Link', 
                          value=get_reaction_parameters().loc['(UMPK_ATP).Km_ADP', 'value'])
            set_parameters('(UMPK_ATP).Km_ADP_Link',
                           type='assignment',
                           expression="(UMPK_ATP).Km_ADP",
                           model=model_object)
            add_parameter('(UMPK_ADP).Km_AMP_Link', 
                          value=get_reaction_parameters().loc['(UMPK_ADP).Km_AMP', 'value'])
            set_parameters('(UMPK_ADP).Km_AMP_Link',
                           type='assignment',
                           expression="(UMPK_ADP).Km_AMP",
                           model=model_object)
            # update mapping of selected reactions with shared parameters via global quantities
            set_reaction_mapping(reaction='UMPK_ATP',
                         mapping={'Km_AMP': '(UMPK_ADP).Km_AMP_Link'},
                         model=model_object)
            set_reaction_mapping(reaction='UMPK_ADP',
                         mapping={'Km_ATP': '(UMPK_ATP).Km_ATP_Link', 'Km_ADP': '(UMPK_ATP).Km_ADP_Link'},
                         model=model_object)
        if (struct_var.loc['UMPK_ATP']['base'] == 3 and struct_var.loc['UMPK_ADP']['base'] == 0 and struct_var.loc['UMPK_ATPP']['base'] == 1):
            add_parameter('(UMPK_ATP).Km_ATP_Link', 
                          value=get_reaction_parameters().loc['(UMPK_ATP).Km_ATP', 'value'])
            set_parameters('(UMPK_ATP).Km_ATP_Link',
                           type='assignment',
                           expression="(UMPK_ATP).Km_ATP",
                           model=model_object)
            add_parameter('(UMPK_ATP).Km_ADP_Link', 
                          value=get_reaction_parameters().loc['(UMPK_ATP).Km_ADP', 'value'])
            set_parameters('(UMPK_ATP).Km_ADP_Link',
                           type='assignment',
                           expression="(UMPK_ATP).Km_ADP",
                           model=model_object)
            add_parameter('(UMPK_ATPP).Km_ATPP_Link', 
                          value=get_reaction_parameters().loc['(UMPK_ATPP).Km_ATPP', 'value'])
            set_parameters('(UMPK_ATPP).Km_ATPP_Link',
                           type='assignment',
                           expression="(UMPK_ATPP).Km_ATPP",
                           model=model_object)
            # update mapping of selected reactions with shared parameters via global quantities
            set_reaction_mapping(reaction='UMPK_ATP',
                         mapping={'Km_ATPP': '(UMPK_ATPP).Km_ATPP_Link'},
                         model=model_object)
            set_reaction_mapping(reaction='UMPK_ATPP',
                         mapping={'Km_ATP': '(UMPK_ATP).Km_ATP_Link' ,'Km_ADP': '(UMPK_ATP).Km_ADP_Link'},
                         model=model_object)
        if (struct_var.loc['UMPK_ATP']['base'] == 4 and struct_var.loc['UMPK_ADP']['base'] == 2 and struct_var.loc['UMPK_ATPP']['base'] == 2):
            add_parameter('(UMPK_ATP).Km_ATP_Link', 
                          value=get_reaction_parameters().loc['(UMPK_ATP).Km_ATP', 'value'])
            set_parameters('(UMPK_ATP).Km_ATP_Link',
                           type='assignment',
                           expression="(UMPK_ATP).Km_ATP",
                           model=model_object)
            add_parameter('(UMPK_ATP).Km_ADP_Link', 
                          value=get_reaction_parameters().loc['(UMPK_ATP).Km_ADP', 'value'])
            set_parameters('(UMPK_ATP).Km_ADP_Link',
                           type='assignment',
                           expression="(UMPK_ATP).Km_ADP",
                           model=model_object)
            add_parameter('(UMPK_ADP).Km_AMP_Link', 
                          value=get_reaction_parameters().loc['(UMPK_ADP).Km_AMP', 'value'])
            set_parameters('(UMPK_ADP).Km_AMP_Link',
                           type='assignment',
                           expression="(UMPK_ADP).Km_AMP",
                           model=model_object)
            add_parameter('(UMPK_ATPP).Km_ATPP_Link', 
                          value=get_reaction_parameters().loc['(UMPK_ATPP).Km_ATPP', 'value'])
            set_parameters('(UMPK_ATPP).Km_ATPP_Link',
                           type='assignment',
                           expression="(UMPK_ATPP).Km_ATPP",
                           model=model_object)
            # update mapping of selected reactions with shared parameters via global quantities
            set_reaction_mapping(reaction='UMPK_ATP',
                         mapping={'Km_AMP': '(UMPK_ADP).Km_AMP_Link', 'Km_ATPP': '(UMPK_ATPP).Km_ATPP_Link'},
                         model=model_object)
            set_reaction_mapping(reaction='UMPK_ADP',
                         mapping={'Km_ATP': '(UMPK_ATP).Km_ATP_Link', 'Km_ADP': '(UMPK_ATP).Km_ADP_Link', 'Km_ATPP': '(UMPK_ATPP).Km_ATPP_Link'},
                         model=model_object)
            set_reaction_mapping(reaction='UMPK_ATPP',
                         mapping={'Km_ATP': '(UMPK_ATP).Km_ATP_Link', 'Km_ADP': '(UMPK_ATP).Km_ADP_Link', 'Km_AMP': '(UMPK_ADP).Km_AMP_Link'},
                         model=model_object)

    # set all species initial concentrations to either a predefined value or a very low value (1e-12) (as zero can cause problems with certain rate laws)
    for species_name in list(get_species(model=model_object).index):
        if species_name in list(init_conc_dict.keys()):
            set_species(name=species_name,
                        initial_concentration=init_conc_dict[species_name],
                        model=model_object)
        else:
            set_species(name=species_name,
                        initial_concentration=1e-12,
                        model=model_object)
    # return COPASI model object
    return model_object, new_reaction_schemes_dict


# core function to set up and run a parameter estimation and calculate model selection criteria using a COPASI model object
def eval_struct_var(struct_var, exp_data_file_names, exp_data_dataframes, fit_params_info=None, num_replicates=3, 
                    PE_algorithm_info={'name': 'Particle Swarm', 
                                       'settings': {'method': {'Iteration Limit': 800,
                                                               'Swarm Size': 50,
                                                               'Std. Deviation': 1e-06}}},
                    model_data=UDP_GalNAc_model_data, term_libs=UDP_GalNAc_model_term_libs):
    """ Evaluate a given structural variant encoded in a data frame of integers which represents combinations of particular base kinetics and regulatory terms for all reactions.

    First, a COPASI model object is created according to the provided structural variant by modifying a base model (see documentation of the function 'create_model_obj_from_struct_var'). Then, a parameter estimation is set up and computed multiple times using this modified model. Different information criteria (AIC, AICc, BIC, etc.) are calculated and returned. The parameter estimation setup (start values, lower and upper boundaries) is first generated internally and then potentially overwritten with partial or complete information on specific start values as well as lower and upper boundaries derived from a given dictionary (optional argument). The estimation uses the experimental data that is provided as a list of file names.

    The information criteria allow a comparison of models (1) with different numbers of parameters (which change depending on the selected combination of base kinetics and regulatory terms) and (2) which had their parameters estimated using different amounts of measurement points (specifically the number of measurements included in the estimation changes if the 'ADP_Decay' reaction is turned on because then - and only then - are AMP measurements included in the fitting process).

    The function can be used as a fitness function for a directed optimization of the AIC (implemented for example as a pygmo UDP, user-defined problem) and also for other approaches such as brute force (test all possible structural variants) or random sampling (test a random selection of possible structural variants).

    :param struct_var: a structural variant matrix
    :type struct_var: pandas.core.frame.DataFrame
    :param exp_data_file_names: a list of experimental data file names (strings)
    :type exp_data_file_names: list
    :param exp_data_dataframes: a list of pandas DataFrames containing experimental data
    :type exp_data_dataframes: list
    :param fit_params_info: an optional dictionary with start values and boundary information of model parameters to be estimated (if not provided then the start values and boundaries are set according to the internal logic of the functions 'create_model_obj_from_struct_var' and 'eval_struct_var' respectively)
    :type fit_params_info: dictionary or None
    :param num_replicates: the amount of parameter estimation replicates (optional; default: 3)
    :type num_replicates: integer
    :param PE_algorithm_info: an optional dictionary with information on the chosen parameter estimation algorithm (default: Particle Swarm with 300 iterations, swarm size of 50 and standard deviation of 1e-06)
    :type PE_algorithm_info: dictionary
    :param model_data: model-specific function that contains information on the reaction schemes, known literature values for reaction parameters, and placeholder values for all non-zero initial concentrations (default: UDP_GalNAc model)
    :type model_data: function
    :param term_libs: model-specific function that contains all potential base rate law equations as well as all potential regulation terms (default: UDP_GalNAc model)
    :type term_libs: function

    :return PE_results: the list of dictionaries which each contain the evaluation result of a parameter estimation replicate
    :type PE_results: list
    """

    ##TODO: very rarely this error happens: "CFitting (7): Insufficient experimental data (rows requested '8', rows found '-1')." -> in those cases the estimation gets empty data, but why? => add print statements to track the flow of experimental data from the input argument to the add_experiment function and beyond to when the estimation is started to check what's going wrong in those cases

    # ----------------------------------------
    # PREPARE EXP DATA AND CREATE MODEL OBJECT
    # ----------------------------------------
    # create list of experiment names from their file names
    exp_data_names = [file_name.partition('_with')[0] for file_name in exp_data_file_names]
    # create COPASI model object according to the kinetic structure defined in specified structural variant; the function also returns the new reaction schemes
    model_object, new_reaction_schemes_dict = create_model_obj_from_struct_var(struct_var, model_data=model_data, term_libs=term_libs)
    # add experiments to the estimation (model selection criteria values are only comparable if every model variant is estimated with the same number of data points, so even if the variant doesn't produce a given species due to the corresponding reaction either being turned off or missing completely, we still load the associated measurements
    for (name, df) in zip(exp_data_names, exp_data_dataframes):
        add_experiment(name=name, data=df,
                       file_name=name + '_dump_from_add_experiment_function.txt',
                       model=model_object)

    # ----------------------------------
    # PREPARE PARAMETERS TO BE ESTIMATED
    # ----------------------------------
    # get the dictionary of known parameter literature values (first and third return values - reaction schemes and initial concentrations of one of the experiments - are not relevant here)
    _, lit_vals_dict, _ = model_data()
    # flatten the dictionary of literature values (from: 'reaction_name': sub_dict{'parameter_name': value} etc. to: '(reaction_name).parameter_name': value etc.; two level dictionary to one level dictionary
    lit_vals_dict_flattened = dict()
    for key, subdict in lit_vals_dict.items():
        for subdict_key, subdict_val in subdict.items():
            new_key_string = '(' + key + ').' + subdict_key
            lit_vals_dict_flattened.update({new_key_string: subdict_val})
    # get only those reaction parameters that are not mapped to an already existing one, i.e., the set of unique reaction parameters
    reaction_params = get_reaction_parameters(model=model_object)
    reaction_params_unique = reaction_params.loc[reaction_params['mapped_to'] == '']
    ## v23b: change the way parameter boundaries are generated for kcat, K_eq, and Km (make it more realistic and stick closer to a literature value if it is available)
    # generate fit parameters object (list of dictionaries; [{'name': str, 'lower': float, 'upper': float}, ...]); start and boundary values are generated according to several criteria (parameter type, substrate, availability of literature values)
    fit_items = list()
    for i in range(len(reaction_params_unique.index)):
        # set start values and lower and upper boundaries for the parameters that are present in the core model (take literature values into account if available)
        if 'kcat_F' in reaction_params_unique.index[i]:
            # kcat values [1/h] are the enzyme specific turnover number and describe the catalytic turn over of the enzyme; their unit is 1/h so the boundaries are chosen such that they resolve to factors of 10 [1/s]
            if reaction_params_unique.index[i] in lit_vals_dict_flattened.keys():
                # if a literature value is available as start value then use it to generate boundaries
                subdict = {'name': reaction_params_unique.index[i],
                           'start': lit_vals_dict_flattened[reaction_params_unique.index[i]],
                           'lower': lit_vals_dict_flattened[reaction_params_unique.index[i]] * 0.5,
                           'upper': lit_vals_dict_flattened[reaction_params_unique.index[i]] * 1.5}
                fit_items.append(subdict)
            else:
                subdict = {'name': reaction_params_unique.index[i],
                           'start': 3.6e+04,     # 10   [1/s]
                           'lower': 3.6e+02,     # 0.1  [1/s]
                           'upper': 3.6e+06}     # 1000 [1/s]
                fit_items.append(subdict)
        elif 'K_eq' in reaction_params_unique.index[i]:
            # K_eq values [-] are equilibrium constants and describe the equilibrium behavior of the  reaction (>1 -> forward reaction is preferred, <1 backwards reaction is preferred); their boundaries are derived from the start value to generate a comparatively narrow interval
            if reaction_params_unique.index[i] in lit_vals_dict_flattened.keys():
                # if a literature value is available as start value then use it to generate boundaries
                subdict = {'name': reaction_params_unique.index[i],
                           'start': lit_vals_dict_flattened[reaction_params_unique.index[i]],
                           'lower': lit_vals_dict_flattened[reaction_params_unique.index[i]] * 0.5,
                           'upper': lit_vals_dict_flattened[reaction_params_unique.index[i]] * 1.5}
                fit_items.append(subdict)
            else:
                subdict = {'name': reaction_params_unique.index[i],
                           'start': 1,
                           'lower': 1e-03,
                           'upper': 1e+03}
                fit_items.append(subdict)
        elif 'Km_' in reaction_params_unique.index[i]:
            # Km values [mM] are the Michaelis-Menten constants and describe the binding affinity of substrates and products (lower Km = higher affinity)
            if reaction_params_unique.index[i] in lit_vals_dict_flattened.keys():
                # if a literature value is available as start value then use it to generate boundaries
                subdict = {'name': reaction_params_unique.index[i],
                           'start': lit_vals_dict_flattened[reaction_params_unique.index[i]],
                           'lower': lit_vals_dict_flattened[reaction_params_unique.index[i]] * 0.5,
                           'upper': lit_vals_dict_flattened[reaction_params_unique.index[i]] * 1.5}
                fit_items.append(subdict)
            else:
                subdict = {'name': reaction_params_unique.index[i],
                           'start': 1e-01,
                           'lower': 1e-02,
                           'upper': 1e+02}
                fit_items.append(subdict)
        elif 'n_Hill' in reaction_params_unique.index[i]:
            # n_Hill is the Hill exponent [-] which describes the degree of cooperativity of substrate/product binding to the enzyme (n<1 negative cooperativity; n=1 no cooperativity; n>1 positive cooperativity)
            if reaction_params_unique.index[i] in lit_vals_dict_flattened.keys():
                # if a literature value is available as start value then use it to generate boundaries
                subdict = {'name': reaction_params_unique.index[i],
                           'start': lit_vals_dict_flattened[reaction_params_unique.index[i]],
                           'lower': lit_vals_dict_flattened[reaction_params_unique.index[i]] * 0.5,
                           'upper': lit_vals_dict_flattened[reaction_params_unique.index[i]] * 1.5}
                fit_items.append(subdict)
            else:
                subdict = {'name': reaction_params_unique.index[i],
                           'start': 1,
                           'lower': 1,
                           'upper': 5}
                fit_items.append(subdict)
        # no literature values are being used for the estimation of ki, ka, and k1 parameter values
        elif 'ki_' in reaction_params_unique.index[i]:
            # ki values [mM] are inhibitor constants and describe the binding affinity of inhibitors  ("((ki)/(ki+[I]))"; lower ki = higher affinity = higher inhibition); the start values are set such that the inhibition is turned off at the beginning of the parameter estimation
            subdict = {'name': reaction_params_unique.index[i],
                       'start': 1e+04,
                       'lower': 1e-02,
                       'upper': 1e+04}
            fit_items.append(subdict)
        elif 'ka_' in reaction_params_unique.index[i]:
            # ka values [1/mM] are a custom definition of activation constants and describe the reciprocal binding affinity of the activator ("(1+(ka*[A])"; higher ka = higher affinity = higher activation); the start values are set such that the activation is turned off at the beginning of the parameter estimation
            subdict = {'name': reaction_params_unique.index[i],
                       'start': 0,
                       'lower': 0,
                       'upper': 1e+03}
            fit_items.append(subdict)
        elif 'k_MA' in reaction_params_unique.index[i]:
            # k is the forward rate constant of an irreversible mass action rate law ("k*[A]*[B]")
            subdict = {'name': reaction_params_unique.index[i],
                       'start': 1,
                       'lower': 1e-06,
                       'upper': 1e+06}
            fit_items.append(subdict)
        # fall back settings for an unknown parameter
        else:
            subdict = {'name': reaction_params_unique.index[i],
                       'start': 0.1,
                       'lower': 1e-06,
                       'upper': 1e+06}
            fit_items.append(subdict)
    # if a dictionary with fit parameters information was provided as optional argument then overwrite the generated list of fit items with this information where applicable using the parameter names as keys; if this optional argument has not been provided (is None) then just move on
    if fit_params_info is not None:
        for param_name in fit_params_info.keys():
            for fit_item in fit_items:
                if param_name == fit_item['name']:
                    # replace the start value, lower and upper boundaries of that parameter in the fit items
                    fit_item['start'] = fit_params_info[param_name]['start']
                    fit_item['lower'] = fit_params_info[param_name]['lower']
                    fit_item['upper'] = fit_params_info[param_name]['upper']
    # add final list of fit parameters to parameter estimation task
    set_fit_parameters(fit_parameters=fit_items, model=model_object)

    # ----------------------------------------------------------
    # ESTIMATE PARAMETERS AND CALCULATE MODEL SELECTION CRITERIA
    # ----------------------------------------------------------
    # define time course settings (the time course task is called during the parameter estimation)
    set_task_settings(model=model_object,
                      task=T.TIME_COURSE,
                      settings={'problem': {'StepNumber': 48.0,
                                            'StepSize': 0.5,
                                            'Duration': 24.0}})
    # store fit setup before running the estimation
    fit_setup = get_fit_parameters(model=model_object)
    # get rid of the 'cn' column of the fit_setup data frame because it contains objects that can not be pickled (COPASI.CCommonName type)
    fit_setup = fit_setup[fit_setup.columns.drop('cn')]
    # band-aid fix for now: get rid of the false (?) parser error "Function (1): Parser error near character position: '12'."
    model_info.get_copasi_messages(0)
    PE_results = list()
    for i in range(num_replicates):
        try:
            # perform parameter estimation task
            run_parameter_estimation(model=model_object, update_model=False,
                                     method=PE_algorithm_info['name'],
                                     settings=PE_algorithm_info['settings'])
            # save estimation (i.e., fitting) results
            estimated_parameters = get_parameters_solution()
            fit_statistics = get_fit_statistic()
        except:
            # delete all experimental data that was added to the parameter estimation task
            remove_experiments(model=model_object)
            # re-add the experimental data to the model object
            for (name, df) in zip(exp_data_names, exp_data_dataframes):
                add_experiment(name=name, data=df,
                               file_name=name + '_dump_from_add_experiment_function.txt',
                               model=model_object)
            # re-run the parameter estimation
            run_parameter_estimation(model=model_object, update_model=False,
                                     method=PE_algorithm_info['name'],
                                     settings=PE_algorithm_info['settings'])
            # save estimation (i.e., fitting) results of this re-run
            estimated_parameters = get_parameters_solution()
            fit_statistics = get_fit_statistic()
        # sometimes the parameter estimation fails (CFitting(7): Insufficient Data; even though the data 
        # is available) and returns math.inf as objective value (resolves to True with math.isinf()) and 
        # all solution values of all parameters to be estimated are set to nan (a numpy.float64 object; 
        # resolves to True with np.isnan()) -> I don't know why this happens but if it does, it needs to 
        # be caught; the while loop will only move on if the objective value of the estimation is no 
        # longer Inf
        while math.isinf(fit_statistics['obj']):
            ##DEBUG-PRINT
            print('\n CFitting(7) error was raised despite all model data available ... re-add the data to the model object and re-run the parameter estimation.')
            ##DEBUG-PRINT
            # delete all experimental data that was added to the parameter estimation task
            remove_experiments(model=model_object)
            # re-add the experimental data to the model object
            for (name, df) in zip(exp_data_names, exp_data_dataframes):
                add_experiment(name=name, data=df,
                               file_name=name + '_dump_from_add_experiment_function.txt',
                               model=model_object)
            # re-run the parameter estimation
            run_parameter_estimation(model=model_object, update_model=False,
                                     method=PE_algorithm_info['name'],
                                     settings=PE_algorithm_info['settings'])
            # save estimation (i.e., fitting) results of this re-run
            estimated_parameters = get_parameters_solution()
            fit_statistics = get_fit_statistic()
            ##DEBUG-PRINT
            print('Objective value of re-run: ', fit_statistics['obj'])
            ##DEBUG-PRINT
        # calculate different model selection criteria; assumption: measurement errors are independent, 
        # identically and normally distributed with the same variance; if this holds then we can use the 
        # least squares cost function (RSS) to compute AIC, AICc, and BIC
        # (1) Akaike Information Criterion (AIC); equation from Portet2020
        # (2) AICc - AIC corrected for small sample sizes (k > (n/40)); equation from Portet2020
        # (3) Bayesian Information Criterion (BIC); equation from Wikipedia
        ##TODO: find good source besides Wikipedia for BIC equation using RSS
        # for all criteria: lower value = better approximating model
        # RSS: residual sum of squares (the objective value of the least squares parameter estimation)
        # k: the number of estimated model parameters + 1 ('for the variance')
        # n: and the number of data points ('observations') included in the estimation
        # literature references:
        #   AIC and AICc: Portet2020 
        #                 (<https://doi.org/10.1016/j.idm.2019.12.010>)
        #   BIC: Wikipedia 
        #        (<https://en.wikipedia.org/wiki/Bayesian_information_criterion#Gaussian_special_case>)
        RSS = fit_statistics['obj']
        k = len(estimated_parameters.index) + 1
        n = fit_statistics['valid_data_points']
        AIC = n * np.log(RSS / n) + 2 * k
        # it is very unlikely but possible that n-k-1 is 0 in which case AICc is undefined; therefore: only attempt to calculate AICc if that is not the case
        if n-k-1 != 0:
            AICc = n * np.log(RSS / n) + ((2 * k * (k+1)) / (n - k - 1))
        else:
            AICc = np.nan
        BIC = n * np.log(RSS / n) + np.log(n) * k
        # CIC: Custom Information Criteria with much more severe parameter penalties that
        #     1) scales with the 'measurement points times log likelihood' term's absolute value
        #     2) is set to a fixed high value
        #     3) is a scaling factor for the log likelihood term ("error per parameter"); only valid as 
        #        long as the log likelihood term is negative (RSS < n), otherwise the criterion will 
        #        decrease for more parameters k
        CIC1 = n * np.log(RSS / n) + (np.abs((n * np.log(RSS / n))) * 0.01) * k
        CIC2 = n * np.log(RSS / n) + 200 * k
        CIC3 = (n * np.log(RSS / n)) * 1 / k
        # collect all relevant outputs in one result dictionary that will be returned; COPASI object 
        # cannot be pickled, and therefore it is no use to return it (TypeError: cannot pickle
        # 'SwigPyObject' object)
        eval_result = {'structural_variant': struct_var,
                       'reaction_schemes': new_reaction_schemes_dict,
                       'fit_setup': fit_setup,
                       'fit_algorithm': PE_algorithm_info['name'],
                       'fit_algorithm_settings': PE_algorithm_info['settings'],
                       'estimated_parameters': estimated_parameters,
                       'fit_statistics': fit_statistics,
                       'information_criteria': {'AIC': AIC,
                                                'AICc': AICc,
                                                'BIC': BIC,
                                                'CIC1': CIC1,
                                                'CIC2': CIC2,
                                                'CIC3': CIC3}}
        PE_results.append(eval_result)

    # --------
    # CLEAN UP
    # --------
    # remove all newly added variant functions from global function list
    for name in sorted(list(get_reactions().index)):
        remove_function(name + '_Variant')
    # remove the Copasi model object
    remove_datamodel(get_current_model())

    return PE_results


# core function to summarize the evaluation results across all repeats (replicates) and identify the best replicate
def get_overall_fitness_across_all_repl(result_list):
    """ Condense the evaluation information of multiple parameter estimation replicates for the same structural variant (as returned by the eval_struct_var function) into a dictionary which adds the fitness across all replicates (based on the best replicate).

    :param result_list: output of the eval_struct_var function (list of parameter estimation replicate 
    evaluation dictionaries)
    :type result_list: list

    :return log_entry_dict: the dictionary with the overall fitness across all replicates
    :type log_entry_dict: dictionary
    """

    # the structural variant data frame is the same for all replicates so just take it from the first one
    struct_var_df = result_list[0]['structural_variant']
    # construct a dictionary of the calculation results for all replicates
    log_entry_dict = dict()
    log_entry_dict.update({'variant': struct_var_df})
    log_entry_dict.update({'estimation_results': result_list})
    # the evaluation result dictionary contains data from multiple estimations, so we identify the replicate with the best RSS and store all of its MSC values; they get the prefix 'min_' to show that they belong to the replicate with the lowest RSS, i.e., the best replicate (with the minimal error)
    fitness_sub_dict = dict()
    replicates_sorted_by_RSS = sorted(result_list, key=lambda d: d['fit_statistics']['obj'])
    for MSC_key, MSC_value in replicates_sorted_by_RSS[0]['information_criteria'].items():
        fitness_sub_dict.update(
            {'min_' + MSC_key: MSC_value})
    log_entry_dict.update({'fitness': fitness_sub_dict})

    return log_entry_dict


# necessary helper function to get rid of files that are created each time experiments are added to a parameter estimation problem in basico
def clean_up_exp_dump_files(exp_data_file_names):
    """Remove any 'dump' text files that are created by the 'add_experiment' function that is called by the eval_struct_var function.
    
    :param exp_data_file_names: a list of experimental data file names (strings)
    :type exp_data_file_names: list
    """

    exp_data_names = [file_name.partition('_with')[0] for file_name in exp_data_file_names]
    for name in exp_data_names:
        # Windows ('nt') and Unix ('posix') systems use different file path separators 
        if os.name == 'nt':
            full_path = os.getcwd() + '\\' + name + '_dump_from_add_experiment_function.txt'
            # check if the file exists before attempting to delete it
            if Path(full_path).is_file():
                os.remove(full_path)
        elif os.name == 'posix':
            full_path = os.getcwd() + '/' + name + '_dump_from_add_experiment_function.txt'
            # check if the file exists before attempting to delete it
            if Path(full_path).is_file():
                os.remove(full_path)


# core function to create an ensemble of n parameter sets for a given model (structural) variant
def create_parameter_ensemble(struct_var, ensemble_size, exp_data_file_names, exp_data_dataframes, 
                              model_variant_name, fit_params_info=None, PE_algorithm_info={'name': 'Particle Swarm', 
                                                                                           'settings': {'method': {'Iteration Limit': 800,
                                                                                                                   'Swarm Size': 50,
                                                                                                                   'Std. Deviation': 1e-06}}},
                              model_data=UDP_GalNAc_model_data, term_libs=UDP_GalNAc_model_term_libs):
    """ Create an ensemble of kinetic parameter sets from repeatedly fitting the model variant to experimental data (i.e., evaluating its fitness/quality). To obtain parameter sets of a certain quality (i.e., of a certain objective value) first calculate a  first ensemble of the given ensemble size. Then, identify the 30% quantile subset of this  first ensemble. The objective value of the worst parameter set in this subset becomes the threshold. Keep calculating parameter estimations and only keep those with an objective value that is equal to or lower than the threshold until the given ensemble size is reached once again.  Return this second ensemble as a Python object, store it as a pickle file and also create a csv file containing all parameter sets of the second ensemble (that can be used in the batch optimization code).
    
    :param struct_var: a structural variant matrix
    :type struct_var: pandas.core.frame.DataFrame
    :param ensemble_size: number of parameter sets in the ensemble (= number of replicates of the evaluation)
    :type ensemble_size: integer
    :param exp_data_file_names: a list of experimental data file names (strings)
    :type exp_data_file_names: list
    :param exp_data_dataframes: a list of pandas DataFrames containing experimental data
    :type exp_data_dataframes: list
    :param model_variant_name: the name used to create the file names of the output pickle and csv files
    :type model_variant_name: string
    :param fit_params_info: an optional dictionary with start values and boundary information of model parameters to be estimated (if not provided then the start values and boundaries are set according to the internal logic of the functions 'create_model_obj_from_struct_var' and 'eval_struct_var' respectively)
    :type fit_params_info: dictionary or None
    :param PE_algorithm_info: an optional dictionary with information on the chosen parameter estimation algorithm (default: Particle Swarm with 300 iterations, swarm size of 50 and standard deviation of 1e-06)
    :type PE_algorithm_info: dictionary
    :param model_data: model-specific function that contains information on the reaction schemes, known literature values for reaction parameters, and placeholder values for all non-zero initial concentrations (default: UDP_GalNAc model)
    :type model_data: function
    :param term_libs: model-specific function that contains all potential base rate law equations as well as all potential regulation terms (default: UDP_GalNAc model)
    :type term_libs: function

    :return evaluated_vars_log: log of evaluation results (only one element as only one structural variant was evaluated; replicates of the evaluation are the ensemble)
    :type evaluated_vars_log: list
    """

    # ######################## #
    # CALCULATE FIRST ENSEMBLE #
    # ######################## #
    # evaluate the model variant (number of evaluation replicates = size of the ensemble to be created)
    print('Calculating first parameter ensemble ...')
    result_list_ens1 = eval_struct_var(struct_var, exp_data_file_names, exp_data_dataframes,
                                       fit_params_info=fit_params_info, num_replicates=ensemble_size, 
                                       PE_algorithm_info=PE_algorithm_info, model_data=model_data, term_libs=term_libs)
    print('Calculation complete.')
    # gather objective values from all parameter sets of the ensemble in a numpy array
    obj_vals_array = np.array([result_list_ens1[i]['fit_statistics']['obj'] for i in range(len(result_list_ens1))])
    # get the 30% quantile of the first ensemble
    quantile_threshold = np.quantile(obj_vals_array, 0.3)
    # only keep fitting dictionary elements in result list which have an objective value equal to or below the quantile threshold
    result_list_ens1_reduced = list()
    for fitting_dict in result_list_ens1:
        if fitting_dict['fit_statistics']['obj'] <= quantile_threshold:
            result_list_ens1_reduced.append(fitting_dict)

    # ########################################### #
    # EXTEND ENSEMBLE BASED ON QUANTILE THRESHOLD #
    # ########################################### #
    # keep calculating parameter sets until the length of the reduced result list is equal to the given ensemble size
    print('Calculating additional parameter sets ...')
    while len(result_list_ens1_reduced) < ensemble_size:
        result_list_add = eval_struct_var(struct_var, exp_data_file_names, exp_data_dataframes,
                                          fit_params_info=fit_params_info, num_replicates=1, 
                                          PE_algorithm_info=PE_algorithm_info, model_data=model_data, term_libs=term_libs)
        if result_list_add[0]['fit_statistics']['obj'] <= quantile_threshold:
            result_list_ens1_reduced.append(result_list_add[0])
    print('Calculation complete.')

    # ################### #
    # CREATE OUTPUT FILES #
    # ################### #
    # summarize fitness metrics (model selection criteria) across all replicates (i.e., across the ensemble) and store resulting object as pickle file
    evaluated_vars_log= list()
    evaluated_vars_log.append(get_overall_fitness_across_all_repl(result_list_ens1_reduced))
    storage_file = open(model_variant_name+'_'+str(ensemble_size)+'runs_evaluated_vars_log', 'wb')
    pickle.dump(evaluated_vars_log, storage_file)
    storage_file.close()
    # create a data frame of the ensemble (rows: different parameter sets, columns: values of the same parameter)
    list_of_series_to_create_df = list()
    for PE_replicate in result_list_ens1_reduced:
        # get the estimated (sol = solution) value (boundaries are not relevant here)
        estim_params_series = PE_replicate['estimated_parameters'].loc[:, 'sol']
        # add the objective value (residual sum of squares of the fit) to the series
        estim_params_series_extended = pd.concat([estim_params_series, pd.Series(PE_replicate['fit_statistics']['obj'], index=['obj'])], axis=0)
        # collect series in a list that can be used to create a data frame
        list_of_series_to_create_df.append(estim_params_series_extended)
    # concatenate series to create the data frame, then transpose the resulting data frame so that each parameter set is a row
    ensemble_df = pd.concat(list_of_series_to_create_df, axis=1).T
    # create name of the csv file
    csv_file_name = 'sampling_output_' + result_list_ens1_reduced[0]['fit_algorithm'].replace(' ', '_') + str(ensemble_size) + 'runs_' + model_variant_name + '.csv'
    # output the data frame as csv file
    ensemble_df.to_csv(csv_file_name, sep=',')

    return evaluated_vars_log


# helper function to generate combinations of individual variable terms (input for the improved extension search)
def extend_vari_terms_list_with_combis(vari_terms, combi_length):
    """ Extend the list of variable terms containing all individual terms by all unique combinations of a specified length. Return the extended list (individual terms first, then all combinations).
    
    :param vari_terms: collection of variable terms to be added to the model rate laws
    :type vari_terms: list
    :param combi_length: length of the combinations to be generated
    :type combi_length: integer
    """
    # unpack term tuples from their inner lists
    vari_terms_tupls = [tupl_list[0] for tupl_list in vari_terms]
    # get all combinations of length n
    vari_terms_combis = list(itertools.combinations(vari_terms_tupls, combi_length))
    # change elements of the list of combinations: convert tuples of tuples to lists of tuples
    vari_terms_combis_lists = [list(inner_tuple) for inner_tuple in vari_terms_combis]
    # append original vari_terms list to the beginning of the list of combinations and return the result
    return vari_terms + vari_terms_combis_lists


# internal function of the improved_extension_search function used to test the impact of including a variable term on the model fitness
def temp_addition_test(term_info, baseline_var_copy, exp_data_file_names, exp_data_dataframes, fit_params_info,
                       PE_algorithm_info={'name': 'Particle Swarm', 
                                          'settings': {'method': {'Iteration Limit': 800,
                                                                  'Swarm Size': 50,
                                                                  'Std. Deviation': 1e-06}}},
                       model_data=UDP_GalNAc_model_data, term_libs=UDP_GalNAc_model_term_libs):
    """Test the effect of including a variable term in the model. Modify the baseline variant by including a variable term in the associated rate law. Then evaluate the impact of that change by multiple parameter estimation repetitions. The best result is chosen and returned.

    :param term_info: list that can contain one or multiple tuples of three values (row coordinate, column coordinate, and index) defining a variable term or a combination of variable terms to be included temporarily in the current kinetic matrix of the model
    :type term_info: list
    :param baseline_var_copy: deep copy of the current kinetic matrix of the model
    :type baseline_var_copy: pandas.core.frame.DataFrame
    :param exp_data_file_names: a list of experimental data file names (strings)
    :type exp_data_file_names: list
    :param exp_data_dataframes: a list of pandas DataFrames containing experimental data
    :type exp_data_dataframes: list
    :param fit_params_info: an optional dictionary with start values and boundary information of model parameters to be estimated (if not provided then the start values and boundaries are set according to the internal logic of the functions 'create_model_obj_from_struct_var' and 'eval_struct_var' respectively)
    :type fit_params_info: dictionary or None
    :param PE_algorithm_info: an optional dictionary with information on the chosen parameter estimation algorithm (default: Particle Swarm with 300 iterations, swarm size of 50 and standard deviation of 1e-06)
    :type PE_algorithm_info: dictionary
    :param model_data: model-specific function that contains information on the reaction schemes, known literature values for reaction parameters, and placeholder values for all non-zero initial concentrations (default: UDP_GalNAc model)
    :type model_data: function
    :param term_libs: model-specific function that contains all potential base rate law equations as well as all potential regulation terms (default: UDP_GalNAc model)
    :type term_libs: function
    """

    # get all tuples contained in term_info, unpack them to get the coordinates and change the copy of the baseline variant at these locations with the corresponding term values; save the original values so that they can be restored later
    original_vals = list()
    for tupl in term_info:
        # unpack the tuple to get the coordinates of the term and its index
        term_coords_tuple = (tupl[0], tupl[1])
        term_idx = tupl[2]
        # store the original value at the inclusion coordinates
        original_vals.append([term_coords_tuple, baseline_var_copy.iloc[term_coords_tuple]])
        # temporarily add the index of the current term for evaluation
        baseline_var_copy.iloc[term_coords_tuple] = term_idx
    # run multiple parameter estimation repetitions
    extended_baseline_eval_result = eval_struct_var(baseline_var_copy.copy(), 
                                                    exp_data_file_names,
                                                    exp_data_dataframes,
                                                    fit_params_info=fit_params_info,
                                                    num_replicates=5,
                                                    PE_algorithm_info=PE_algorithm_info,
                                                    model_data=model_data,
                                                    term_libs=term_libs)
    # get model parameter values of the best replicate
    replicates_sorted_by_RSS = sorted(extended_baseline_eval_result,
                                      key=lambda d: d['fit_statistics']['obj'])
    best_replicate_estim_params = replicates_sorted_by_RSS[0]['estimated_parameters']
    # restore the baseline variant by setting the original value again
    for entry in original_vals:
        # the first element in each entry is the coordinate tuple, the second element in each entry is the associated original value
        baseline_var_copy.iloc[entry[0]] = entry[1]
    # get fitness across all replicates
    log_entry_dict = get_overall_fitness_across_all_repl(extended_baseline_eval_result)
    # return result dictionary
    return {'term_info': term_info,
            'term_coords': [(tupl[0], tupl[1]) for tupl in term_info],
            'term_index': [(tupl[2]) for tupl in term_info],
            'selected_PE_stage': 1,
            'eval_result': log_entry_dict,
            'best_replicate_estimated_parameters': best_replicate_estim_params,
            'full_PE_output': extended_baseline_eval_result}


# main SSO function that implements the complete model selection algorithm
def improved_extension_search(start_var, vari_terms, exp_data_file_names, exp_data_dataframes,
                              PE_algorithm_info={'name': 'Particle Swarm', 
                                                 'settings': {'method': {'Iteration Limit': 800,
                                                                         'Swarm Size': 50,
                                                                         'Std. Deviation': 1e-06}}},
                              selected_MSC=['AIC', 'AICc', 'BIC', 'CIC1', 'CIC2', 'CIC3'],
                              start_var_ensemble=None, model_data=UDP_GalNAc_model_data, term_libs=UDP_GalNAc_model_term_libs):
    """ Evaluate the performance of the start model variant (or load an existing parameter ensemble of the start model variant) and then successively extend it with variable terms chosen from a list until all terms are added (greedy strategy). The variable terms to be permanently added to the rate laws of the model in each outer iteration of the search loop are selected based on the results of temporary addition tests where each possible variable term (possibly also including all combinations of variable terms of length n; this depends on the contents of the list of variable terms) is temporarily added to an associated model rate law and the resulting change of the model fitness is calculated and recorded (exhaustive approach). The term (or combination of terms) that leads to the highest increase of the models fitness (ranked and compared across all model selection criteria) is chosen for permanent addition. Every time the model is changed via the permanent addition of a variable term the model variant and all associated information are stored in a list. At the end, the best model variants (best mean and median fitness across all different model selection criteria) is selected from that list. An output dictionary that contains the full list of model variants (evaluated_vars_log) and the analysis result featuring the constructed fitness data frame and the identified log entries of model variants with the best mean and the best median fitness is returned.

    :param start var: structural variant matrix of the start model variant
    :type start_var: pandas.core.frame.DataFrame
    :param vari_terms: collection of variable terms to be added to the model rate laws
    :type vari_terms: list
    :param exp_data_file_names: a list of experimental data file names (strings)
    :type exp_data_file_names: list
    :param exp_data_dataframes: a list of pandas DataFrames containing experimental data
    :type exp_data_dataframes: list
    :param PE_algorithm_info: an optional dictionary with information on the chosen parameter estimation algorithm (default: Particle Swarm with 300 iterations, swarm size of 50 and standard deviation of 1e-06)
    :type PE_algorithm_info: dictionary
    :param selected_MSC: list of model selection criteria (any selection of the strings 'AIC', 'AICc', 'BIC', 'CIC1', 'CIC2', 'CIC3'; a list of all of them is used as default)
    :type selected_MSC: list
    :param start var_ensemble: available evaluation log of the start variant
    :type start_var_ensemble: list or None
    :param model_data: model-specific function that contains information on the reaction schemes, known literature values for reaction parameters, and placeholder values for all non-zero initial concentrations (default: UDP_GalNAc model)
    :type model_data: function
    :param term_libs: model-specific function that contains all potential base rate law equations as well as all potential regulation terms (default: UDP_GalNAc model)
    :type term_libs: function

    :return output_dict: contains full list of model variants and the results of the final model selection
    :type output_dict: dictionary
    """

    # ##### #
    # SETUP #
    # ##### #
    # initialize the five logs of the improved extension search loop 
    #           (1) var_df_record: record of permanently changed variants (only the 
    #               data frames) to be used by each next iteration during the run 
    #               time of the loop; the first entry is the unaltered start variant 
    #           (2) best_estim_params_record: record of the best estimated parameters 
    #               (best across all parameter estimation replicates) from the term 
    #               that was permanently added in each iteration; they serve as 
    #               starting values for the estimations in the next iterations
    #           (3) evaluated_vars_log: complete log containing the evaluation 
    #               results of all permanently changed variants
    #           (4) full_temp_addition_test_results: for debugging purposes store full temp. addition 
    #               test results
    #           (5) list_of_permanently_included_terms: each time a variable term is permanently added 
    #               to the model it is added to this list; this way it is possible to check if a variable
    #               term was already included at some point in the model (even if it is not present in 
    #               the current model variant)
    var_df_record = [start_var]
    var_df_record_idx = 0
    best_estim_params_record = []
    best_estim_params_record_idx = 0
    evaluated_vars_log = list()
    full_temp_addition_test_results = list()
    list_of_permanently_included_terms = list()
    # check if an existing evaluation log of the start variant was passed as argument
    if start_var_ensemble is not None:
        # use the existing evaluation log (parameter ensemble of the start variant)
        print('Using existing parameter ensemble of the initial core model (start variant).')
        evaluated_vars_log.append(start_var_ensemble[0])
    else:
        # evaluate start variant/initial core model (negative control) and store results
        print('Evaluating initial core model ...')
        baseline_eval_result = eval_struct_var(var_df_record[var_df_record_idx], 
                                               exp_data_file_names, exp_data_dataframes, 
                                               fit_params_info=None, num_replicates=3, 
                                               PE_algorithm_info=PE_algorithm_info,
                                               model_data=model_data, term_libs=term_libs)
        log_entry_dict = get_overall_fitness_across_all_repl(baseline_eval_result)
        evaluated_vars_log.append(log_entry_dict)
        print('Evaluation of initial core model complete.')
    # store parameter set of best baseline estimation replicate (best RSS) as first 
    # entry in the record of the best estimated parameters
    replicates_sorted_by_RSS = sorted(evaluated_vars_log[0]['estimation_results'], key=lambda d: d['fit_statistics']['obj'])
    best_estim_params_df = replicates_sorted_by_RSS[0]['estimated_parameters']
    best_estim_params_df.columns = best_estim_params_df.columns.str.replace('sol', 'start')
    best_estim_params_record.append(best_estim_params_df.T.to_dict())

    # ############## #
    # EXTENSION LOOP #
    # ############## #
    # generate list of model variants by applying the greedy extension strategy
    for i in tqdm(range(len(vari_terms)), desc='Extension Loop'):
        # set start values to best estimated parameters from previous iteration
        fit_params_info = best_estim_params_record[best_estim_params_record_idx]
        # deep copy of current start (baseline) variant
        baseline_var_copy = var_df_record[var_df_record_idx].copy()

        # TEMPORARY ADDITION TESTS WITH TWO-STAGE PARAMETER ESTIMATIONS
        # temporarily add variable terms one by one and test resulting extended variants; if the term that the coordinates point to has already permanently been added in a previous iteration then it can be found in the log of permanently included terms and it is not necessary to test its temporary addition again; term_info is always a list of either an individual term (one tuple) or of a combination of multiple terms (multiple tuples)
        addition_tests_results = Parallel(n_jobs=-1)(delayed(temp_addition_test)(term_info,
                                                                                 baseline_var_copy,
                                                                                 exp_data_file_names,
                                                                                 exp_data_dataframes,
                                                                                 fit_params_info,
                                                                                 PE_algorithm_info,
                                                                                 model_data,
                                                                                 term_libs)
                                                             for term_info in vari_terms
                                                             if term_info not in list_of_permanently_included_terms)

        # only proceed if there are any temporary addition test results at all (for example it could be that only combinations of terms remain to be tested but all of their individual parts have already been permanently included and so they will all be skipped and no temporary addition tests will be performed leading to an empty list making it impossible to select a term for permanent addition)
        if len(addition_tests_results) == 0:
            # break will terminate the extension loop and skip to the final model selection
            break

        # SELECT TERM FOR PERMANENT ADDITION
        # rank the results of the temporary addition test by each of the different model selection criteria; first: create data frame of addition test results
        list_of_dicts_to_create_df = list()
        for entry in addition_tests_results:
            # the coordinates of the term refer to its position in the structural matrix, the index of the term refers to the type of regulation term that it points to
            entry_term_coords = entry['term_coords']
            entry_term_index = entry['term_index']
            entry_eval_result_fitness = entry['eval_result']['fitness']
            # add coordinates and index of the term that was added as entries to the dictionary
            entry_eval_result_fitness.update({'term_coords': entry_term_coords}) 
            entry_eval_result_fitness.update({'term_index': entry_term_index})
            # append dictionary to overall list
            list_of_dicts_to_create_df.append(entry_eval_result_fitness)
        # convert list of dictionaries into a data frame (dictionary keys become column headers)
        addition_tests_results_df = pd.DataFrame(list_of_dicts_to_create_df)
        # rank all entries in each column from best (lowest) to worst (highest) fitness value (only include the selected model selection criteria columns)
        addition_tests_results_df_ranked = addition_tests_results_df.loc[:, ['min_'+MSC for MSC in selected_MSC]].rank()
        # identify row with the best average rank across all fitness criteria ranks
        addition_tests_results_df['mean_rank'] = addition_tests_results_df_ranked.mean(axis=1)
        row_with_best_mean_rank = addition_tests_results_df[addition_tests_results_df['mean_rank'] == 
                                                           addition_tests_results_df['mean_rank'].min()]
        # store temporary test evaluation result of the regulation term that was selected to be permanently added
        evaluated_vars_log.append(addition_tests_results[row_with_best_mean_rank.index[0]]['eval_result'])
        # store full results of all temporary addition tests
        full_temp_addition_test_results.append(addition_tests_results)

        # PERMANENTLY ADD SELECTED TERM TO CREATE NEW BASELINE VARIANT
        # select coordinates and index of the term that was selected for permanent addition
        add_coords = row_with_best_mean_rank['term_coords'].iloc[0]
        add_index = row_with_best_mean_rank['term_index'].iloc[0]
        # permanently add the selected term index to the copy of the baseline variant
        for tupl_coords, tupl_index_val in zip(add_coords, add_index):
            baseline_var_copy.iloc[tupl_coords] = tupl_index_val
        # add the full list of tuples (vari_terms_list) that was permanently added to a log
        list_of_permanently_included_terms.append(addition_tests_results[row_with_best_mean_rank.index[0]]['term_info'])
        # add the permanently changed copy of the baseline variant to the record so that it is available as starting point (i.e., new baseline variant) for the next iteration
        var_df_record.append(baseline_var_copy)
        var_df_record_idx = var_df_record_idx + 1
        # add the best estimated parameters (across all replicates) of the addition test that was selected for permanent addition; those will be used as the starting values for the next iteration; the optional fit_params_info parameter for the eval_struct_var function needs to be a dictionary that contains start values as well as upper and lower boundaries for all parameters that are to be set as fit items so the data frame of the best estimation result is converted accordingly
        best_estim_params_df = addition_tests_results[row_with_best_mean_rank.index[0]]['best_replicate_estimated_parameters']
        best_estim_params_df.columns = best_estim_params_df.columns.str.replace('sol', 'start')
        best_estim_params_record.append(best_estim_params_df.T.to_dict())
        best_estim_params_record_idx = best_estim_params_record_idx + 1

    # ##################### #
    # FINAL MODEL SELECTION #
    # ##################### #
    # build fitness data frame (rows: model variants, columns: different fitness criteria and index of the model variant in evaluated_vars_log)
    imp_ext_search_fitness_df = pd.DataFrame([var_dict['fitness'] for var_dict in evaluated_vars_log])
    imp_ext_search_fitness_df['var_idx'] = imp_ext_search_fitness_df.index
    # rank all entries in each column from best (lowest) to worst (highest) fitness value (only include the model selection criteria columns)
    imp_ext_search_fitness_df_ranked = imp_ext_search_fitness_df.loc[:, ['min_'+MSC for MSC in selected_MSC]].rank()
    # find row (i.e., model variant) that has the best average rank across all columns (i.e., across all fitness criteria)
    imp_ext_search_fitness_df['mean_rank'] = imp_ext_search_fitness_df_ranked.mean(axis=1)
    row_with_best_mean_rank = imp_ext_search_fitness_df[imp_ext_search_fitness_df['mean_rank'] == 
                                                        imp_ext_search_fitness_df['mean_rank'].min()]
    # find row (i.e., model variant) that has the best median rank across all columns (i.e., across all fitness criteria)
    imp_ext_search_fitness_df['median_rank'] = imp_ext_search_fitness_df_ranked.median(axis=1)
    row_with_best_median_rank = imp_ext_search_fitness_df[imp_ext_search_fitness_df['median_rank'] == 
                                                          imp_ext_search_fitness_df['median_rank'].min()]
    # return a dictionary of the fitness data frame and the identified log entries of the model variants with the best mean and median ranks; as multiple entries can have the same best mean or median rank a map with the getitem method is used to access multiple elements of the log of evaluated variants at once if necessary
    analysis_result = {'fitness_dataframe': imp_ext_search_fitness_df,
                       'fitness_dataframe_ranked': imp_ext_search_fitness_df_ranked,
                       'best_mean_rank': list(map(evaluated_vars_log.__getitem__, list(row_with_best_mean_rank.var_idx))), 
                       'best_median_rank': list(map(evaluated_vars_log.__getitem__, list(row_with_best_median_rank.var_idx)))}
    # create output directory
    output_dict = {'evaluated_vars_log': evaluated_vars_log, 'analysis_result': analysis_result, 
                   'full_temp_addition_test_results': full_temp_addition_test_results}

    return output_dict

# ------------------------------------------------------------------------------------------------------ #

# ANALYSIS FUNCTIONS

def sim_model_obj_with_estim_params(model_obj, estim_params, init_conc_dict, duration=24, step_number=48):
    """ Return the result of a time course simulation for a given model object. Before running the simulation, the kinetic parameter values of the model are set to those of the provided pandas data frame (return object of basico's 'run_parameter_estimation' function) or dictionary. Also, the initial species concentrations are set according to the provided dictionary.

    :param model_obj: a Copasi model object
    :type model_obj: COPASI.CDataModel
    :param estim_params: a data frame of estimated parameters (return object of Copasi's parameter estimation task)
    :type estim_params: pandas.core.frame.DataFrame
    :param init_conc_dict: a dictionary of initial concentrations per metabolite
    :type init_conc_dict: dictionary

    :return tc_result: the result of the time course simulation (return object of Copasi's time course simulation task)
    :type tc_result: pandas.core.frame.DataFrame
    """

    # unpack pandas data frame or dictionary containing the result of a parameter estimation and apply estimated parameter values to model parameters
    if isinstance(estim_params, pd.DataFrame):
        param_names = list(estim_params.index)
        # the name of the column with the estimated parameter values is usually 'sol' but sometimes it can also be 'start' (this happens in algorithms like the improved extension search where a solution of a previous iteration is used as start values for the next iteration and so the column name is changed)
        if 'sol' in estim_params.columns:
            param_estim_vals = list(estim_params.loc[:, 'sol'])
        elif 'start' in estim_params.columns: 
            param_estim_vals = list(estim_params.loc[:, 'start'])
        for i in range(len(param_estim_vals)):
            set_reaction_parameters(name=param_names[i], value=param_estim_vals[i])
    elif isinstance(estim_params, dict):
        for param_name, param_value in estim_params.items():
            set_reaction_parameters(name=param_name, value=param_value)

    # set initial species concentrations
    for species_name, init_conc in init_conc_dict.items():
        set_species(name=species_name, initial_concentration=init_conc)

    # simulate model
    tc_result = run_time_course(model=model_obj, duration=duration, step_number=step_number)

    return tc_result


##TODO: this function is outdated (8/12/2025) ... needs to be updated
def gen_data_from_struct_var(ref_var, ref_eval_result_file_name, exp_data_file_names, exp_data_dataframes, 
                             output_data_file_names, ref_var_ensemble=None, fit_params_info=None, 
                             model_data=UDP_GalNAc_model_data, term_libs=UDP_GalNAc_model_term_libs):
    """ Generates data (with and without added noise) from simulating a reference model variant. The kinetic parameter values required for the simulation are taken either from an existing parameter ensemble of the reference model or from repeated evaluations of the reference model (the best parameter set is chosen based on the RSS).

    :param ref_var: a structural variant matrix that is used as reference to generate data from
    :type ref_var: pandas.core.frame.DataFrame
    :param ref_eval_result_file_name: the name of the file that contains the evaluation result of the reference variant
    :type ref_eval_result_file_name: string
    :param exp_data_file_names: a list of experimental data file names (strings)
    :type exp_data_file_names: list
    :param exp_data_dataframes: a list of pandas DataFrames containing experimental data
    :type exp_data_dataframes: list
    :param output_data_file_names: a list of file name strings for the generated data files
    :type output_data_file_names: list
    :param ref var_ensemble: available evaluation log of the reference variant
    :type ref_var_ensemble: list or None
    :param fit_params_info: an optional dictionary with start values and boundary information of model parameters to be estimated (if not provided then the start values and boundaries are set according to the internal logic of the functions 'create_model_obj_from_struct_var' and 'eval_struct_var' respectively)
    :type fit_params_info: dictionary or None
    :param model_data: model-specific function that contains information on the reaction schemes, known literature values for reaction parameters, and placeholder values for all non-zero initial concentrations (default: UDP_GalNAc model)
    :type model_data: function
    :param term_libs: model-specific function that contains all potential base rate law equations as well as all potential regulation terms (default: UDP_GalNAc model)
    :type term_libs: function
    """

    # evaluate reference model to get different sets of kinetic parameter values from the replicates then pick the best one; check if an existing evaluation log of the reference variant was passed as argument
    if ref_var_ensemble is not None:
        # use the existing evaluation log (parameter ensemble of the reference variant)
        print('Using existing parameter ensemble of the reference model variant.')
        ref_model_eval_result = ref_var_ensemble[0]
    else:
        # evaluate reference model variant and store results
        print('Evaluating reference model variant ...')
        ref_model_eval_result = eval_struct_var(ref_var, exp_data_file_names, exp_data_dataframes, 
                                                fit_params_info=fit_params_info, num_replicates=3, 
                                                PE_algorithm_info=PE_algorithm_info,
                                                model_data=model_data, term_libs=term_libs)
        print('Evaluation of reference model variant complete.')
        storage_file = open(ref_eval_result_file_name, 'wb')
        pickle.dump(ref_model_eval_result, storage_file)
        storage_file.close()
    # sort the repeats of the reference model evaluation according to their RSS (all replicates have the same number of model parameters and are estimated with the same number of experimental data points so they can be compared directly via the residual sum of squares)
    ref_model_eval_result_sorted = sorted(ref_model_eval_result['estimation_results'], key=lambda d: d['fit_statistics']['obj'])
    best_ref_model_params = ref_model_eval_result_sorted[0]['estimated_parameters']
    # create COPASI model object of the reference model variant
    ref_model_obj, _ = create_model_obj_from_struct_var(ref_var, model_data=model_data, term_libs=term_libs)
    # extract initial species concentrations from experimental data
    init_conc_dicts = list()
    for df in exp_data_dataframes:
        # only take those columns that are marked with '_0' to signify initial values
        df_init_cols = df.filter(regex='_0')
        # only the first value at time point 0 is relevant
        init_conc_dict = dict(df_init_cols.iloc[0, :])
        # remove the trailing '_0' from all dictionary keys (subsequent simulations expect species names without them)
        corrected_dict = { k.replace('_0', ''): v for k, v in init_conc_dict.items() }
        init_conc_dicts.append(corrected_dict)
    # run time course simulations (and show progress bar)
    ##TODO: time course data frames use square brackets in column names ...
    ref_model_tcs = []
    for init_conc_dict in tqdm(init_conc_dicts, desc='Ref. Model Time Course Simulations (with best parameter set)'):
        ref_model_tc = sim_model_obj_with_estim_params(ref_model_obj, best_ref_model_params, init_conc_dict)
        ref_model_tcs.append(ref_model_tc)
    # add multiplicative noise (signal * (1 + noise)) to the simulation results; the noise is drawn from a (Gaussian) normal distribution
    ref_model_tcs_noisy = []
    for tc_df in ref_model_tcs:
        mu, sigma = 0, 0.1 # mu: mean; sigma: standard deviation
        gaussian_noise = np.random.normal(mu, sigma, [len(tc_df.index), len(tc_df.columns)])
        ref_model_tcs_noisy.append(tc_df * (1 + gaussian_noise))
    ref_model_tcs_data = [ref_model_tcs, ref_model_tcs_noisy]
    # change time course output data frames to match the format of the experimental data frames (only selected time points and species columns, different column names, take initial substrate concentrations from new time courses and initial enzyme concentrations from experimental data) - repeat this for both normal and noise time course data lists
    output_list = list()
    for df_list in ref_model_tcs_data:
        i = 0
        for init_conc_dict in init_conc_dicts:
            ref_model_tc_reduced = df_list[i].loc[[0, 0.5, 1, 2, 4, 7, 20, 24], :]
            output_df = pd.DataFrame(columns=["Uri", "UMP", "UDP_GalNAc", "UDP", "UTP", "AMP", "ADP",
                                              "ATP", "GalNAc_0", "ATP_0", "ADP_0", "Uri_0", "E_PPA_0", "E_PPK3_0",
                                              "E_UDK_0", "E_NAHK_0", "E_GLMU_0", "E_UMPK_0"],
                                     index=ref_model_tc_reduced.index)
            ##TODO: ... the keys used for ref_model_tc_reduced don't use square brackets ('ADP' instead of '[ADP]' which is the column naming scheme that Copasi is using for time course results) ... why does this work?
            output_df["ADP"] = pd.concat([ref_model_tc_reduced.loc[:, "ADP"]], axis=1)
            output_df["UMP"] = pd.concat([ref_model_tc_reduced.loc[:, "UMP"]], axis=1)
            output_df["UTP"] = pd.concat([ref_model_tc_reduced.loc[:, "UTP"]], axis=1)
            output_df["ATP"] = pd.concat([ref_model_tc_reduced.loc[:, "ATP"]], axis=1)
            output_df["UDP"] = pd.concat([ref_model_tc_reduced.loc[:, "UDP"]], axis=1)
            output_df["UDP_GalNAc"] = pd.concat([ref_model_tc_reduced.loc[:, "UDP_GalNAc"]], axis=1)
            output_df["Uri"] = pd.concat([ref_model_tc_reduced.loc[:, "Uri"]], axis=1)
            output_df["AMP"] = pd.concat([ref_model_tc_reduced.loc[:, "AMP"]], axis=1)
            # add initial species concentrations as first entries in the X_0 columns; GalNAc is not measured, so we assume that is starts with the same concentration as Uri since this is what the experimentalists aim for when setting up the process
            output_df.loc[0, "GalNAc_0"] = ref_model_tc_reduced.loc[0, 'Uri']
            output_df.loc[0, "ATP_0"] = ref_model_tc_reduced.loc[0, 'ATP']
            output_df.loc[0, "ADP_0"] = ref_model_tc_reduced.loc[0, 'ADP']
            output_df.loc[0, "Uri_0"] = ref_model_tc_reduced.loc[0, 'Uri']
            # get initial enzyme concentrations from experimental data
            output_df.loc[0, "E_PPA_0"] = init_conc_dict['E_PPA']
            output_df.loc[0, "E_PPK3_0"] = init_conc_dict['E_PPK3']
            output_df.loc[0, "E_UDK_0"] = init_conc_dict['E_UDK']
            output_df.loc[0, "E_NAHK_0"] = init_conc_dict['E_NAHK']
            output_df.loc[0, "E_GLMU_0"] = init_conc_dict['E_GLMU']
            output_df.loc[0, "E_UMPK_0"] = init_conc_dict['E_UMPK']
            output_list.append(output_df)
            i = i + 1
    # export generated normal (index 0) and noisy (index 1) output data
    for df, file_name in zip(output_list, output_data_file_names):
        df.to_csv(file_name, na_rep='', sep="\t")


def reduce_log_to_unique_entries(evaluated_vars_log):
    """ Go through the evaluation log and merge all duplicate entries (i.e, entries that share the same variant integer matrix). This is necessary since the same variant can be chosen randomly and evaluated at multiple times during the directed optimization process. The resulting reduced log will only contain one entry per variant (with all information from all detected duplicate entries in the original log merged together). The new overall fitness across all duplicates is stored as first entry in a list that is assigned to the key 'fitness'.
    
    :param evaluated_vars_log: log of evaluation results (each of which is a return object from eval_struct_var())
    :type evaluated_vars_log: list

    :return evaluated_vars_log_uniques: new version of the log reduced to only unique entries
    :type evaluated_vars_log_uniques: list
    """

    # get dimensions of the variant data frames that are found in the provided log; as all data frames have the same dimensions in one log it is sufficient to take these dimensions from the first entry
    var_shape_tuple = evaluated_vars_log[0]['variant'].shape
    var_n_rows = var_shape_tuple[0]
    var_n_cols = var_shape_tuple[1]

    # reduce evaluation log to unique variants only; all multiples of each variant will be merged so that only unique entries remain in the log
    # 1) create a version of the evaluation log which only contains string representations of the integer matrices of all variants (reordered so that each row of reg. terms of each variant is monotonically decreasing)
    evaluated_vars_log_sorted_var_mats_str = list()
    for var_dict in evaluated_vars_log:
        var_matrix = var_dict['variant'].values
        # store column of base terms
        var_matrix_base_col = var_dict['variant'].iloc[:, 0].values
        # order reg. term row integers from highest to lowest
        new_var_reg_terms_sub_matrix = list()
        for reg_row in var_matrix[0:var_n_rows, 1:var_n_cols]:
            sorted_reg_row = sorted(reg_row, reverse=True)
            new_var_reg_terms_sub_matrix.append(sorted_reg_row)
        # add base terms column, build complete array and convert to string
        sorted_var_matrix_str = str(np.column_stack((var_matrix_base_col, np.array(new_var_reg_terms_sub_matrix))))
        # store variant matrix string in a list
        evaluated_vars_log_sorted_var_mats_str.append(sorted_var_matrix_str)
    # 2) now that all integer matrices of all variants are sorted according to the same logic and are converted to an object type (here a string) with an unambiguous truth value, they can be compared and the indices of equivalent matrices can be collected
    all_vars_duplicates_indices = list()
    for var_mat_str_obj in evaluated_vars_log_sorted_var_mats_str:
        var_mat_str_obj_duplicate_idxs = list()
        for idx in range(len(evaluated_vars_log_sorted_var_mats_str)):
            if var_mat_str_obj == evaluated_vars_log_sorted_var_mats_str[idx]:
                var_mat_str_obj_duplicate_idxs.append(idx)
        all_vars_duplicates_indices.append(var_mat_str_obj_duplicate_idxs)
    # reduce list of duplicate indices sublists to set of unique entries (because the same sublist of the same duplicates is listed for each variant index that occurs in that sublist)
    all_vars_duplicates_indices_uniques = list()
    for x in map(list, OrderedDict.fromkeys(map(tuple, all_vars_duplicates_indices)).keys()):
        all_vars_duplicates_indices_uniques.append(x)
    # 3) with lists of indices that point to the all duplicates of each variant we can now merge all duplicate entries for each different variant and thereby create a set of unique variant entries
    evaluated_vars_log_uniques = list()
    for sublist in all_vars_duplicates_indices_uniques:
        duplicate_entries = list()
        for idx in sublist:
            duplicate_entries.append(evaluated_vars_log[idx])
        # construct new entry dictionary by merging information from all duplicates
        #   - variant matrix is the same for all
        #   - estimation results are all kept unaltered and stored in a new list
        #   - fitness results are all kept unaltered and stored in a new list
        merged_entry_variant = duplicate_entries[0]['variant']
        merged_entry_estim_res = list()
        merged_entry_fitness = list()
        for entry_dict in duplicate_entries:
            # get estimation results from each parameter estimation replicate
            for PE_replicate in entry_dict['estimation_results']:
                merged_entry_estim_res.append(PE_replicate)
            # get fitness once per duplicate as it was already calculatd across all replicates
            merged_entry_fitness.append(entry_dict['fitness'])
        merged_entry_dict = {'variant': merged_entry_variant,
                             'estimation_results': merged_entry_estim_res,
                             'fitness': merged_entry_fitness}
        # append to list of unique entries
        evaluated_vars_log_uniques.append(merged_entry_dict)
    return evaluated_vars_log_uniques


def search_var_in_log(struct_var, evaluated_vars_log):
    """ Check all variants in a given log of evaluated variants and test if any of them is equivalent to a provided target variant. Return all associated entries (and their indices) where this is true.

    :param struct_var: a structural variant matrix
    :type struct_var: pandas.core.frame.DataFrame
    :param evaluated_vars_log: log of evaluation results (each of which is a return object from eval_struct_var())
    :type evaluated_vars_log: list

    :return identified_entries: list of tuples of the shape (index, entry dictionary)
    :type identified_entries: list
    """

    # unpack variant matrix data frame of the target variant x; store base terms column and rows of regulation terms in separate lists
    x = struct_var
    x_base_terms = list(x.base)
    x_reg_terms = list()
    for i in range(len(x.index)):
        x_reg_terms.append(list(x.iloc[i, 1::]))
    # iterate over all entries of the provided log object and compare the variant matrix data frame of each entry y to the target variant x; sort the reg. term row lists of both to enable a proper comparison of the contents (as the order of the reg. terms is not important); save all identified entries
    identified_entries = list()
    for idx, res_dict in enumerate(evaluated_vars_log):
        # unpack variant matrix data frame of current entry
        y = res_dict['variant']
        # store base terms column and rows of regulation terms in separate arrays
        y_base_terms = list(y.base)
        y_reg_terms = list()
        for i in range(len(y.index)):
            y_reg_terms.append(list(y.iloc[i, 1::]))
        # compare y to x
        if (y_base_terms == x_base_terms and [sorted(y_reg_row) for y_reg_row in y_reg_terms] == [sorted(x_reg_row) for x_reg_row in x_reg_terms]):
            # x and y are equivalent and the corresponding entry and its index are saved
            identified_entries.append((idx, res_dict))
        else:
            # x and y are not equivalent, move on to the next entry
            pass
    return identified_entries


def build_copasi_model_file_from_log_entry(log_entry_dict, output_file_name, model_data=UDP_GalNAc_model_data, term_libs=UDP_GalNAc_model_term_libs):
    """ This function takes all relevant information from a given entry dictionary of a log of evaluated variants and uses it to build a COPASI model (.cps) file. First, a COPASI model object is created according to the kinetic rate laws defined in the structural variant matrix found in the given log entry and all kinetic parameter values are set to the estimated values of the best parameter estimation replicate. Then, the completed model object is exported as a COPASI model file using the given output file name.

    :param log_entry_dict: an entry (dictionary) of a log of evaluated variants (list)
    :type log_entry_dict: dictionary
    :param exp_data_file_names: a list of experimental data file names (strings)
    :type exp_data_file_names: list
    :param exp_data_dataframes: a list of pandas DataFrames containing experimental data
    :type exp_data_dataframes: list
    :param output_file_name: the name of the COPASI model (.cps) file
    :type output_file_name: string
    :param model_data: model-specific function that contains information on the reaction schemes, known literature values for reaction parameters, and placeholder values for all non-zero initial concentrations (default: UDP_GalNAc model)
    :type model_data: function
    :param term_libs: model-specific function that contains all potential base rate law equations as well as all potential regulation terms (default: UDP_GalNAc model)
    :type term_libs: function
    """

    # create initial model object from structural variant matrix found in the log entry dictionary
    model_object, _ = create_model_obj_from_struct_var(log_entry_dict['variant'], model_data=model_data, term_libs=term_libs)
    
    # overwrite all kinetic parameter values with the estimated values of the best replicate
    PE_replicates_sorted_by_obj = sorted(log_entry_dict['estimation_results'],
                                               key=lambda d: d['fit_statistics']['obj'])
    best_estimated_parameters = PE_replicates_sorted_by_obj[0]['estimated_parameters']
    for row in best_estimated_parameters.iterrows():
        if row[0] in get_reaction_parameters().index:
            # the column that contains the estimated parameter values is sometimes called 'sol' and sometimes misleadingly 'start' (depending on the SSO algorithm that was used); just test for both
            if 'sol' in row[1].index:
                set_reaction_parameters(name = row[0], value=row[1].sol)
            elif 'start' in row[1].index:
                set_reaction_parameters(name = row[0], value=row[1].start)

    # export the finished model object as a COPASI (.cps) file
    save_model(output_file_name, model=model_object)


def get_imp_ext_search_trace_for_selected_var(start_var, selected_var, evaluated_vars_log, term_libs=UDP_GalNAc_model_term_libs):
    """ For a log of evaluated variants from an Improved Extension Search run, return the trace which documents the addition of variable terms from a given start to a given selected variant. This trace lists all added terms in the order of their addition with their common names (instead of their term library indices) and also shows which reactions the terms have been added to (again the common names of the reactions are shown instead of their row indices in the variant matrices).

    :param start_var: structural variant matrix of the start variant
    :type start_var: pandas.core.frame.DataFrame
    :param selected_var: structural variant matrix of the selected variant
    :type selected_var: pandas.core.frame.DataFrame
    :param evaluated_vars_log: log of evaluation results specifically of an Improved Extension Search run
    :type evaluated_vars_log: list
    :param term_libs: model-specific function that contains all potential base rate law equations as well as all potential regulation terms
    :type term_libs: function

    :return selected_var_trace:
    :type selected_var_trace:
    """

    # get the log indices of the start and the selected variants
    start_var_idx_in_log = search_var_in_log(start_var, evaluated_vars_log)[0][0]
    selected_var_idx_in_log = search_var_in_log(selected_var, evaluated_vars_log)[0][0]
    # get trace of term additions between start and selected variant; the trace is a list where the order of the elements reflects the order of additions between the start and the selected variant (additions are 
    # given with their common names instead of indices and the reactions these terms are added to are also identified by their common name instead of by their row index in the variant matrix)
    selected_var_trace = list()
    # the loop starts at 'start_variant_idx+1' because we are interested in the first addition which was added to the variant at index 0 to create the variant at index 1 and the loop ends at 'selected_variant_idx+1' because the range() function stops one before the value of its second argument (list(range(0, 3))) => [0, 1, 2]) so in order for the output of the range() function to include the index of the selected variant the range needs to be extended by 1
    for i in range(start_var_idx_in_log+1, selected_var_idx_in_log+1):
        # the information about which variable term was added at which position is part of the sub-dictionary accessed by the 'fitness' key
        term_coords = evaluated_vars_log[i]['fitness']['term_coords']
        term_index = evaluated_vars_log[i]['fitness']['term_index']
        # look up the name of the term in the term library via its index; for the regulation terms look them up in the library of regulation terms; for the ADP_Decay variants (here identified by their coordinates (0, 0) and (1, 0)) set their names accordingly; the second entry of each tuple in the list of term coords is checked to differentiate between the 'base' column and the regulation term columns of the structural variant matrix (if the sum of the column indices of all tuples is zero then the term applies only to the base column)
        if term_coords not in [[(0, 0)], [(1, 0)]] and sum([tpl[1] for tpl in term_coords]) != 0:
            _, regulatory_terms = term_libs()
            ##INFO: there can be multiple terms in the term_coords and term_index lists (as multiple 
            ##      terms can be added at once) - at the moment this is not accounted for here; we just 
            ##      take the first term via [0]
            term_name = regulatory_terms[term_index[0]]['name']
        elif term_coords == [(0, 0)]:
            term_name = 'ADP_Decay.Active.Variant1'
        elif term_coords == [(1, 0)]:
            term_name = 'ADP_Decay.Active.Variant2'
        elif sum([tpl[1] for tpl in term_coords]) == 0:
            # all column indices of the term tuples are 0 so they are pointing to the first column of the structural variant matrix which specifies the base term
            # check the term index; it's a list that contains the IDs of the new base terms; possible combinations: [1,0,0], [2,1,0], [3,0,1], [4,2,2] (NAHK, UDK, UMPK) or [1,1,0], [2,2,1] (PPK3) -> the position of the 1 tells us which new reaction has been turned on -> this is the name that needs to be shown in the trace
            if 1 in term_index:
                if term_index.index(1) == 0:
                    # the first 1 is at position 0 so nothing has changed; [1,0,0] or [1,1,0]
                    term_name = 'No_Change'
                elif term_index.index(1) == 1:
                    # the first 1 is at position 1 which means the ADP variant has been turned on [2,1,0]
                    term_name = 'ADP_Variant'
                elif term_index.index(1) == 2:
                    # the first 1 is at position 2 which means the ATPP variant has been turned on [3,0,1] or [2,2,1]
                    term_name = 'ATPP_Variant'
            else:
                # there's no 1 in the term index, this is only the case if the term index is [4,2,2] in which case both the ADP and the ATPP reaction variants are turned on
                term_name = 'ADP_and_ATPP_Variant'
        # convert the coordinates of the term to a reaction name via the row names of the variant matrix; the 
        # first value of the coordinates tuples specifies the row
            ##INFO: there can be multiple terms in the term_coords and term_index lists (as multiple 
            ##      terms can be added at once) - at the moment this is not accounted for here; we just 
            ##      take the first entry of the first term via [0][0]
        reaction_name = evaluated_vars_log[i]['variant'].index[term_coords[0][0]]
        # append addition information to trace list
        selected_var_trace.append((i, reaction_name, term_name))
    return selected_var_trace


# helper function for generate_fedbatch_model_info(); it is used to change the names of parameters in rate laws strings to their full names (which include the name of the associated reaction)
def replace_substrings(original_string, replacement_dict):
    """ Iterates through the original string and check for matches with the substrings in the replacement_dict. When a match is found, replace the matched substring and skip ahead in the string to avoid overlapping matches. This ensures that each substring is only replaced once and prevents double insertions.

    :param original_string: string with substrings that need to be replaced
    :type original_string: string
    :param replacement_dict: dictionary that maps the new versions (keys) and the old versions (values) of the substrings that need to be replaced
    :type replacement_dict: dictionary

    :return result_string: updated string with all replacements
    :type result_string: string
    """

    # sort the keys by length in descending order to handle longer substrings first
    sorted_keys = sorted(replacement_dict.keys(), key=len, reverse=True)
    # create a list to hold the updated string characters
    result = list()
    # iterate over the original string (character by character)
    i = 0
    n = len(original_string)
    while i < n:
        matched = False
        for old_str in sorted_keys:
            len_old_str = len(old_str)
            if original_string.startswith(old_str, i):
                # replace the matched substring
                result.append(replacement_dict[old_str])
                # the substring was found so we use its length to skip multiple characters at once
                i = i + len_old_str
                matched = True
                break
        # no match is found so we add the current single character to the result and move to the next
        if not matched:
            result.append(original_string[i])
            i = i + 1
    # merge the updated characters of the result list
    result_string = ''.join(result)

    return result_string


"""#TEST
# load an ImpExtSearch result
reading_file = open('rep1\\SSO_ImpExtSearch_v22_rep1_output_dict', 'rb')
output_dict = pickle.load(reading_file)
reading_file.close()

# analyze single output log
analysis_result = output_dict['analysis_result']
# check if there is more than one model variant that achieves the best (lowest) median rank
if len(analysis_result['best_median_rank']) == 1:
    # there's only one model variant that achieves the best median rank so it is saved
    best_ranking_var_log_dict = analysis_result['best_median_rank'][0]
else:
    # there's more than one model variant that achieves the best median rank -> select the model that reaches the best (median+mean) rank
    sum_of_mean_and_median_ranks = analysis_result['fitness_dataframe']['mean_rank']+analysis_result['fitness_dataframe']['median_rank']
    best_ranking_var_log_dict = output_dict['evaluated_vars_log'][sum_of_mean_and_median_ranks.idxmin()]
# best median: v22_rep1_MV5 (median rank 6.0)
# best mean: v22_rep1_MV6 (mean rank 6.3)
# => select MV5 as it reaches the best median rank

model_name = 'ImpExtSearch_v22_rep1_MV5'
inputs = ['F_uri', 'F_gn', 'F_polyp'] # feeding substrates

log_entry_dict = best_ranking_var_log_dict
"""


def generate_fedbatch_model_info(model_name, log_entry_dict, inputs, model_data=UDP_GalNAc_model_data, term_libs=UDP_GalNAc_model_term_libs):
    """ Create a Copasi model object from a structural variant and extract information on its structure (species names, stoichiometry, reaction rate laws, kinetic parameter names). Extend this information to build a fed-batch version of the model which includes defined input feeds, dilution terms for all ODE's and an extra ODE for the change of the volume over time. Export this information in a format that can be used by the dynamic fed-batch optimization code (which requires a separate virtual environment and an older Python version to run HILO-MPC).

    :param model_name: name of the model (will be used to name the output pickle file)
    :type model_name: string
    :param log_entry_dict: an entry (dictionary) of a log of evaluated variants (list)
    :type log_entry_dict: dictionary
    :param inputs: a list with the names of the input flows (['F_x1', 'F_x2', ...] where x1, x2, ... are the names of the feeding substrates/inputs)
    :type inputs: list
    :param model_data: model-specific function that contains information on the reaction schemes, known literature values for reaction parameters, and placeholder values for all non-zero initial concentrations (default: UDP_GalNAc model)
    :type model_data: function
    :param term_libs: model-specific function that contains all potential base rate law equations as well as all potential regulation terms (default: UDP_GalNAc model)
    :type term_libs: function

    :return output_dict: collection of fed-batch model information (names of species, inputs, and parameters; names and equations of rate laws and ODEs)
    :type output_dict: dictionary
    """

    # [Model Object]
    # create initial model object from structural variant matrix found in the log entry dictionary
    model_object, _ = create_model_obj_from_struct_var(log_entry_dict['variant'], model_data=model_data, term_libs=term_libs)

    # [Species]
    # get a list of all species names (x), add an entry for the volume ('V') and sort the list alphabetically
    x_names = list(get_species(model=model_object).index)
    x_names.insert(0, 'V')
    x_names.sort()

    # [Input Feeds]
    # store provided list of the names of the inputs (u) in a new variable
    u_names = inputs

    # [Parameters]
    # get names of kinetic model parameters (p) and clean up the formating of the parameter names; only store names of the 'local' reaction parameters, i.e.,those parameters that are directly mapped to numerical values (and which are being estimated) -> remove parameters that are mapped to global quantities as they only act as pointers to other local parameters
    react_params_df = get_reaction_parameters()
    p_names = [p_name.translate(str.maketrans('', '', '()')).replace('.','_') for p_name in list(react_params_df.index) if react_params_df.loc[p_name, 'type'] == 'local']

    # [Rate Laws]
    # get reaction rate law equations (all rate laws that were constructed from the kinetic information of the structural variant have the suffix '_Variant')
    function_lib = get_functions()
    rate_law_df = function_lib[function_lib.index.str.contains('_Variant')]
    # create dictionary of reaction rate law equations
    rate_law_eq_dict = dict()
    for row in rate_law_df.iterrows():
        # remove the _Variant suffix from the rate law name
        rate_law_name = row[0].replace('_Variant','')
        rate_law_eq_string = row[1].loc['formula']
        # remove all '*1' entries in the rate law equation string to get a reduced version
        rate_law_eq_string_red = rate_law_eq_string.replace('*1','')
        # if an exponent is part of the rate law, change '^' to '**' (the Python syntax for exponents)
        if '^' in rate_law_eq_string_red:
            rate_law_eq_string_red = rate_law_eq_string_red.replace('^','**')
        # get the names of the parameters as they appear in the rate law equation
        param_types = ['k_MA', 'kcat_F', 'Km_', 'K_eq', 'ki_', 'ka_']
        rate_law_params = list()
        for param_type in param_types:
            # split rate law equation string by mathematical operator symbols to get a list of all the variables and parameters of the equation
            list_of_rate_law_eq_parts = re.split('[( ) + - * / ^]', rate_law_eq_string_red)
            # get only the parameters
            rate_law_params.append(list(set([eq_part for eq_part in list_of_rate_law_eq_parts if param_type in eq_part])))
        # flatten list
        rate_law_params_flat = [param for param_sublist in rate_law_params for param in param_sublist]
        # get the full names of the rate law parameters (full names include the name of the reaction; e.g., 'GLMU_kcat_F' instead of just 'kcat_F')
        params_full_names = [p_name for p_name in p_names if f'{rate_law_name}_' in p_name]
        # check if both lists have the same length; if this is not true it is because there are some global parameters represented in the rate law equation string which are not part of the p_names list that is used to get the full parameter names
        if len(rate_law_params_flat) != len(params_full_names):
            # the lengths are different because params_full_names so far only includes the local parameters (that map to numerical values which are estimated) while the rate law equation string also includes references of the global parameters: identify the parameters that are missing in the list of full names
            missing_params = list()
            for param_name in rate_law_params_flat:
                if f'{rate_law_name}_{param_name}' not in params_full_names:
                    missing_params.append(param_name)
            # get a data frame of all global parameters with information which local parameter they are pointing to and get their full names
            global_params_df = react_params_df.loc[react_params_df['type'] == 'global']
            for index, mapping_str in zip(list(global_params_df.index), global_params_df['mapped_to']):
                global_params_df.loc[index, 'mapped_to'] = mapping_str.translate(str.maketrans('','','()')).replace('.','_').replace('_Link', '')
            missing_params_full_names = list()
            for missing_param in missing_params:
                if f'({rate_law_name}).{missing_param}' in global_params_df.index:
                    missing_params_full_names.append(global_params_df.loc[f'({rate_law_name}).{missing_param}', 'mapped_to'])
            # add the full names of the missing parameters to the list of full parameter names
            params_full_names.extend(missing_params_full_names)
        # replace all parameter names with their full names in the rate law equation string; for that a mapping between short and full names of all parameters of the rate law is needed
        # create a dictionary to map each short name (substring) to its full name; the substrings are at the end of the full strings
        params_replacement_dict = {sub: full for sub in rate_law_params_flat for full in params_full_names if full.endswith(sub)}
        rate_law_eq_string_red_fullparamnames = replace_substrings(rate_law_eq_string_red, params_replacement_dict)
        # create dictionary entry
        rate_law_eq_dict.update({rate_law_name: rate_law_eq_string_red_fullparamnames})
    # add the dilution term to the rate law dictionary which is derived from the sum of the specified input flows
    plus_signs = ['+'] * (len(u_names)-1)
    plus_sign_indices = list(range(1, len(u_names)+2, 2)) # plus signs will be added at uneven indices (i.e., between the flow variables)
    for index, obj in zip(plus_sign_indices, plus_signs):
        u_names.insert(index, obj)
    rate_law_eq_dict.update({'D_l': f'1/V*({''.join(u_names)})'})

    # [ODEs]
    # create dictionary of ordinary differential equations (extended with dilution terms); add an ODE for the volume
    active_reaction_names = [f'({react_name})' for react_name in list(rate_law_eq_dict.keys())]
    # get the stoichiometric matrix; the rows only cover non-enzyme species (construct and add the ODE's of the enzymes and the volume later)
    stoich_matrix = get_stoichiometry_matrix(model=model_object)
    # get a subset of the stoichiometric matrix that only includes columns of active reactions; 
    stoich_matrix_subset = stoich_matrix.loc[:,active_reaction_names[0:-1]]
    # iterate over this subset and build the ODEs from it
    ode_dict = dict()
    for row in stoich_matrix_subset.iterrows():
        species_name = row[0].translate(str.maketrans('','','()')) # remove parentheses
        stoich_data = row[1]
        # only get the non-zero entries (i.e., the stoichiometric coefficients that are associated to reactions where that particular species participates)
        stoich_data_nonzero = stoich_data[stoich_data!=0]
        stoich_data_nonzero_react_names = stoich_data_nonzero.index
        # turn stoichiometric coefficients to strings and insert '+' symbols in front of each positive coefficient
        stoich_coeff_strings_list = [f'+{str(coeff)}' if coeff>0 else str(coeff) for coeff in list(stoich_data_nonzero.to_numpy())]
        # combine stoichiometric coefficients with the associated reaction names to a create string that represent the right hand sides of the ODE
        ode_rhs = ''.join([f'{stoich_coeff}*{reaction_name.translate(str.maketrans('','','()'))}' for stoich_coeff, reaction_name in zip(stoich_coeff_strings_list, list(stoich_data_nonzero_react_names))])
        # update the ode dictionary; create the key based on the name of the species (add a 'd' in front of it as that's what HILO-MPC expects for the name of an ODE); also put parentheses around the rhs string and add the dilution term
        ode_dict.update({f'd{species_name}': f'({ode_rhs})-D_l*{species_name}'})
    # add ODE's for the enzymes to the dictionary (species that start with 'E_')
    enzyme_names = [x_name for x_name in x_names if 'E_' in x_name]
    for enzyme_name in enzyme_names:
        # enzymes are not produced or consumed; their concentration is only reduced over time in the fed-batch scenario due to dilution (since the input feeds increase the volume)
        ode_dict.update({f'd{enzyme_name}': f'-D_l*{enzyme_name}'})
    # add an ODE for the volume; the list of inputs was already extended with '+' signs between the input feeds to we can use that here
    ode_dict.update({'dV': ''.join(u_names)})

    # [Export]
    # add extra entries to the list of parameter names that describe the concentrations of the substrate feed inflows (u)
    p_names.extend([f'{u_name}_in' for u_name in u_names if u_name != '+'])
    # gather all of the fed-batch model information and export it using pickle
    output_dict = dict()
    output_dict.update({'model_name': model_name})
    output_dict.update({'species_names': x_names})
    output_dict.update({'input_names':  [u_name for u_name in u_names if u_name != '+']})
    output_dict.update({'parameter_names': p_names})
    output_dict.update({'rate_law_equations': rate_law_eq_dict})
    output_dict.update({'ode_equations': ode_dict})
    # store result as pickle file
    storage_file = open(f'{model_name}.pkl', 'wb')
    pickle.dump(output_dict, storage_file)
    storage_file.close()

    return output_dict


##TODO: build further convenience functions to analyze logs of evaluated variants:
##      1) check all variants in the log for whether any of them includes a provided dictionary of 
##         {reaction: regulation, ...} (variant fragment) and return all associated entries and their 
##         indices in the log

# ------------------------------------------------------------------------------------------------------ #

# PLOTTING FUNCTIONS (for results of the UDP-GalNAc Model)

##TODO: is there a good way to make these function model-independent?

# plot a single model variant with experimental data
def plot_struct_var_tcs(var, exp_data_dataframes, exp_ID, exp_name, plot_name, model_data=UDP_GalNAc_model_data, term_libs=UDP_GalNAc_model_term_libs):
    """ Visualize the performance of a selected structural variant. First, simulate the structural variant for all estimated parameter sets using the initial concentrations of the selected experimental data set. Then, plot the time course trajectories (average +- 95% CI) and add all data points of the selected experimental data sets to the plots.
    
    :param var: log entry of the first variant (top plot)
    :type var: dictionary
    :param exp_data_dataframes: a list of pandas DataFrames containing experimental data
    :type exp_data_dataframes: list
    :param exp_ID: identifier of the selected experiment in the provided list
    :type exp_ID: integer
    :param exp_name: name of the selected experiment in the provided list
    :type exp_name: string
    :param plot_name: string that is used to create the name of the generated plot
    :type plot_name: string
    :param model_data: model-specific function that contains information on the reaction schemes, known literature values for reaction parameters, and placeholder values for all non-zero initial concentrations (default: UDP_GalNAc model)
    :type model_data: function
    :param term_libs: model-specific function that contains all potential base rate law equations as well as all potential regulation terms (default: UDP_GalNAc model)
    :type term_libs: function
    """
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

    # SIMULATION
    # create COPASI model object based on variant 1
    var_model_obj, _ = create_model_obj_from_struct_var(var['variant'], model_data=model_data, term_libs=term_libs)
    # run time course simulation with model object set to all sets of estimated parameters
    var_TC_results_list = []
    PE_replicate_idx = 0
    for PE_replicate in tqdm(var['estimation_results']):
        PE_replicate_tc = sim_model_obj_with_estim_params(var_model_obj, PE_replicate['estimated_parameters'],
                                                          exp_init_conc_dicts[exp_ID])
        # add a 'Time' column (generated from data frame index)
        PE_replicate_tc['Time'] = PE_replicate_tc.index
        # add a 'PE_replicate' column
        PE_replicate_tc['PE_replicate'] = [PE_replicate_idx for i in range(len(list(PE_replicate_tc.index)))]
        # remove 'Values[X]' columns (global quantities that I don't need here)
        PE_replicate_tc.drop(list(PE_replicate_tc.filter(regex='Values')),
                             axis=1, inplace=True)
        # collect simulation data frame
        var_TC_results_list.append(PE_replicate_tc)
        PE_replicate_idx = PE_replicate_idx + 1
    # concatenate all data frames
    var_TC_results_merged_df = pd.concat(var_TC_results_list)
    # define the column order (the column order of the time course data frames returned by basico is not reliable)
    var_TC_results_merged_df = var_TC_results_merged_df.reindex(
        columns=[ "Uri", "UMP", "UDP", "UTP", "UDP_GalNAc", 
                  "AMP", "ADP", "ATP", "ATPP", 
                  "P", "PP", "PolyP", "GalNAc", "GalNAc1P",
                  "Time", "PE_replicate"])

    # MULTIPLOT
    # create seaborn line plots from simulated result data frame (separate plots for uridine-based, GalNAc and GalNAc1p, and adenosine-based species)
    fig, axs = plt.subplots(1, 3, figsize=(16,9), layout='constrained')
    scatterplot_marker_size = 20
    # 1) plots of variant 1 (multiplot row index 0); reshape data frame from wide to long format for plotting of simulated results with error bands
    var_TC_results_list_merged_longform = pd.DataFrame(columns=['PE_replicate', 'Time', 'Species', 'Concentration'])
    var_TC_results_list_merged_longform['PE_replicate'] = pd.concat(
        [pd.Series(var_TC_results_merged_df['PE_replicate']).repeat(len(var_TC_results_merged_df.columns) - 2)],
        axis=0)
    var_TC_results_list_merged_longform = var_TC_results_list_merged_longform.reset_index(drop=True)
    # since the merged data frame in wide format has two extra columns ('Time' and 'PE_replicate'): repeat(n-2), we only want the length of the species columns
    var_TC_results_list_merged_longform['Time'] = pd.concat(
        [pd.Series(var_TC_results_merged_df.index).repeat(len(var_TC_results_merged_df.columns) - 2)], axis=0,
        ignore_index=True)
    var_TC_results_list_merged_longform['Species'] = pd.concat(
        [pd.Series(list(var_TC_results_merged_df.columns[0:-2]) * len(var_TC_results_merged_df.index))], axis=0,
        ignore_index=True)
    var_TC_results_list_merged_longform['Concentration'] = pd.concat([var_TC_results_merged_df.iloc[i, 0:-2]
                                                                       for i in
                                                                       range(len(var_TC_results_merged_df.index))],
                                                                      axis=0, ignore_index=True)
    # plot the mean and 95% confidence interval by aggregating over PE replicates (at each time point)
    sns.lineplot(data=var_TC_results_list_merged_longform[
        var_TC_results_list_merged_longform.Species.isin(['Uri', 'UMP', 'UDP', 'UTP', 'UDP_GalNAc'])],
                 x="Time", y="Concentration", hue="Species", ax=axs[0],
                 palette=[sns.color_palette("Set1")[0], sns.color_palette("Set1")[1], sns.color_palette("Set1")[2],
                          sns.color_palette("Set1")[3], sns.color_palette("Set1")[4]])
    sns.lineplot(data=var_TC_results_list_merged_longform[
        var_TC_results_list_merged_longform.Species.isin(['GalNAc', 'GalNAc1P'])],
                 x="Time", y="Concentration", hue="Species", ax=axs[1],
                 palette=[sns.color_palette("Set1")[0], sns.color_palette("Set1")[1]])
    sns.lineplot(data=var_TC_results_list_merged_longform[
        var_TC_results_list_merged_longform.Species.isin(['AMP', 'ADP', 'ATP', 'ATPP'])],
                 x="Time", y="Concentration", hue="Species", ax=axs[2],
                 palette=[sns.color_palette("Set1")[0], sns.color_palette("Set1")[1], sns.color_palette("Set1")[2],
                          sns.color_palette("Set1")[3]])
    # add experimental data points
    sns.scatterplot(data=exp_data_dataframes[exp_ID], x='Time', y='[Uri]', ax=axs[0],
                    color=sns.color_palette("Set1")[0], s=scatterplot_marker_size)
    sns.scatterplot(data=exp_data_dataframes[exp_ID], x='Time', y='[UMP]', ax=axs[0],
                    color=sns.color_palette("Set1")[1], s=scatterplot_marker_size)
    sns.scatterplot(data=exp_data_dataframes[exp_ID], x='Time', y='[UDP]', ax=axs[0],
                    color=sns.color_palette("Set1")[2], s=scatterplot_marker_size)
    sns.scatterplot(data=exp_data_dataframes[exp_ID], x='Time', y='[UTP]', ax=axs[0],
                    color=sns.color_palette("Set1")[3], s=scatterplot_marker_size)
    sns.scatterplot(data=exp_data_dataframes[exp_ID], x='Time', y='[UDP_GalNAc]', ax=axs[0],
                    color=sns.color_palette("Set1")[4], s=scatterplot_marker_size)
    sns.scatterplot(data=exp_data_dataframes[exp_ID], x='Time', y='[AMP]', ax=axs[2],
                    color=sns.color_palette("Set1")[0], s=scatterplot_marker_size)
    sns.scatterplot(data=exp_data_dataframes[exp_ID], x='Time', y='[ADP]', ax=axs[2],
                    color=sns.color_palette("Set1")[1], s=scatterplot_marker_size)
    sns.scatterplot(data=exp_data_dataframes[exp_ID], x='Time', y='[ATP]', ax=axs[2],
                    color=sns.color_palette("Set1")[2], s=scatterplot_marker_size)
    sns.scatterplot(data=exp_data_dataframes[exp_ID], x='Time', y='[ATPP]', ax=axs[2],
                    color=sns.color_palette("Set1")[3], s=scatterplot_marker_size)
    # overwrite axis labels which were automatically generated from data frame column names
    axs[0].set(xlabel='Time [h]', ylabel='Concentration [mM]')
    axs[1].set(xlabel='Time [h]', ylabel='Concentration [mM]')
    axs[2].set(xlabel='Time [h]', ylabel='Concentration [mM]')
    # save the plot and the associated time course simulation data 
    file_name = plot_name + '_Data' + exp_name
    plt.savefig(file_name + '_tc_plot.png', dpi=200)
    storage_file = open(file_name + '_TC_results.pkl', 'wb')
    pickle.dump(var_TC_results_merged_df, storage_file)
    storage_file.close()


# plot variants in separate rows of a multi plot with experimental data
def compare_struct_var_tcs(var1, var2, exp_data_dataframes, exp_ID, exp_name, plot_name, model_data=UDP_GalNAc_model_data, term_libs=UDP_GalNAc_model_term_libs):
    """ Visualize the performance of two selected structural variants. First, simulate both structural variants for all estimated parameter sets using the initial concentrations of the selected experimental data set. Then, plot the time course trajectories (average +- 95% CI) and add all data points of the selected experimental data sets to the plots.
    
    :param var1: log entry of the first variant (top plot)
    :type var1: dictionary
    :param var2: log entry of the second variant (bottom plot)
    :type var2: dictionary
    :param exp_data_dataframes: a list of pandas DataFrames containing experimental data
    :type exp_data_dataframes: list
    :param exp_ID: identifier of the selected experiment in the provided list
    :type exp_ID: integer
    :param exp_name: name of the selected experiment in the provided list
    :type exp_name: string
    :param plot_name: string that is used to create the name of the generated plot
    :type plot_name: string
    :param model_data: model-specific function that contains information on the reaction schemes, known literature values for reaction parameters, and placeholder values for all non-zero initial concentrations (default: UDP_GalNAc model)
    :type model_data: function
    :param term_libs: model-specific function that contains all potential base rate law equations as well as all potential regulation terms (default: UDP_GalNAc model)
    :type term_libs: function
    """
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

    # -###########################################################################
    # VARIANT 1 SIMULATION
    # create COPASI model object based on variant 1
    var1_model_obj, _ = create_model_obj_from_struct_var(var1['variant'], model_data=model_data, term_libs=term_libs)
    # run time course simulation with model object set to all sets of estimated parameters
    var1_TC_results_list = []
    PE_replicate_idx = 0
    for PE_replicate in tqdm(var1['estimation_results']):
        PE_replicate_tc = sim_model_obj_with_estim_params(var1_model_obj, PE_replicate['estimated_parameters'],
                                                          exp_init_conc_dicts[exp_ID])
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
    # define the column order (the column order of the time course data frames 
    # returned by basico are not reliable)
    var1_TC_results_merged_df = var1_TC_results_merged_df.reindex(
        columns=[ "Uri", "UMP", "UDP", "UTP", "UDP_GalNAc", 
                  "AMP", "ADP", "ATP", "ATPP", 
                  "P", "PP", "PolyP", "GalNAc", "GalNAc1P",
                  "Time", "PE_replicate"])
    # -###########################################################################
    # VARIANT 2 SIMULATION
    # create COPASI model object based on variant 2
    var2_model_obj, _ = create_model_obj_from_struct_var(var2['variant'], model_data=model_data, term_libs=term_libs)
    # run time course simulation with model object set to all sets of estimated parameters
    var2_TC_results_list = []
    PE_replicate_idx = 0
    for PE_replicate in tqdm(var2['estimation_results']):
        PE_replicate_tc = sim_model_obj_with_estim_params(var2_model_obj, PE_replicate['estimated_parameters'],
                                                          exp_init_conc_dicts[exp_ID])
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
    # define the column order (the column order of the time course data frames returned by basico are not reliable)
    var2_TC_results_merged_df = var2_TC_results_merged_df.reindex(
        columns=[ "Uri", "UMP", "UDP", "UTP", "UDP_GalNAc", 
                  "AMP", "ADP", "ATP", "ATPP", 
                  "P", "PP", "PolyP", "GalNAc", "GalNAc1P",
                  "Time", "PE_replicate"])
    # -###########################################################################
    # MULTIPLOT
    # create seaborn line plots from simulated result data frame (separate plots for uridine-based and adenosine-based species) for variant 1 (multi plot row index 0) and variant 2 (multi plot row index 1)
    fig, axs = plt.subplots(2, 3, figsize=(16,9), layout='constrained')
    scatterplot_marker_size = 20
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
    sns.lineplot(data=var1_TC_results_list_merged_longform[
        var1_TC_results_list_merged_longform.Species.isin(['Uri', 'UMP', 'UDP', 'UTP', 'UDP_GalNAc'])],
                 x="Time", y="Concentration", hue="Species", ax=axs[0, 0],
                 palette=[sns.color_palette("Set1")[0], sns.color_palette("Set1")[1], sns.color_palette("Set1")[2],
                          sns.color_palette("Set1")[3], sns.color_palette("Set1")[4]])
    sns.lineplot(data=var1_TC_results_list_merged_longform[
        var1_TC_results_list_merged_longform.Species.isin(['GalNAc', 'GalNAc1P'])],
                 x="Time", y="Concentration", hue="Species", ax=axs[0, 1],
                 palette=[sns.color_palette("Set1")[0], sns.color_palette("Set1")[1]])
    sns.lineplot(data=var1_TC_results_list_merged_longform[
        var1_TC_results_list_merged_longform.Species.isin(['AMP', 'ADP', 'ATP', 'ATPP'])],
                 x="Time", y="Concentration", hue="Species", ax=axs[0, 2],
                 palette=[sns.color_palette("Set1")[0], sns.color_palette("Set1")[1], sns.color_palette("Set1")[2],
                          sns.color_palette("Set1")[3]])
    # add experimental data points
    sns.scatterplot(data=exp_data_dataframes[exp_ID], x='Time', y='[Uri]', ax=axs[0, 0],
                    color=sns.color_palette("Set1")[0], s=scatterplot_marker_size)
    sns.scatterplot(data=exp_data_dataframes[exp_ID], x='Time', y='[UMP]', ax=axs[0, 0],
                    color=sns.color_palette("Set1")[1], s=scatterplot_marker_size)
    sns.scatterplot(data=exp_data_dataframes[exp_ID], x='Time', y='[UDP]', ax=axs[0, 0],
                    color=sns.color_palette("Set1")[2], s=scatterplot_marker_size)
    sns.scatterplot(data=exp_data_dataframes[exp_ID], x='Time', y='[UTP]', ax=axs[0, 0],
                    color=sns.color_palette("Set1")[3], s=scatterplot_marker_size)
    sns.scatterplot(data=exp_data_dataframes[exp_ID], x='Time', y='[UDP_GalNAc]', ax=axs[0, 0],
                    color=sns.color_palette("Set1")[4], s=scatterplot_marker_size)
    sns.scatterplot(data=exp_data_dataframes[exp_ID], x='Time', y='[AMP]', ax=axs[0, 1],
                    color=sns.color_palette("Set1")[0], s=scatterplot_marker_size)
    sns.scatterplot(data=exp_data_dataframes[exp_ID], x='Time', y='[ADP]', ax=axs[0, 1],
                    color=sns.color_palette("Set1")[1], s=scatterplot_marker_size)
    sns.scatterplot(data=exp_data_dataframes[exp_ID], x='Time', y='[ATP]', ax=axs[0, 1],
                    color=sns.color_palette("Set1")[2], s=scatterplot_marker_size)
    sns.scatterplot(data=exp_data_dataframes[exp_ID], x='Time', y='[ATPP]', ax=axs[0, 1],
                    color=sns.color_palette("Set1")[3], s=scatterplot_marker_size)
    # 2) plots of variant 2 (multiplot row index 1); reshape data frame from wide to long format for plotting of simulated results with error bands
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
    sns.lineplot(data=var2_TC_results_list_merged_longform[
        var2_TC_results_list_merged_longform.Species.isin(['Uri', 'UMP', 'UDP', 'UTP', 'UDP_GalNAc'])],
                 x="Time", y="Concentration", hue="Species", ax=axs[1, 0],
                 palette=[sns.color_palette("Set1")[0], sns.color_palette("Set1")[1], sns.color_palette("Set1")[2],
                          sns.color_palette("Set1")[3], sns.color_palette("Set1")[4]])
    sns.lineplot(data=var2_TC_results_list_merged_longform[
        var2_TC_results_list_merged_longform.Species.isin(['GalNAc', 'GalNAc1P'])],
                 x="Time", y="Concentration", hue="Species", ax=axs[1, 1],
                 palette=[sns.color_palette("Set1")[0], sns.color_palette("Set1")[1]])
    sns.lineplot(data=var2_TC_results_list_merged_longform[
        var2_TC_results_list_merged_longform.Species.isin(['AMP', 'ADP', 'ATP', 'ATPP'])],
                 x="Time", y="Concentration", hue="Species", ax=axs[1, 2],
                 palette=[sns.color_palette("Set1")[0], sns.color_palette("Set1")[1], sns.color_palette("Set1")[2],
                          sns.color_palette("Set1")[3]])
    # add experimental data points
    sns.scatterplot(data=exp_data_dataframes[exp_ID], x='Time', y='[Uri]', ax=axs[1, 0],
                    color=sns.color_palette("Set1")[0], s=scatterplot_marker_size)
    sns.scatterplot(data=exp_data_dataframes[exp_ID], x='Time', y='[UMP]', ax=axs[1, 0],
                    color=sns.color_palette("Set1")[1], s=scatterplot_marker_size)
    sns.scatterplot(data=exp_data_dataframes[exp_ID], x='Time', y='[UDP]', ax=axs[1, 0],
                    color=sns.color_palette("Set1")[2], s=scatterplot_marker_size)
    sns.scatterplot(data=exp_data_dataframes[exp_ID], x='Time', y='[UTP]', ax=axs[1, 0],
                    color=sns.color_palette("Set1")[3], s=scatterplot_marker_size)
    sns.scatterplot(data=exp_data_dataframes[exp_ID], x='Time', y='[UDP_GalNAc]', ax=axs[1, 0],
                    color=sns.color_palette("Set1")[4], s=scatterplot_marker_size)
    sns.scatterplot(data=exp_data_dataframes[exp_ID], x='Time', y='[AMP]', ax=axs[1, 1],
                    color=sns.color_palette("Set1")[0], s=scatterplot_marker_size)
    sns.scatterplot(data=exp_data_dataframes[exp_ID], x='Time', y='[ADP]', ax=axs[1, 1],
                    color=sns.color_palette("Set1")[1], s=scatterplot_marker_size)
    sns.scatterplot(data=exp_data_dataframes[exp_ID], x='Time', y='[ATP]', ax=axs[1, 1],
                    color=sns.color_palette("Set1")[2], s=scatterplot_marker_size)
    sns.scatterplot(data=exp_data_dataframes[exp_ID], x='Time', y='[ATPP]', ax=axs[1, 1],
                    color=sns.color_palette("Set1")[3], s=scatterplot_marker_size)
    # overwrite axis labels which were automatically generated from data frame column names
    axs[0, 0].set(xlabel='Time [h]', ylabel='Concentration [mM]')
    axs[0, 1].set(xlabel='Time [h]', ylabel='Concentration [mM]')
    axs[0, 2].set(xlabel='Time [h]', ylabel='Concentration [mM]')
    axs[1, 0].set(xlabel='Time [h]', ylabel='Concentration [mM]')
    axs[1, 1].set(xlabel='Time [h]', ylabel='Concentration [mM]')
    axs[1, 2].set(xlabel='Time [h]', ylabel='Concentration [mM]')
    # save the plot and the associated time course simulation data 
    file_name = plot_name + '_Data' + exp_name
    plt.savefig(file_name + '_Comparison.png', dpi=200)
    list_of_TC_results = [var1_TC_results_merged_df, var2_TC_results_merged_df]
    storage_file = open(file_name + '_TC_results.pkl', 'wb')
    pickle.dump(list_of_TC_results, storage_file)
    storage_file.close()


# plot variants on top of each other in the same plot (differentiated by line style) with experimental data
def compare_struct_var_tcs_merged(var1, var2, exp_data_dataframes, exp_ID, exp_name, plot_name, model_data=UDP_GalNAc_model_data, term_libs=UDP_GalNAc_model_term_libs):
    """ Visualize the performance of two selected structural variants. First, simulate both structural variants for all estimated parameter sets using the initial concentrations of the selected experimental data set. Then, plot the time course trajectories (average +- 95% CI) and add all data points of the selected experimental data sets to the plots.
    
    :param var1: log entry of the first variant (top plot)
    :type var1: dictionary
    :param var2: log entry of the second variant (bottom plot)
    :type var2: dictionary
    :param exp_data_dataframes: a list of pandas DataFrames containing experimental data
    :type exp_data_dataframes: list
    :param exp_ID: identifier of the selected experiment in the provided list
    :type exp_ID: integer
    :param exp_name: name of the selected experiment in the provided list
    :type exp_name: string
    :param plot_name: string that is used to create the name of the generated plot
    :type plot_name: string
    :param model_data: model-specific function that contains information on the reaction schemes, known literature values for reaction parameters, and placeholder values for all non-zero initial concentrations (default: UDP_GalNAc model)
    :type model_data: function
    :param term_libs: model-specific function that contains all potential base rate law equations as well as all potential regulation terms (default: UDP_GalNAc model)
    :type term_libs: function
    """
    # -###########################################################################
    # SIMULATION
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
    # check if the simulation was already calculated and its results saved as pickled file; if not: calculate it now (for both variants)
    try: 
        with open(f'{plot_name}_Data{exp_name}_TC_results_merged.pkl', 'rb') as file:
            list_of_TC_results = pickle.load(file)
            print(f'Simulation results loaded from {plot_name}_Data{exp_name}_TC_results_merged.pkl')
        var1_TC_results_merged_df = list_of_TC_results[0]
        var2_TC_results_merged_df = list_of_TC_results[1]
    except:
        # VARIANT 1 SIMULATION
        # create COPASI model object based on variant 1
        var1_model_obj, _ = create_model_obj_from_struct_var(var1['variant'], model_data=model_data, term_libs=term_libs)
        # run time course simulation with model object set to all sets of estimated parameters
        var1_TC_results_list = []
        PE_replicate_idx = 0
        for PE_replicate in tqdm(var1['estimation_results']):
            PE_replicate_tc = sim_model_obj_with_estim_params(var1_model_obj, PE_replicate['estimated_parameters'],
                                                              exp_init_conc_dicts[exp_ID])
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
        # define the column order (the column order of the time course data frames returned by basico is not reliable)
        var1_TC_results_merged_df = var1_TC_results_merged_df.reindex(
            columns=[ "Uri", "UMP", "UDP", "UTP", "UDP_GalNAc", 
                      "AMP", "ADP", "ATP", "ATPP", 
                      "P", "PP", "PolyP", "GalNAc", "GalNAc1P",
                      "Time", "PE_replicate"])
        # -###########################################################################
        # VARIANT 2 SIMULATION
        # create COPASI model object based on variant 2
        var2_model_obj, _ = create_model_obj_from_struct_var(var2['variant'], model_data=model_data, term_libs=term_libs)
        # run time course simulation with model object set to all sets of estimated parameters
        var2_TC_results_list = []
        PE_replicate_idx = 0
        for PE_replicate in tqdm(var2['estimation_results']):
            PE_replicate_tc = sim_model_obj_with_estim_params(var2_model_obj, PE_replicate['estimated_parameters'],
                                                              exp_init_conc_dicts[exp_ID])
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
        # define the column order (the column order of the time course data frames returned by basico is not reliable)
        var2_TC_results_merged_df = var2_TC_results_merged_df.reindex(
            columns=[ "Uri", "UMP", "UDP", "UTP", "UDP_GalNAc", 
                      "AMP", "ADP", "ATP", "ATPP", 
                      "P", "PP", "PolyP", "GalNAc", "GalNAc1P",
                      "Time", "PE_replicate"])
    # -###########################################################################
    # DATA PREPARATION
    # reshape data frame from wide to long format for plotting of simulated results with error bands
    # 1) results of variant 1
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
    # 2) results of variant 2
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
    # -###########################################################################
    # MULTIPLOT
    # create seaborn line plots from simulated result data frame (separate plots for uridine-based, GalNAc and GalNAc1p, and adenosine-based species) for variant 1 (multi plot row index 0, line style solid) and variant 2 (multi plot row index 0, line style dashed); for the simulation results: plot the mean and 95% confidence interval by aggregating over PE replicates (at each time point)
    fig, axs = plt.subplots(1, 3, figsize=(16,7), layout='constrained')
    scatterplot_marker_size = 20
    # define the species groups for each subplot
    species_groups = [['Uri', 'UMP', 'UDP', 'UTP', 'UDP_GalNAc'],
                      ['GalNAc', 'GalNAc1P'],
                      ['AMP', 'ADP', 'ATP', 'ATPP']]
    # define the colors for each species group
    colors = [[sns.color_palette("Set1")[0], sns.color_palette("Set1")[1], sns.color_palette("Set1")[2], sns.color_palette("Set1")[3], sns.color_palette("Set1")[4]],
              [sns.color_palette("Set1")[0], sns.color_palette("Set1")[1]],
              [sns.color_palette("Set1")[0], sns.color_palette("Set1")[1], sns.color_palette("Set1")[2], sns.color_palette("Set1")[3]]]
    # plot the data
    for i, ax in enumerate(axs):
        # plot variant 1 (dashed lines)
        sns.lineplot(data=var1_TC_results_list_merged_longform[var1_TC_results_list_merged_longform.Species.isin(species_groups[i])],
                     x="Time", y="Concentration", hue="Species", ax=ax, palette=colors[i], linestyle='dashed', legend=False)
        # plot variant 2 (solid lines)
        sns.lineplot(data=var2_TC_results_list_merged_longform[var2_TC_results_list_merged_longform.Species.isin(species_groups[i])],
                     x="Time", y="Concentration", hue="Species", ax=ax, palette=colors[i], linestyle='solid', legend=False)
        # plot experimental data (scatter points)
        for j, species in enumerate(species_groups[i]):
            # check if experimental data is available for this species
            if f'[{species}]' in exp_data_dataframes[exp_ID].columns:
                sns.scatterplot(data=exp_data_dataframes[exp_ID], x='Time', y=f'[{species}]', ax=ax,
                                color=colors[i][j], s=scatterplot_marker_size, label=f'{species}')
        # create custom legend
        legend_elements = list()
        for j, species in enumerate(species_groups[i]):
            # replace underscores with dashes if they show up in the species name 
            if '_' in species:
                species = species.replace('_', '-')
            legend_elements.append(Patch(facecolor=colors[i][j], label=species))
        ax.legend(handles=legend_elements, frameon=True)
    # overwrite axis labels which were automatically generated from data frame column names
    axs[0].set(xlabel='Time [h]', ylabel='Concentration [mM]')
    axs[1].set(xlabel='Time [h]', ylabel='Concentration [mM]')
    axs[2].set(xlabel='Time [h]', ylabel='Concentration [mM]')
    # save the plot and the associated time course simulation data 
    plt.savefig(f'{plot_name}_Data{exp_name}_Comparison_merged.png', dpi=200)
    list_of_TC_results = [var1_TC_results_merged_df, var2_TC_results_merged_df]
    storage_file = open(f'{plot_name}_Data{exp_name}_TC_results_merged.pkl', 'wb')
    pickle.dump(list_of_TC_results, storage_file)
    storage_file.close()


##TODO: update the logic of this function (more and different combinations of indices are possible now)
def compare_occurrence_heatmaps(top_quartile, bottom_quartile, bestX_struct_vars, worstX_struct_vars):
    """ Visualize the occurrence of base and regulation terms of the best and worst variants which are selected from provided top and a bottom quartiles of the fitness distribution.

    :param top_quartile: top quartile cut-off (0% - X%)
    :type top_quartile: float
    :param bottom_quartile: bottom quartile starting point (X% - 100%)
    :type bottom_quartile: float
    :param bestX_struct_vars: top slice of the sorted log of evaluated variants (based on the top quartile cut-off)
    :type bestX_struct_vars: list
    :param worstX_struct_vars: bottom slice of the sorted log of evaluated variants (based on the bottom quartile starting point)
    :type worstX_struct_vars: list
    """

    # -###########################################################################
    # PREPARATION OF DATA FRAMES AND COUNTING OF TERMS
    # initialize occurrence data frames (basis of heat map)
    reg_term_occurrence_df_best = pd.DataFrame(np.zeros([9, 23]))
    reg_term_occurrence_df_best.columns = ['Base_active', '(1)ADP_Act', '(2)ADP_Inhib', '(3)AMP_Act', '(4)AMP_Inhib',
                                           '(5)ATP_Act', '(6)ATP_Inhib', '(7)GalNAc_Act',
                                           '(8)GalNAc_Inhib', '(9)P_Act', '(10)P_Inhib', '(11)PP_Act', '(12)PP_Inhib',
                                           '(13)UDP_Act', '(14)UDP_Inhib', '(15)UDP_GalNAc_Act',
                                           '(16)UDP_GalNAc_Inhib', '(17)UMP_Act', '(18)UMP_Inhib', '(19)Uri_Act',
                                           '(20)Uri_Inhib', '(21)UTP_Act', '(22)UTP_Inhib']
    reg_term_occurrence_df_best.set_index(
        pd.Index(['ADP_Decay', 'GLMU', 'NAHK', 'PPA', 'PPK3_A', 'PPK3_U', 'UDK_UMP', 'UDK_Uri', 'UMPK']), inplace=True)
    reg_term_occurrence_df_worst = pd.DataFrame(np.zeros([9, 23]))
    reg_term_occurrence_df_worst.columns = ['Base_active', '(1)ADP_Act', '(2)ADP_Inhib', '(3)AMP_Act', '(4)AMP_Inhib',
                                            '(5)ATP_Act', '(6)ATP_Inhib', '(7)GalNAc_Act',
                                            '(8)GalNAc_Inhib', '(9)P_Act', '(10)P_Inhib', '(11)PP_Act', '(12)PP_Inhib',
                                            '(13)UDP_Act', '(14)UDP_Inhib', '(15)UDP_GalNAc_Act',
                                            '(16)UDP_GalNAc_Inhib', '(17)UMP_Act', '(18)UMP_Inhib', '(19)Uri_Act',
                                            '(20)Uri_Inhib', '(21)UTP_Act', '(22)UTP_Inhib']
    reg_term_occurrence_df_worst.set_index(
        pd.Index(['ADP_Decay', 'GLMU', 'NAHK', 'PPA', 'PPK3_A', 'PPK3_U', 'UDK_UMP', 'UDK_Uri', 'UMPK']), inplace=True)
    # concatenate all integer arrays of the top best X and worst X structural variants horizontally
    bestX_struct_var_arrays_merged = pd.concat([bestX_struct_vars[i]['variant'] for i in range(len(bestX_struct_vars))],
                                               axis=1)
    worstX_struct_var_arrays_merged = pd.concat(
        [worstX_struct_vars[i]['variant'] for i in range(len(worstX_struct_vars))], axis=1)
    # count occurrence of all reg. terms for each row (= reaction) of the merged array of the regulation term columns all top X and worst X structural variants
    best_df_row_idx = 0
    for idx, row in bestX_struct_var_arrays_merged.iterrows():
        # count occurrence of all base terms; we only care about whether they are on (index 1)
        reg_term_occurrence_df_best.iloc[best_df_row_idx, 0] = list(row).count(1)
        # count occurrence of all reg. terms; 22 possible reg. terms (IDs: 1 to 22)
        for i in range(1, 23):
            # store count directly in occurrence data frame
            reg_term_occurrence_df_best.iloc[best_df_row_idx, i] = list(row.drop(labels=['base'])).count(i)
        best_df_row_idx = best_df_row_idx + 1
    worst_df_row_idx = 0
    for idx, row in worstX_struct_var_arrays_merged.iterrows():
        # count occurrence of all base terms; we only care about whether they are (index 1)
        reg_term_occurrence_df_worst.iloc[worst_df_row_idx, 0] = list(row).count(1)
        # count occurrence of all reg. terms; 22 possible reg. terms (IDs: 1 to 22)
        for i in range(1, 23):
            # store count directly in occurrence data frame
            reg_term_occurrence_df_worst.iloc[worst_df_row_idx, i] = list(row.drop(labels=['base'])).count(i)
        worst_df_row_idx = worst_df_row_idx + 1
    # create data frame of absolute differences (best-worst)
    reg_term_diffs_df = reg_term_occurrence_df_best - reg_term_occurrence_df_worst
    # create data frame of percentage differences (abs(best-worst)/avg(best, worst))
    reg_term_percentage_diffs_df = (np.abs(reg_term_occurrence_df_best - reg_term_occurrence_df_worst) / (
                (reg_term_occurrence_df_best + reg_term_occurrence_df_worst) / 2)) * 100
    # -###########################################################################
    # MULTIPLOT
    # visualize occurrence data frame as seaborn heat map with annotations (color 'k' = black, from 'key')
    fig, axs = plt.subplots(2, 2, layout='constrained')
    bestX_occurrence_hm = sns.heatmap(reg_term_occurrence_df_best, cmap='Reds',
                                      cbar=False, xticklabels=1, yticklabels=1,
                                      vmin=0, vmax=len(bestX_struct_vars),
                                      annot=True, annot_kws={'size': 9},
                                      fmt='.0f', square=True, ax=axs[0, 0])
    worstX_occurrence_hm = sns.heatmap(reg_term_occurrence_df_worst, cmap='Reds',
                                       cbar=False, xticklabels=1, yticklabels=1,
                                       vmin=0, vmax=len(worstX_struct_vars),
                                       annot=True, annot_kws={'size': 9},
                                       fmt='.0f', square=True, ax=axs[0, 1])
    # for differences heatmap the boundaries of the color map need to be set to the highest value that can be found in both bestX and worstX data frames
    diffs_hm_cmap_boundary = max(reg_term_occurrence_df_best.values.max(), reg_term_occurrence_df_worst.values.max())
    diffs_hm = sns.heatmap(reg_term_diffs_df, cmap='vlag',
                           cbar=False, xticklabels=1, yticklabels=1,
                           vmin=-diffs_hm_cmap_boundary, vmax=diffs_hm_cmap_boundary,
                           annot=True, annot_kws={'size': 9},
                           fmt='.0f', square=True, ax=axs[1, 0])
    # percentage differences between 0 and 0 are NaN and therefore not shown in the heat map (we could replace all NaN values in the data frame with .fillna(0))
    percentage_diffs_hm = sns.heatmap(reg_term_percentage_diffs_df, cmap='vlag',
                                      cbar=False, xticklabels=1, yticklabels=1,
                                      vmin=reg_term_percentage_diffs_df.fillna(0).values.min(),
                                      vmax=reg_term_percentage_diffs_df.fillna(0).values.max(),
                                      annot=True, annot_kws={'size': 9, 'color': 'k'},
                                      fmt='.1f', square=True, ax=axs[1, 1])
    # set sub plot titles
    axs[0, 0].set_title(
        'Occurrence of Base and Regulation Terms (Top ' + str(int(top_quartile * 100)) + str("%") + ', ' + str(
            len(bestX_struct_vars)) + ' Variants)')
    axs[0, 1].set_title(
        'Occurrence of Base and Regulation Terms (Worst ' + str(int(np.round(1 - bottom_quartile, 2) * 100)) + str(
            "%") + ', ' + str(len(worstX_struct_vars)) + ' Variants)')
    axs[1, 0].set_title('Differences between Best and Worst Variants')
    axs[1, 1].set_title('Percentage Differences between Best and Worst Variants')
    # set font size of x-axis and y-axis tick labels
    axs[0, 0].set_xticklabels(axs[0, 0].get_xticklabels(), size=7)
    axs[0, 0].set_yticklabels(axs[0, 0].get_yticklabels(), size=7)
    axs[0, 1].set_xticklabels(axs[0, 1].get_xticklabels(), size=7)
    axs[0, 1].set_yticklabels(axs[0, 1].get_yticklabels(), size=7)
    axs[1, 0].set_xticklabels(axs[1, 0].get_xticklabels(), size=7)
    axs[1, 0].set_yticklabels(axs[1, 0].get_yticklabels(), size=7)
    axs[1, 1].set_xticklabels(axs[1, 1].get_xticklabels(), size=7)
    axs[1, 1].set_yticklabels(axs[1, 1].get_yticklabels(), size=7)
    # add shared horizontal color bar
    # mappable = bestX_occurrence_hm.get_children()[0]
    # plt.colorbar(mappable, ax = [axs[0], axs[1]], orientation='horizontal')
    plt.show()
