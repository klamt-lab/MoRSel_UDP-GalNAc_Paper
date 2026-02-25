#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Visualizations of the estimated parameter ensembles of the NegCtrl model variant (first run) and of the best selected model variant rep2_MV5 as multiplots of histograms."""

from func_lib import *
import matplotlib.ticker as tckr
import os

# read the data of the parameter estimation sampling (the first column is all zeros, parameter columns start at .iloc[:,1]); both optimization variants (a and b) use the same parameter ensembles so it is sufficient to load them from one of the two directories
neg_ctrl_fits_data = pd.read_csv('sampling_output_Particle_Swarm50runs_v23bNegControl_Exp36.37.39.csv')
rep2_MV5_fits_data = pd.read_csv(f'{os.getcwd()}\\rep2\\sampling_output_Particle_Swarm50runs_ImpExtSearch_v23b_rep2_ModelVar5.csv')
# change the 'obj' column name; the objective value of a least squares estimation is the residual sum of squares (RSS)
neg_ctrl_fits_data = neg_ctrl_fits_data.rename(columns={'obj': 'RSS'})
rep2_MV5_fits_data = rep2_MV5_fits_data.rename(columns={'obj': 'RSS'})

# get literature values that were used for the parameter estimations
_, param_dict, _ = UDP_GalNAc_model_data()

def hist_multiplot(fits_data, param_dict, num_rows, num_cols, file_name):
    """ Visualize the data of a parameter ensemble as a multiplot of histograms.
    
    :param fits_data: data of the parameter ensemble
    :type fits_data: pandas.core.frame.DataFrame
    :param param_dict: information on literature values that were used as start values for the estimations
    :type param_dict: dictionary
    :param num_rows: number of rows of the multiplot
    :type num_rows: integer
    :param num_cols: number of columns of the multiplot
    :type num_cols: integer
    :param file_name: name used for the created image files
    :type file_name: string
    """

    # set multiplot size and define layout of the subplots
    fig, ax = plt.subplots(num_rows, num_cols, constrained_layout=True, figsize=(num_cols*3, num_rows*3), dpi=200)
    # create histograms
    n = 0
    for i in range(num_rows):
        for j in range(num_cols):
            # len-2 because two cols in fits_data are no param cols (first and last)
            if n <= len(fits_data.columns)-2:
                # get the full name of the parameter column (contains reaction and parameter name; pattern: "(reaction_name).parameter_name")
                name = fits_data.columns[n+1]
                # for all kinetic model parameters, isolate the names of the associated reaction and the parameter itself from the full name
                if name != 'RSS':
                    name_assoc_react = name.split('.')[0].split('(')[1].split(')')[0]
                    name_param = name.split('.')[1]
                else:
                    name_assoc_react = ''
                    name_param = ''
                # figure out if the histogram needs a log x-axis; if the maximum value is greater equal than 1e+02 times the minimum value than the range is too large for normal plotting; exclude the first column because it is the index column
                vals = list(fits_data.iloc[:,n+1])
                if max(vals) >= min(vals)*1e2:
                    log_flag = True
                else:
                    log_flag = False
                # plot estimated parameter values as histogram (kde = kernel density estimate, provides information about the shape of the distribution)
                sns.histplot(data=fits_data, x=name, ax=ax[i,j],
                             kde=False, log_scale=log_flag)
                # draw vertical line for literature value (if one was available)
                if name_assoc_react in param_dict.keys() and name_param in param_dict[name_assoc_react].keys():
                    ax[i,j].axvline(x=param_dict[name_assoc_react][name_param], ymin=0, ymax=len(fits_data.index),
                                    color='r', alpha=0.75, linestyle='-')
                # remove axis labels
                ax[i,j].set(xlabel='',ylabel='')
                # change size of x axis and y axis tick labels
                ax[i,j].tick_params(axis='x', labelsize=12)
                ax[i,j].tick_params(axis='y', labelsize=12)
                # set parameter name as title
                ax[i,j].set_title(name, fontsize=16)
                n = n+1
            else:
                # only use the necessary subplots
                ax[i, j].axis('off')

    # save plot as .pdf
    plt.savefig(f'{file_name}.pdf')
    # save plot as .png (for Word)
    plt.savefig(f'{file_name}.png')
    # save plot as .svg (for Word)
    plt.savefig(f'{file_name}.svg')

# create plots
neg_ctrl_file_name = 'v23bNegCtrl_PartSwarm50runs_hist_fig'
rep2_MV5_file_name = 'v23brep2MV5_PartSwarm50runs_hist_fig'
hist_multiplot(neg_ctrl_fits_data, param_dict, 5, 9, neg_ctrl_file_name)
hist_multiplot(rep2_MV5_fits_data, param_dict, 6, 9, rep2_MV5_file_name)
