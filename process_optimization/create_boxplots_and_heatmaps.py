#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import pickle
import seaborn as sns


# define plotting functions

def visualize_opt_results(result_file_name, opt_variant_name):
    """Plot the optimization results (i.e., the predicted optimal initial enzyme and substrate concentrations for all parameter sets) as combination of box and strip plot. "The box shows the quartiles of the dataset while the whiskers extend to show the rest of the distribution, except for points that are determined to be “outliers” using a method that is a function of the inter-quartile range." (taken from https://seaborn.pydata.org/generated/seaborn.boxplot.html).

    :param result_file_name: path of the calculated optimization result (.pkl file)
    :type fits_data_path: string
    :param opt_variant_name: optimization variant ('a' or 'b')
    :type opt_variant_name: string
    """

    # LOAD AND PREPARE DATA
    # import an opt_results list object via pickle
    with open(f'{result_file_name}.pkl', 'rb') as f:
        opt_results = pickle.load(f)
    # convert enzyme concentrations [mM] of each optimization result to [g/l] with molecular weights: (mM value/1000) [mol/l] * molecular weight [g/mol] = lab value [g/l]
    molecular_weights = pd.Series({'E_GLMU': 49190, 'E_NAHK': 39903, 'E_PPA': 19313, 'E_PPK3': 34740, 'E_UDK': 24353, 'E_UMPK': 22405})
    for opt_result in opt_results:
        # unpack data frame
        opt_result_df = opt_result[0]
        # iterate over the enzyme rows
        for enzyme in molecular_weights.index:
            if enzyme in opt_result_df.index:
                # convert the 'start' and 'sol' enzyme concentration values from mM to g/l
                opt_result_df.loc[enzyme, 'start'] = (opt_result_df.loc[enzyme, 'start'] / 1000) * molecular_weights[enzyme]
                opt_result_df.loc[enzyme, 'sol'] = (opt_result_df.loc[enzyme, 'sol'] / 1000) * molecular_weights[enzyme]
    # prepare opt result data for plotting -> seaborn plots need pandas data frame:
    # cols: ['name', 'sol']
    # rows: amount of rows = amount of different species (enzymes and substrates) * amount of optimizations
    # method: create pandas data frame from pandas series; then use data frame for plotting
    plot_data = []
    for sublist in opt_results:
        # unpack sublist
        result_df = sublist[0]    # pandas data frame
        result_stats = sublist[1] # dictionary
        # get series of recommended initial concentrations for all enzymes and substrates
        rec_init_concs = result_df.iloc[:,3]
        # append series to plot_data list
        plot_data.append(rec_init_concs)
    # stack all series vertically
    plot_data = pd.concat(plot_data, axis = 0)
    # reset index to create data frame (former index is now names column)
    plot_df = plot_data.reset_index()
    # check which optimization variables where used (three groups: all enzymes, main substrates Uri and GalNAc, and ATP) and define necessary variables and specific options for plotting
    plot_info_dict = dict()
    fig_num_cols = 0
    fig_width_ratios = list()
    if set(['E_GLMU', 'E_NAHK', 'E_PPA', 'E_PPK3', 'E_UDK', 'E_UMPK']).issubset(list(opt_results[0][0].index)):
        # create data frames with enzyme concentrations used as starting points for the optimization; the start concentrations were the same for each optimization run so we can take them from the first element of the opt_results list
        pre_opt_df = pd.DataFrame({'name': ['E_GLMU', 'E_NAHK', 'E_PPA', 'E_PPK3', 'E_UDK', 'E_UMPK'],
                                   'start': opt_results[0][0].loc[['E_GLMU', 'E_NAHK', 'E_PPA', 'E_PPK3', 'E_UDK', 'E_UMPK'], 'start']}).reset_index(drop=True)
        # define specific plotting info
        plot_info_dict.update({'enzymes_subplot': {'species_names': ['E_GLMU', 'E_NAHK', 'E_PPA', 'E_PPK3', 'E_UDK', 'E_UMPK'],
                                                   'pre_opt_conc_df': pre_opt_df,
                                                   'boxplot_width': 0.4,
                                                   'x_axis_tick_labels': ['GLMU', 'NAHK', 'PPA', 'PPK3', 'UDK', 'UMPK'],
                                                   'y_axis_label': 'Concentration [g/l]'}})
        fig_num_cols += 1
        fig_width_ratios.append(len(['E_GLMU', 'E_NAHK', 'E_PPA', 'E_PPK3', 'E_UDK', 'E_UMPK']))
    if set(['Uri', 'GalNAc']).issubset(list(opt_results[0][0].index)):
        # create data frames with main substrate concentrations used as starting points for the optimization; the start concentrations were the same for each optimization run so we can take them from the first element of the opt_results list
        pre_opt_df = pd.DataFrame({'name': ['Uri', 'GalNAc'],
                                   'start': opt_results[0][0].loc[['Uri', 'GalNAc'], 'start']}).reset_index(drop=True)
        # define specific plotting info
        plot_info_dict.update({'mainS_subplot': {'species_names': ['Uri', 'GalNAc'],
                                                 'pre_opt_conc_df': pre_opt_df,
                                                 'boxplot_width': 0.4,
                                                 'x_axis_tick_labels': ['Uridine', 'GalNAc'],
                                                 'y_axis_label': 'Concentration [mM]'}})
        fig_num_cols += 1
        fig_width_ratios.append(len(['Uri', 'GalNAc']))
    if 'ATP' in list(opt_results[0][0].index):
        # create data frames with ATP concentration used as starting point for the optimization; the start concentrations were the same for each optimization run so we can take them from the first element of the opt_results list
        pre_opt_df = pd.DataFrame({'name': ['ATP'],
                                   'start': opt_results[0][0].loc[['ATP'], 'start']}).reset_index(drop=True)
        # define specific plotting info
        plot_info_dict.update({'ATP_subplot': {'species_names': ['ATP'], 
                                               'pre_opt_conc_df': pre_opt_df,
                                               'boxplot_width': 0.4,
                                               'x_axis_tick_labels': ['ATP'],
                                               'y_axis_label': 'Concentration [mM]'}})
        fig_num_cols += 1
        fig_width_ratios.append(len(['ATP']))

    # PLOTTING
    # visualize the optimization result with a box plot+ strip plot: recommended initial enzyme ans substrate concentration values are used to create box plots + are drawn as scattered dots per species category; the number of subplot columns depends on which optimization variables were used
    fig, axs = plt.subplots(1, fig_num_cols, gridspec_kw={'width_ratios': fig_width_ratios}, figsize=(8, 4))
    # adjust spacing between subplots
    plt.subplots_adjust(wspace=0.3)
    # create each subplot
    for ax, (subplot_name, subplot_info) in zip(axs, plot_info_dict.items()):
        # unpack some plot information
        names = subplot_info['species_names']
        box_width = subplot_info['boxplot_width']
        # 1) Box Plot
        # create boxplot using the matplotlib default function because it allows defining the width of each box plot; by default matplotlib uses 1-based indexing for x-ticks while seaborn uses 0-based indexing => therefore: fix box plot positions manually; also: make the widths of the box plots depend on how many are drawn together in  subplot (fewer box plots -> lower width)
        positions = np.arange(len(names))
        ax.boxplot([plot_df.loc[plot_df['name'] == name, 'sol'].values for name in names],
                    showfliers=False,
                    widths=box_width,
                    patch_artist=True,
                    positions=positions)
        # change transparency of boxes
        for patch in ax.patches:
            r, g, b, a = patch.get_facecolor()
            patch.set_facecolor((r, g, b, .3))
        # 2) Strip Plot
        # create seaborn strip plot (with custom x-positions); to align with box plots, repeat the x-positions for each point
        stripplot_data = plot_df.loc[plot_df['name'].isin(names)].copy()
        stripplot_data['x_pos'] = stripplot_data['name'].apply(lambda x: positions[list(names).index(x)])
        sns.stripplot(x='x_pos',  # use custom x-positions
                      y='sol',
                      data=stripplot_data,
                      s=2,
                      jitter=True,  # spread points within the box width
                      dodge=False,  # avoid shifting points
                      ax=ax)
        # 3) Scatter Plot
        # create seaborn scatter plot to add markers for start concentrations (with custom x-positions)
        start_df = subplot_info['pre_opt_conc_df'].copy()
        start_df['x_pos'] = start_df['name'].apply(lambda x: positions[list(names).index(x)])
        sns.scatterplot(x='x_pos',
                        y='start',
                        data=start_df,
                        color='r',
                        s=1000,
                        marker='_',
                        ax=ax)
        # set custom x-axis tick positions and labels
        ax.set_xticks(range(len(names)))
        ax.set_xticklabels(subplot_info['x_axis_tick_labels'])
        # add padding to the x-axis limits
        ax.set_xlim(-0.5, len(names) - 0.5)
        # set axis labels
        ax.set(xlabel='', ylabel=f'{subplot_info["y_axis_label"]} \n')
    # save complete figure
    fig.tight_layout()
    fig.savefig(f'{result_file_name}_{opt_variant_name}_boxplots.png', bbox_inches='tight', dpi=200)

def visualize_OP01_cross_val_tables(cross_val_result, scoring_result, file_name):
    """Plot the results of the cross-validation of OP01 optimization results (titer heatmap and score bar plot). The rows of the heatmap correspond to the parameter sets of the ensemble and the columns correspond to the different optimization results (sets of initial enzyme and substrate concentrations). The value of each cell is then derived from simulating the model with a particular combination of a parameter set and an optimization result. The central diagonal corresponds to cases where the parameter set was the one used to calculate the optimization in the first place whereas the off-diagonal elements are cases where an optimization result is combined with a new parameter set that was not used to solve it. The scores are calculated as described in the relevant scoring function. The best scoring optimization result is highlighted in both plots.

    :param cross_val_result: cross-validation result table (titers)
    :type cross_val_result: list
    :param scoring_result: scoring results of all optimization results
    :type scoring_result: pandas.core.frame.DataFrame
    :param file_name: name of the generated image file
    :type file_name: string
    """

    # get index of best scoring optimization result
    best_total_score_index_name = scoring_result.loc[:,'total_score'].idxmax()
    best_total_score_index = list(scoring_result.index).index(best_total_score_index_name)

    # setup multiplot
    fig, ((ax1, cbar_ax), (ax2, dummy_ax)) = plt.subplots(nrows=2, ncols=2, figsize=(7, 6), sharex='col',
                                                          gridspec_kw={'height_ratios': [10, 1], 'width_ratios': [20, 1]})

    # HEATMAP
    hm = sns.heatmap(pd.DataFrame(cross_val_result),
                     cmap=sns.cm.rocket_r,
                     xticklabels=False,
                     yticklabels=True,
                     cbar_ax=cbar_ax,
                     ax=ax1)
    hm.tick_params(left=True, bottom=False)
    hm.set_ylabel('Parameter Sets')
    cbar_ax.set_ylabel('\n Concentration [mM]')
    # set y-axis tick labels: only display labels for multiples of 5 and the first parameter set
    num_rows = len(cross_val_result)
    yticks_pos = []
    yticks_labels = []
    for i in range(num_rows):
        if (i + 1) % 5 == 0 or (i + 1) == 1:
            yticks_pos.append(i + ((i+1-i)*0.5))  # place tick mark in the middle of the row
            yticks_labels.append(f'P{i + 1}')
    # set the y-ticks and labels
    ax1.set_yticks(yticks_pos)
    ax1.set_yticklabels(yticks_labels)
    # highlight best scoring column of the heatmap
    x_start = best_total_score_index
    x_end = best_total_score_index + 1
    y_start = 0
    y_end = len(cross_val_result)
    rect = plt.Rectangle((x_start, y_start),
                         x_end - x_start,
                         y_end - y_start,
                         linewidth=2,
                         edgecolor='red',
                         facecolor='none',
                         transform=ax1.transData)
    ax1.add_patch(rect)

    # BAR PLOT
    # create the bar plot
    bars = ax2.bar([i + 0.5 for i in range(len(scoring_result))],
                   scoring_result['total_score'],
                   align='center')
    # enable tick marks on both axis
    ax2.tick_params(axis='both', which='both', length=7, direction='out', bottom=True, left=True)
    # set the x-ticks and labels to only display x tick labels that are multiples of 5 (and the one of the very first opt. result)
    x_tick_labels = list()
    x_tick_pos = list()
    for i, opt_name in enumerate(scoring_result.index):
        opt_num = int(opt_name[1::])
        if opt_num % 5 == 0 or opt_num == 1:
            x_tick_labels.append(opt_name)
            x_tick_pos.append(i + 0.5)  # +0.5 to center the tick under the bar
    ax2.set_xticks(x_tick_pos)
    ax2.set_xticklabels(x_tick_labels)
    # set y-axis ticks and labels to show only 0 and 3 (lowest and highest possible scoring values)
    ax2.set_yticks([0, 3])
    ax2.set_yticklabels(['0', '3'])
    # set axis labels
    ax2.set_xlabel('Optimization Results')
    ax2.set_ylabel('Score [-] \n')
    # highlight the bar of the best scoring optimization result
    bars[best_total_score_index].set_color('red')
    # remove interior grid lines
    ax2.grid(False, axis='both')

    # disable dummy subplot at [2, 2]
    dummy_ax.axis('off')

    # save plot
    plt.tight_layout()
    plt.savefig(file_name, dpi=200)
    plt.close()

def visualize_OP05withS_cross_val_tables(cross_val_result_Etot, cross_val_result_titers, cross_val_result_yields, scoring_result, file_name):
    """Plot the results of the cross-validation of OP05_withS optimization results (E_tot bar plot, titer heatmap, yield heatmap and score bar plot). The rows of the heatmaps correspond to the parameter sets of the ensemble and the columns correspond to the different optimization results (sets of initial enzyme and substrate concentrations). The value of each cell is then derived from simulating the model with a particular combination of a parameter set and an optimization result. The central diagonals correspond to cases where the parameter set was the one used to calculate the optimization in the first place whereas the off-diagonal elements are cases where an optimization result is combined with a new parameter set that was not used to solve it. The scores are calculated as described in the relevant scoring function. The best scoring optimization result is highlighted in all plots.

    :param cross_val_result_Etot: cross-validation result table (enzyme loads)
    :type cross_val_result_Etot: list
    :param cross_val_result_titers: cross-validation result table (titers)
    :type cross_val_result_titers: list
    :param cross_val_result_yields: cross-validation result table (yields)
    :type cross_val_result_yields: list
    :param scoring_result: scoring results of all optimization results
    :type scoring_result: pandas.core.frame.DataFrame
    :param file_name: name of the generated image file
    :type file_name: string
    """

    # get index of best scoring optimization result
    best_total_score_index_name = scoring_result.loc[:,'total_score'].idxmax()
    best_total_score_index = list(scoring_result.index).index(best_total_score_index_name)

    # setup multiplot
    fig, ((ax1, dummy_ax1), (ax2, cbar_ax2), (ax3, cbar_ax3), (ax4, dummy_ax4)) = plt.subplots(nrows=4, ncols=2, figsize=(7, 9), sharex='col',
                                                                                               gridspec_kw={'height_ratios': [1, 10, 10, 1], 'width_ratios': [20, 1]})

    # BAR PLOT - ENZYME LOAD
    # create the bar plot
    bars = ax1.bar([i + 0.5 for i in range(len(cross_val_result_Etot[0]))],
                   cross_val_result_Etot[0],
                   align='center')
    # disable x-axis tick marks and labels
    ax1.tick_params(axis='x',
                    which='both',
                    bottom=False,
                    top=False,
                    labelbottom=False,
                    labeltop=False)
    # enable y-axis tick marks
    ax1.tick_params(axis='y',
                    which='both',
                    left=True,
                    length=7,
                    direction='out')
    # set the axis labels
    ax1.set_xlabel('')
    ax1.set_ylabel('Enzyme \n Load [g/l]')
    # highlight the bar of the best scoring optimization result
    bars[best_total_score_index].set_color('red')
    # remove interior grid lines
    ax1.grid(False, axis='both')

    # disable dummy subplot at [1, 2]
    dummy_ax1.axis('off')

    # HEATMAP - TITERS
    hm = sns.heatmap(pd.DataFrame(cross_val_result_titers),
                     cmap=sns.cm.rocket_r,
                     xticklabels=False,
                     yticklabels=True,
                     cbar_ax=cbar_ax2,
                     ax=ax2)
    hm.tick_params(left=True, bottom=False)
    hm.set_ylabel('Parameter Sets')
    cbar_ax2.set_ylabel('\n Concentration [mM]')
    # set y-axis tick labels: only display labels for multiples of 5 and the first parameter set
    num_rows = len(cross_val_result_titers)
    yticks_pos = []
    yticks_labels = []
    for i in range(num_rows):
        if (i + 1) % 5 == 0 or (i + 1) == 1:
            yticks_pos.append(i + ((i+1-i)*0.5))  # place tick mark in the middle of the row
            yticks_labels.append(f'P{i + 1}')
    # set the y-ticks and labels
    ax2.set_yticks(yticks_pos)
    ax2.set_yticklabels(yticks_labels)
    # highlight best scoring column of the heatmap
    x_start = best_total_score_index
    x_end = best_total_score_index + 1
    y_start = 0
    y_end = len(cross_val_result_titers)
    rect = plt.Rectangle((x_start, y_start),
                         x_end - x_start,
                         y_end - y_start,
                         linewidth=2,
                         edgecolor='red',
                         facecolor='none',
                         transform=ax2.transData)
    ax2.add_patch(rect)

    # HEATMAP - YIELDS
    hm = sns.heatmap(pd.DataFrame(cross_val_result_yields),
                     cmap=sns.cm.mako_r,
                     xticklabels=False,
                     yticklabels=True,
                     cbar_ax=cbar_ax3,
                     ax=ax3)
    hm.tick_params(left=True, bottom=False)
    hm.set_ylabel('Parameter Sets')
    cbar_ax3.set_ylabel('\n Yield [-]')
    # set y-axis tick labels: only display labels for multiples of 5 and the first parameter set
    num_rows = len(cross_val_result_yields)
    yticks_pos = []
    yticks_labels = []
    for i in range(num_rows):
        if (i + 1) % 5 == 0 or (i + 1) == 1:
            yticks_pos.append(i + ((i+1-i)*0.5))  # place tick mark in the middle of the row
            yticks_labels.append(f'P{i + 1}')
    # set the y-ticks and labels
    ax3.set_yticks(yticks_pos)
    ax3.set_yticklabels(yticks_labels)
    # highlight best scoring column of the heatmap
    x_start = best_total_score_index
    x_end = best_total_score_index + 1
    y_start = 0
    y_end = len(cross_val_result_yields)
    rect = plt.Rectangle((x_start, y_start),
                         x_end - x_start,
                         y_end - y_start,
                         linewidth=2,
                         edgecolor='red',
                         facecolor='none',
                         transform=ax3.transData)
    ax3.add_patch(rect)

    # BAR PLOT - SCORES
    # create the bar plot
    bars2 = ax4.bar([i + 0.5 for i in range(len(scoring_result))],
                   scoring_result['total_score'],
                   align='center')
    # enable tick marks on both axis
    ax4.tick_params(axis='both', which='both', length=7, direction='out', bottom=True, left=True)
    # set the x-ticks and labels to only display x tick labels that are multiples of 5 (and the one of the very first opt. result)
    x_tick_labels = list()
    x_tick_pos = list()
    for i, opt_name in enumerate(scoring_result.index):
        opt_num = int(opt_name[1::])
        if opt_num % 5 == 0 or opt_num == 1:
            x_tick_labels.append(opt_name)
            x_tick_pos.append(i + 0.5)  # +0.5 to center the tick under the bar
    ax4.set_xticks(x_tick_pos)
    ax4.set_xticklabels(x_tick_labels)
    # set y-axis ticks and labels to show only 0 and 3 (lowest and highest possible scoring values)
    ax4.set_yticks([0, 3])
    ax4.set_yticklabels(['0', '3'])
    # set axis labels
    ax4.set_xlabel('Optimization Results')
    ax4.set_ylabel('Score [-] \n')
    # highlight the bar of the best scoring optimization result
    bars2[best_total_score_index].set_color('red')
    # remove interior grid lines
    ax4.grid(False, axis='both')

    # disable dummy subplot at [4, 2]
    dummy_ax4.axis('off')

    # save plot
    plt.tight_layout()
    plt.savefig(file_name, dpi=200)
    plt.close()

def create_boxplots_and_heatmaps(model_variant_name, opt_variant_name, path_to_files):
    """Visualize (1) the results of the repeated optimization for each parameter set of the ensemble (distributions of optimal enzyme and substrate concentrations) and (2) the results of the cross-validation (tables of simulated target values - titers, yields, enzyme loads - for all combinations of parameter sets and optimization results; enables selection of the best-scoring column, i.e., best-scoring optimization result).
    
    :param model_variant_name: name of the model variant that was used to calculate the optimizations
    :type model_variant_name: string
    :param opt_variant_name: optimization variant ('a' or 'b')
    :type opt_variant_name: string
    :param path_to_files: path to the result files of the optimization (this is also where the generated figures will be stored)
    :type path_to_files:
    """

    # set working directory to the provided path so that result files are found and generated figures are saved to the correct directory
    os.chdir(path_to_files)

    # 1) Repeated Optimization
    # ------------------------

    # show horizontal grid lines
    sns.set_theme(style='whitegrid')

    # OP01withS_7h
    visualize_opt_results(f'{model_variant_name}_PartSwarm50x_OP01withS_7h', opt_variant_name)
    plt.close()

    # OP01withS_12h
    visualize_opt_results(f'{model_variant_name}_PartSwarm50x_OP01withS_12h', opt_variant_name)
    plt.close()

    # OP05withS
    visualize_opt_results(f'{model_variant_name}_PartSwarm50x_OP05withS', opt_variant_name)
    plt.close()


    # 2) Cross-Validation
    # -------------------

    # OP01withS_7h
    # load results
    with open(f'{model_variant_name}_PartSwarm50x_OP01withS_7h_CrossVal.pkl', 'rb') as file:
        cross_val_result_OP01withS_7h = pickle.load(file)
    with open(f'{model_variant_name}_PartSwarm50x_OP01withS_7h_CrossVal_Scoring_Result_full.pkl', 'rb') as file:
        scoring_result_OP01withS_7h = pickle.load(file)
    # visualize results (titers heatmap and scores bar plot; highlight best scoring result)
    visualize_OP01_cross_val_tables(cross_val_result_OP01withS_7h, scoring_result_OP01withS_7h,
                                    f'{model_variant_name}_PartSwarm50x_OP01withS_7h_{opt_variant_name}_CrossVal_heatmap.png')
    plt.close()

    # OP01withS_12h
    # load results
    with open(f'{model_variant_name}_PartSwarm50x_OP01withS_12h_CrossVal.pkl', 'rb') as file:
        cross_val_result_OP01withS_12h = pickle.load(file)
    with open(f'{model_variant_name}_PartSwarm50x_OP01withS_12h_CrossVal_Scoring_Result_full.pkl', 'rb') as file:
        scoring_result_OP01withS_12h = pickle.load(file)
    # visualize results (titers heatmap and scores bar plot; highlight best scoring result)
    visualize_OP01_cross_val_tables(cross_val_result_OP01withS_12h, scoring_result_OP01withS_12h,
                                    f'{model_variant_name}_PartSwarm50x_OP01withS_12h_{opt_variant_name}_CrossVal_heatmap.png')
    plt.close()

    # OP05withS
    # load results
    with open(f'{model_variant_name}_PartSwarm50x_OP05withS_CrossVal_allEtot.pkl', 'rb') as file:
        cross_val_result_OP05withS_allEtot = pickle.load(file)
    with open(f'{model_variant_name}_PartSwarm50x_OP05withS_CrossVal_allTit.pkl', 'rb') as file:
        cross_val_result_OP05withS_allTit = pickle.load(file)
    with open(f'{model_variant_name}_PartSwarm50x_OP05withS_CrossVal_allYields.pkl', 'rb') as file:
        cross_val_result_OP05withS_allYields = pickle.load(file)
    with open(f'{model_variant_name}_PartSwarm50x_OP05withS_CrossVal_Scoring_Result_full.pkl', 'rb') as file:
        scoring_result_OP05withS = pickle.load(file)
    # visualize results (E_tot bar plot, titers heatmap, yields heatmap and scores bar plot; highlight best scoring result)
    visualize_OP05withS_cross_val_tables(cross_val_result_OP05withS_allEtot, cross_val_result_OP05withS_allTit,
                                         cross_val_result_OP05withS_allYields, scoring_result_OP05withS,
                                         f'{model_variant_name}_PartSwarm50x_OP05withS_{opt_variant_name}_CrossVal_heatmaps.png')
    plt.close()


# create boxplots and heatmaps of the results of both optimization variants (each for the selected model variant)

working_dir = os.getcwd()

# optimization variant a: PolyP fixed at 32 mM; ATP as optimization variable with its upper bound set to 2.5 mM
create_boxplots_and_heatmaps('ImpExtSearch_v23b_rep2_MV5', 'a', f'{working_dir}\\a_constPolyP32mM_ubATP2.5mM\\rep2_MV5')

# optimization variant b: PolyP fixed at 32 mM and ATP fixed at 0.5 mM
create_boxplots_and_heatmaps('ImpExtSearch_v23b_rep2_MV5', 'b', f'{working_dir}\\b_constPolyP32mMATP0.5mM\\rep2_MV5')
