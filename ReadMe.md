# Combining Automated Kinetic Model Selection with Robust Optimization and Application to Maximize Cell-free Production of UDP-GalNAc

## Introduction
This repository contains all code and experimental data used to generate the results presented in [Huber et al., 2026](#publication). The implementation and application of *MoRSel*, the **Mo**del **R**efinement and **Sel**ection method that we introduce and apply in this work, can be found in the [eponymous directory](/morsel). All files of the robust ensemble-based process optimization can be found in the [process optimization directory](/process_optimization).

<details>

<summary>Python Requirements</summary>

All code was developed using Python version 3.12.2 as well as the following packages:

| Python Package |  Version |
| -------------- | -------- |
| copasi_basico  | 0.58     |
| joblib         | 1.4.0    |
| matplotlib     | 3.10.7   |
| numpy          | 2.3.4    |
| pandas         | 2.3.3    |
| python-copasi  | 4.42.284 |
| seaborn        | 0.13.2   |
| tqdm           | 4.67.1   |

</details>

## *MoRSel* User Guide

### 1. Implementing a model
The [*MoRSel* function library](/morsel/func_lib.py) contains the definitions of all necessary functions grouped by type. In the first section, the model-specific functions are defined. Two different models are already implemented: the UDP-GalNAc model and the toymodel. Each is represented by two functions: *\<model-name\>_data* and *\<model-name\>_term_libs*. In order to add a new model, both of these functions need to be defined.

#### *\<model-name\>_data* function
Defines and returns dictionaries with the following information: (1) stoichiometry of all model reactions (using [Copasi syntax](https://copasi.org/Support/User_Manual/Model_Creation/Reactions/)), (2) any known values of kinetic model parameters, and (3) initial substrate concentrations. The latter are set as placeholder values to fallback on if the model is simulated without specifying any initial conditions.

<details>

<summary>Python Implementation</summary>

```python
def toymodel_data():

    # all reactions are defined as reversible ('=')
    reaction_scheme_dict = {'r1': 'S = P',
                            'r2': 'S = X',
                            'r3': 'P = X',
                            'r4': 'S + P = X'}

    # here, no prior knowlegde of kinetic parameter values is available
    param_dict = {'r1': {},
                  'r2': {},
                  'r3': {},
                  'r4': {}}

    init_conc_dict = {'S': 10}  # initial concentration of the substrate S [mmol/l]

    return reaction_scheme_dict, param_dict, init_conc_dict
```

</details>

#### *\<model-name\>_term_libs* function
Defines and returns (1) a dictionary that for each reaction contains all possible kinetic rate law variants and (2) a list of all possible multiplicative terms that can be included in the rate laws.

<details>

<summary>Python Implementation</summary>

```python
def toymodel_term_libs():

    # two different rate laws were defined for each reaction: a zero rate law at index 0 (so that the reaction can be turned off)
    # and a mass action rate law at index 1
    base_kinetics = {
        'r1': [{'equation': '0', 
                'mapping': {}},
               {'equation': '(k_MAforward*S-k_MAreverse*P)', 
                'mapping': {'k_MAforward': 'parameter', 'k_MAreverse': 'parameter', 'S': 'substrate', 'P': 'product'}}],
        'r2': [{'equation': '0', 
                'mapping': {}},
               {'equation': '(k_MAforward*S-k_MAreverse*X)', 
                'mapping': {'k_MAforward': 'parameter', 'k_MAreverse': 'parameter', 'S': 'substrate', 'X': 'product'}}],
        'r3': [{'equation': '0', 
                'mapping': {}},
               {'equation': '(k_MAforward*P-k_MAreverse*X)', 
                'mapping': {'k_MAforward': 'parameter', 'k_MAreverse': 'parameter', 'P': 'substrate', 'X': 'product'}}],
        'r4': [{'equation': '0', 
                'mapping': {}},
               {'equation': '(k_MAforward*S*P-k_MAreverse*X)', 
                'mapping': {'k_MAforward': 'parameter', 'k_MAreverse': 'parameter', 'S': 'substrate', 'P': 'substrate', 'X': 'product'}}]
    }

    # the list of multiplicative terms covers all possible allosteric inhibitions and activations
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
```

</details>

### 2. Adding experimental data
Only time-resolved data can be used with entries being separated by tabs. Column headers are constructed from species names enclosed by square brackets. Initial conditions can be added as extra columns with a "\_0" added to the header. Multiple files following this basic structure can be added.

<details>

<summary>Example Data</summary>

| Time | \[P\] | \[S\] | \[X\] | \[P\]\_0 | \[S\]\_0 | \[X\]\_0 |
| ---- | ----- | ----  | ----- | -------- | -------- | -------- |
| 0.0  | 0.0   | 10.0  | 0.0   | 0.0      | 10.0     | 0.0      |
| 2.0  | 0.9   |  3.2  | 2.8   |          |          |          |
| 4.0  | 1.4   |  2.3  | 3.1   |          |          |          |

</details>

Check the [data](/morsel/toymodel/true_model_data.txt) generated by the toymodel used for the rediscovery test as another example.

### 3. Parameterizing a model variant
The parameters of a model variant can be estimated using the added experimental data following a least-squares approach. The parameterization can be repeated to handle parametric uncertainty resulting in an ensemble of parameter sets.

First, the structure of the model variant needs to be defined. This is done by combining certain selected kinetic rate laws and multiplicative terms. The structure of the model variant is represented by a so called "kinetic table". Its elements are the indices of the selected kinetic rate law variants ("base") and multiplicative terms ("terms") as they are stored in the objects defined by the *\<model-name\>_term_libs* function.

<details>

<summary>Structure of the Kinetic Table</summary>

|            | base                           | term 1 | ... | term j |
| :--------: | :----------------------------: | :----: | :--:| :----: |
| reaction 1 | selected base index reaction 1 | selected multiplicative term 1 of reaction 1 | ... | selected multiplicative term j of reaction 1 |
| ...        | ...                            | ...                                          | ... | ...                                          |
| reaction i | selected base index reaction i | selected multiplicative term 1 of reaction i | ... | selected multiplicative term j of reaction i |

</details>

Then, the experimental data is loaded and the least-squares parameter estimation is computed following the logic of the *create_parameter_ensemble* function defined in the [function library](/morsel/func_lib.py). Alternatively, the start and boundary values of some or all parameters to be estimated can be set manually in a dictionary that is then passed to the *create_parameter_ensemble* function.

<details>

<summary>Python Implementation</summary>

```python
from func_lib import *

# set a name for the model variant
model_name = 'model_variant'

# define the structure of the model variant
# here, only reaction 1 is active and no multiplicative terms are included in the rate laws
model_var = [1,    # r1: S = P
             0,    # r2: S = X
             0,    # r3: P = X
             0]    # r4: S + P = X
model_var = pd.DataFrame(np.array(model_var).reshape(4, 1),
                         columns=['base'],
                         index=['r1', 'r2', 'r3', 'r4'])

# load experimental data
exp_data_file_names = ['path/to/data.txt']
exp_data_dataframes = [pd.read_csv(file_name, sep='\t') for file_name in exp_data_file_names]

# example of an optional dictionary providing a specific parameter estimation setup which contains 
# start values, upper and lower boundaries of selected parameters to be estimated - the start and 
# boundary values of all remaining parameters that are not defined here are set according to the 
# logic implemented in the eval_struct_var function called by create_parameter_ensemble (this is 
# the case for all parameters if this dictionary is not supplied); the syntax of the parameter keys 
# is '(reaction_name).parameter_name'
PE_setup = {'(r1).k_MAforward': {'start': 0.01, 'lower' 1e-2:, 'upper': 10},
            '(r3).k_MAforward': {'start': 0.5, 'lower': 0.13, 'upper': 1e1},
            '(r3).k_MAreverse': {'start': 1.9, 'lower': 1.75, 'upper': 2.1}}

# repeat the evaluation of the same model variant 10 times and store the results; the start and 
# boundary values of three selected parameters are defined according to the dictionary that is 
# provided, the rest are defined according to the internal logic of eval_struct_var called by 
# the function create_parameter_ensemble
n_runs = 10
evaluated_vars_log = create_parameter_ensemble(model_var, n_runs, exp_data_file_names, exp_data_dataframes, model_name,
                                               fit_params_info=PE_setup, model_data=toymodel_data, term_libs=toymodel_term_libs)
```

</details>

The calculated parameter ensemble is stored as a pickled object and as a csv file (with each parameter set in a new row and separate columns for each parameter).

### 4. Running *MoRSel*'s model refinement and selection loop
For a given initial core model variant and a library of candidate model changes, MoRSel can explore the resutling solution space to find the modified model variant that best describes the experimental data. This search method is implemented in the *improved_extension_search* function defined in the [*MoRSel function library*](/morsel/func_lib.py). The name "improved extension search" was used during development to refer to *MoRSel*'s refinement and selection loop.

<details>

<summary>Python Implementation</summary>

```python
# define initial core model (start model variant)
#            base term1
start_var = [1,   0,    # r1: S = P
             0,   0,    # r2: S = X
             0,   0,    # r3: P = X
             0,   0]    # r4: S + P = X
start_var = pd.DataFrame(np.array(start_var).reshape(4, 2),
                         columns=['base', 'term1'],
                         index=['r1', 'r2', 'r3', 'r4'])

# define list of candidate model changes (the first two elements of each tuple are the row and 
# column coordinates pointing to elements in the start variant data frame and the third element 
# is the index of the selected kinetic rate law or multiplicative term) - also: tuples are placed
# inside inner lists because combinations of multiple tuples are also valid candidate model changes
vari_terms = [[(1,0,1)],  # r2.on
              [(2,0,1)],  # r3.on
              [(3,0,1)],  # r4.on
              [(0,1,5)],  # r1.X_Act
              [(0,1,6)]]  # r1.X_Inhib

# load experimental data
exp_data_file_names = ['true_model_data.txt']
exp_data_dataframes = [pd.read_csv(file_name, sep='\t') for file_name in exp_data_file_names]

# define which model selection criteria are to be used
selected_MSC = ['AIC', 'AICc', 'BIC', 'CIC1', 'CIC2', 'CIC3']

# use the existing parameter ensemble of the initial core model (here also referred to as the negative control)
# alternatively, the initial core model can be evaluated inside the improved_extension_search function
reading_file = open('parameter_ensemble_pickle_file_name', 'rb')
neg_ctrl_log = pickle.load(reading_file)
reading_file.close()
start_var_ensemble = neg_ctrl_log

# run the MoRSel model refinement and selection loop
# (here referred to by its old name as "improved extension search")
output_dict = improved_extension_search(start_var, vari_terms,
                                        exp_data_file_names, exp_data_dataframes,
                                        PE_algorithm_info=PE_algorithm_info,
                                        selected_MSC=selected_MSC,
                                        start_var_ensemble=start_var_ensemble,
                                        model_data=toymodel_data,
                                        term_libs=toymodel_term_libs)

# store the result 
model_selection_result_name = 'result_name'
storage_file = open(f'{model_selection_result_name}_output_dict', 'wb')
pickle.dump(output_dict, storage_file)
storage_file.close()
```

</details>

### 5. Selecting the overall best model variant
The result of MoRSel's refinement and selection loop is a list of evaluated model variants (the best of each iteration). A final selection step can be performed to select the overall best model variant from all elements of that list.

<details>

<summary>Python Implementation</summary>

```python
# load a result of a MoRSel refinement and selection loop
reading_file = open(f'{model_selection_result_name}_output_dict', 'rb')
output_dict = pickle.load(reading_file)
reading_file.close()

# analyze result and identify best model variant
analysis_result = output_dict['analysis_result']
# get the index of the model variant with the best median rank and the fewest additions
best_median_rank_model_idx = analysis_result["fitness_dataframe"]['median_rank'].idxmin()
# get the index of the model variant with the best mean rank and the fewest additions
best_mean_rank_model_idx = analysis_result["fitness_dataframe"]['mean_rank'].idxmin()
# compare both and select the model with the lowest index (= lowest number of additions)
selected_idx = 0
if best_median_rank_model_idx < best_mean_rank_model_idx:
    best_ranking_var_log_dict = output_dict['evaluated_vars_log'][best_median_rank_model_idx]
    print(f"Model variant {best_median_rank_model_idx} with best median rank {analysis_result['fitness_dataframe'].iloc[best_median_rank_model_idx, :]['median_rank']} was selected.")
    selected_idx = best_median_rank_model_idx
elif best_median_rank_model_idx > best_mean_rank_model_idx:
    best_ranking_var_log_dict = output_dict['evaluated_vars_log'][best_mean_rank_model_idx]
    print(f"Model variant {best_mean_rank_model_idx} with best mean rank {analysis_result['fitness_dataframe'].iloc[best_mean_rank_model_idx, :]['mean_rank']} was selected.")
    selected_idx = best_mean_rank_model_idx
elif best_median_rank_model_idx == best_mean_rank_model_idx:
    # both indices are the same so they both point to the same model variant 
    # -> therefore it doesn't matter which one is used to look it up in the ImpExtSearch result
    best_ranking_var_log_dict = output_dict['evaluated_vars_log'][best_median_rank_model_idx]
    print(f"Model variant {best_median_rank_model_idx} with best median rank {analysis_result['fitness_dataframe'].iloc[best_median_rank_model_idx, :]['median_rank']} and best mean rank {analysis_result['fitness_dataframe'].iloc[best_mean_rank_model_idx, :]['mean_rank']} was selected.")
    selected_idx = best_median_rank_model_idx

# create parameter ensemble of the selected model variant 
model_variant_name = 'selected_model_variant_name'
evaluated_vars_log = create_parameter_ensemble(best_ranking_var_log_dict['variant'], n_runs, 
                                               exp_data_file_names, exp_data_dataframes, model_variant_name,
                                               fit_params_info=fit_params_info, model_data=toymodel_data, term_libs=toymodel_term_libs)

# get best parameter set of the ensemble
final_ensemble = pd.read_csv('selected_model_variant_ensemble.csv', sep=',')
final_ensemble_best_params = final_ensemble.loc[final_ensemble['obj'] == final_ensemble['obj'].min()]

# save result
result_log = {'best_ranking_model_variant_idx': selected_idx,
              'best_ranking_model_variant_structure': best_ranking_var_log_dict['variant'],
              'best_param_set': final_ensemble_best_params}
with open('rediscovery_test_result.txt', 'w') as file:
    for key,value in result_log.items():
        file.write("\n%s:\n%s\n" % (key,value))

```

</details>

## Applications of *MoRSel* in this Work

### Rediscovery Test (Toymodel)
The *MoRSel* method was evaluated by performing a so called "rediscovery test". Its aim is to test whether, given data generated from a hidden true model, *MoRSel* is able to rediscover the structure (and parameterization) of this true model. The rediscovery test was performed using the toymodel as defined in the [*MoRSel* function library](/morsel/func_lib.py) and was implemented in [this Python script file](/morsel/toymodel/rediscovery_test.py) with all results stored in the same directory.

### UDP-GalNAc Model
*MoRSel* was applied to refine an initial model of the UDP-GalNAc process in order to generate a modified model variant that better describes the data. The complete framework was repeated five times. The results of the first run are located directly in the [MoRSel directory](/morsel). The results of the other four repetitions are located in the directories: [rep1](/morsel/rep1), [rep2](/morsel/rep2), [rep3](/morsel/rep3), and [rep4](/morsel/rep4). A [report](/morsel/IESv23b_SSO_report.csv) that summarizes the results of all five runs was generated by running [this script](/morsel/create_SSO_report.py).

All experimental data sets of the UDP-GalNAc process were stored in text files labeled "UDP-GalNAc\<experiment-name\>_with_initConcColumns" ([example](/morsel/UDP-GalNAc36_50_with_initConcColumns.txt)). For further information on the experimental data used for *MoRSel* check the text and supplementary material of the publication.

Parameter ensembles were stored as csv files with the prefix "sampling_output_Particle_Swarm50runs" ([example](/morsel/sampling_output_Particle_Swarm50runs_v23bNegControl_Exp36.37.39.csv)). Here, "sampling" refers to repeating the same parameter estimation problem with a global solver (Particle Swarm [as implemented in Copasi](https://copasi.org/Support/User_Manual/Methods/Optimization_Methods/Particle_Swarm/)), which has some stochastic elements leading to the estimation of different parameter sets each run. The parameter ensembles of the [initial core model](/morsel/v23bNegCtrl_PartSwarm50runs_hist_fig.png) and the [selected best-ranking model](/morsel/v23brep2MV5_PartSwarm50runs_hist_fig.png) were visualized as histograms by [this script](/morsel/plot_param_ensembles.py).

For a full run of *MoRSel* two Python script files were executed:
1. "calc_\<run-index\>_ImpExtSearch_Particle_Swarm.py" ([example](/morsel/calc_v23b_SSO_ImpExtSearch_Particle_Swarm.py)) where the parameter ensemble of the initial core model and the refinement and selection loop are calculated. The numerical solver Particle Swarm was used for all parameter estimations (both for the initial ensemble and during the loop). The result was saved in a file with the suffix "_output_dict" ([example](/morsel/SSO_ImpExtSearch_v23b_output_dict)). This is a pickled Python dictionary which contains information on the best model variants of each iteration of the refinement and selection loop.
2. "calc_\<run-index\>_ImpExtSearchBestMedianResult_MV\<final-model-variant-index\>_Ensemble_PS50" ([example](/morsel/calc_v23b_SSO_ImpExtSearchBestMedianResult_MV5_Ensemble_PS50.py)) where the final selection (primarily based on the median rank across all model selection criteria) and the parameter ensemble of the selected model variant were calculated. 50 parameter sets were estimated using the Particle Swarm solver, hence the label "PS50".

The selected best-ranking model variant across all five runs is called "rep2_MV5" refering to the second repetition of *MoRSel* (third run in total) where model variant number 5 was chosen in the final selection step. Figures comparing the fit quality of the initial core model (dashed lines) to the fit quality of the selected model variant (solid lines) were created with [this script](/morsel/plot_NegCtrl_and_rep2MV5_time_courses.py) for all experimental data sets ([example](/morsel/v23bNegCtrldashed_v23IESrep2MV5solid_DataUDP-GalNAc36_50_Comparison_merged.png)).

A [Copasi model file of the selected best-ranking model variant]((/morsel/rep2/ImpExtSearch_v23b_rep2_MV5.cps) for subsequent process optimization was created by [this script](/morsel/rep2/create_Copasi_model_file.py).

## Process Optimization of the UDP-GalNAc Process
The selected model variant was used to optimize the cell-free production of UDP-GalNAc compared to the baseline behavior ([experiment 36, 50mM](/process_optimization/exp_val_data/UDP-GalNAc36_50_with_initConcColumns.txt)). Two different optimization scenarios were calculated:

1. Titer optimization: Maximize the UDP-GalNAc titer (at t = 7 hours) while limiting the required enzyme load (sum of all initial enzyme concentrations) to that of the baseline. ID: "OP01withS_7h".
2. Enzyme load optimization: Minimize the enzyme load that is required to reach the UDP-GalNAc titer and yield which was reach by the baseline process (at t = 24 hours). ID: "OP05withS".

There were some slight differences between the setups of the different optimization scenarios. For the titer optimization, polyphosphate (PolyP) is set to a constant value while adenosine triphosphate (ATP) is included as optimization variable with an upper bound of 2.5 mM ([setup A](/process_optimization/a_constPolyP32mM_ubATP2.5mM)). For the enzyme load optimization, both PolyP and ATP are set to constant values ([setup B](/process_optimization/b_constPolyP32mMATP0.5mM)).

In both cases the optimizations were performed according to the robust ensemble-based optimization methodology introduced by us in [Huber et al. 2024](https://doi.org/10.1016/j.ymben.2023.10.007).

The results of both optimizations were analyzed and visualized by two Python scripts: [*create_boxplots_and_heatmaps*](/process_optimization/create_boxplots_and_heatmaps.py) and [*visualize_selected_predictions*](/process_optimization/visualize_selected_predictions.py). The generated figures were saved in the model-specific subdirectories of the optimization setups A and B.

The data of the validation experiments is available in [this directory](/process_optimization/exp_val_data). In the file names, "L" refers to the enzyme load optimization and "T" refers to the titer optimization. "O" indicates that the selected model variant (with the **o**ptimized model structure) was used.

## Publication
Huber, Son, Espinel-Ríos, Rexer, Reichl, Klamt (2026), Combining Kinetic Model Selection and Robust Optimization Strategies to Maximize Production of UDP-GalNAc in a Cell-free Batch Process, in preparation