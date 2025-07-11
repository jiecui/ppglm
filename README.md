# Forecasting Seizures Over Days using PPGLM

This is the modified code for the article

Proix T, Truccolo W, Leguia M, King-Stephens D, Rao V, Baud M (2021). Forecasting seizure risk in adults with focal epilepsy: a development and validation study. The Lancet Neurology, 20(1):127-135. https://doi.org/10.1016/S1474-4422(20)30396-3.

## Code website

[Original code](https://zenodo.org/records/4274624)

## Installation

Software:

- MATLAB R2017b (for the processing)
- R 3.3.3 (for the forecasting).

Additionally, you need the following R libraries for the forecasting:

- tscount 1.4.1
- pROC 1.10.0
- R.utils 2.10.1
- feather 0.3.5
- pracma 2.3.3
- ncdf4 1.17
- reshape2 1.4.4

## Content of folders

- __plot__: code to reproduce figures of the paper
- __process__: processing code
- __run__: code to obtain the results of the paper
- __src__: sources, main code
- __stats__: code to reproduce stats of the paper
- __utils__: utility functions

## Order in which to run the files

### Processing (__process__ folder)

1. process_data
2. process_data_surrogates

### Simulations (__run__ folder)

1. test_no_cross_auto
2. PSE
3. select_optimal_history
4. ahead_prediction
5. minimum_training_duration
6. online_retraining

### Surrogates (__run__ folder)

1. test_no_cross_auto_surrogates 
2. concatenate_no_cross_auto_surrogates 
3. PSE_surrogates
4. select_optimal_history_surrogates 
5. concatenate_PSE_surrogates

### Figures (__plot__ folder)

1. plot_AUC
2. plot_BSS
3. plot_ahead_prediction
4. plot_minimal_training
5. plot_online_retraining

### Stats (__stats__ folder)

1. stats_AUC
2. stats_BSS
