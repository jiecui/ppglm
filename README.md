# Forecasting over days

This is the companion code for the article:
Proix T, Truccolo W, Leguia M, King-Stephens D, Rao V, Baud M (2021). Forecasting seizure risk in adults with focal epilepsy: a development and validation study. The Lancet Neurology, 20(1):127-135. https://doi.org/10.1016/S1474-4422(20)30396-3

### Installation:
Software:
- Matlab R2017b (for the processing).
- R 3.3.3 (for the forecasting).

Additionally, you need the following R libraries for the forecasting:
- tscount 1.4.1
- pROC 1.10.0
- R.utils 2.10.1
- feather 0.3.5
- pracma 2.3.3
- ncdf4 1.17
- reshape2 1.4.4

### Content
- plot: code to reproduce figures of the paper
- process: processing code
- run: code to obtain the results of the paper
- src: sources, main code
- stats: code to reproduce stats of the paper
- utils: utility functions

### Order in which to run the files 

#### Processing (process folder)
1. process_data
2. process_data_surrogates

#### Simulations (run folder)
1. test_no_cross_auto 
2. PSE 
3. select_optimal_history 
4. ahead_prediction 
5. minimum_training_duration 
6. online_retraining 
 
#### Surrogates (run folder)
1. test_no_cross_auto_surrogates 
2. concatenate_no_cross_auto_surrogates 
3. PSE_surrogates
4. select_optimal_history_surrogates 
5. concatenate_PSE_surrogates

#### Figures (plot folder)
1. plot_AUC
2. plot_BSS
3. plot_ahead_prediction
4. plot_minimal_training
5. plot_online_retraining

#### Stats (stats folder)
1. stats_AUC
2. stats_BSS

#### Citation

Please cite the following reference in your publications if you have used our code:
Proix T, Truccolo W, Leguia M, King-Stephens D, Rao V, Baud M (2021). Forecasting seizure risk in adults with focal epilepsy: a development and validation study. The Lancet Neurology, 20(1):127-135. https://doi.org/10.1016/S1474-4422(20)30396-3

#### Licensing

Copyright (c) 2017-2021, Timothée Proix, Maxime Baud. All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, version 3 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.