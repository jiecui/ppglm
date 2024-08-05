# Author: Timothée Proix
# License: GPL-3.0-only

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import os
import feather
import numpy as np
import ptitprince as pp

sns.set_context('paper')
sns.set(font_scale=1.1)
sns.set_style("ticks")

rootResults = 'YourPath'
dataset = 'UCSF_patients' 
path_res = os.path.join(rootResults, 'forecasting_over_days', dataset)

sampling = 'day' # day, hour
cov_type = 'acausal4'

if sampling=='day':
    lcov_choice = ['only_auto', 'only_cross', 'week_hist', 'phase_multi']
    lcov_choice_nca = ['week_hist', 'phase_multi']
elif sampling=='hour':
    lcov_choice = ['only_auto', 'circ_hist', 'only_cross', 'phase_circ', 'phase_multi', 'multivariate']
    lcov_choice_nca = ['circ_hist', 'phase_circ', 'phase_multi', 'multivariate']

try:
    all_scores_ca = feather.read_dataframe(os.path.join(path_res, 'all_patients_optimum_test_' + sampling + '_' + cov_type + '.feather'))
    all_scores_nca = feather.read_dataframe(os.path.join(path_res, 'all_patients_no_cross_auto_optimum_test_' + sampling + '_' + cov_type + '.feather'))
except FileNotFoundError as e:
    print(e)
    print("file not found")

score_selection_ca = all_scores_ca.loc[all_scores_ca['cov_choice'].isin(['only_auto', 'only_cross'])].copy()
if sampling=='day':
    score_selection_nca = all_scores_nca.loc[all_scores_nca['cov_choice'].isin(lcov_choice_nca)].copy()
elif sampling=='hour':
    score_selection_nca = all_scores_nca.loc[all_scores_nca['cov_choice'].isin(lcov_choice_nca)].copy()

score_selection = pd.concat([score_selection_ca, score_selection_nca])
score_selection['dataset'] = dataset
score_selection = score_selection.loc[score_selection['score_method']=='prop_time_corr_1']

# removing NA and reset index to guarantee uniqueness
score_selection = score_selection.dropna()
score_selection = score_selection.reset_index(drop=True)

# merging patients with two blocks
for index, iRow in score_selection[score_selection['patient'].str.contains('2')].iterrows():
    jRow = score_selection[(score_selection['patient']==iRow['patient'][:2]) & (score_selection['cov_choice']==iRow['cov_choice']) & (score_selection['score_method']==iRow['score_method'])]
    if len(jRow)>0:
        score_selection.at[jRow.index.values[0], 'score'] = (jRow['score'] + iRow['score'])/2
        score_selection = score_selection.drop([iRow.name])

# raincloud part
fig, ax = plt.subplots(figsize=(8,5))
if sampling=='day':
    g = pp.RainCloud(x="cov_choice", y="score", data=score_selection, color='g', order=lcov_choice, width_box=0.2, ax=ax, alpha=0.55, dodge=True)
elif sampling=='hour':
    lcolor = [plt.get_cmap('Set2')(0)]*8
    g = pp.RainCloud(x="cov_choice", y="score", data=score_selection, order=lcov_choice, palette=lcolor)
ax.set_xlabel('')
ax.set_ylim([0.3, 1])
plt.ylabel('AUC')
sns.despine(offset=0, trim=True);
plt.show()

