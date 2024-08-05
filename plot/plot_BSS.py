# Author: Timothée Proix
# License: GPL-3.0-only

import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import os
import feather
from scipy import stats 
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
    all_scores = feather.read_dataframe(os.path.join(path_res, 'all_patients_optimum_test_' + sampling + '_' + cov_type + '.feather'))
    all_ts = feather.read_dataframe(os.path.join(path_res, 'all_patients_ts_test_' + sampling + '_' + cov_type + '.feather'))
    all_scores_nca = feather.read_dataframe(os.path.join(path_res, 'all_patients_no_cross_auto_optimum_test_' + sampling + '_' + cov_type + '.feather'))
    all_ts_nca = feather.read_dataframe(os.path.join(path_res, 'all_patients_no_cross_auto_ts_test_' + sampling + '_' + cov_type + '.feather'))
except FileNotFoundError as e:
    print(e)
    print("file not found")

all_scores = pd.concat([all_scores, all_scores_nca])
all_ts = pd.concat([all_ts, all_ts_nca])
all_scores['dataset'] = dataset
all_ts['dataset'] = dataset

# removing NA
all_scores = all_scores.dropna()
all_scores = all_scores.reset_index(drop=True)
all_ts = all_ts.dropna()
all_ts = all_ts.reset_index(drop=True)

# merging patients with two blocks
for index, iRow in all_scores[all_scores['patient'].str.contains('2')].iterrows():
    jRow = all_scores[(all_scores['patient']==iRow['patient'][:2]) & (all_scores['cov_choice']==iRow['cov_choice']) & (all_scores['score_method']==iRow['score_method'])]
    if len(jRow)>0:
        all_scores.at[jRow.index.values[0], 'score'] = (jRow['score'] + iRow['score'])/2
        all_scores = all_scores.drop([iRow.name])

for ipatient in all_ts[all_ts['patient'].str.contains('2')]['patient'].unique():
    print(ipatient)
    for icov_choice in all_ts['cov_choice'].unique():
        for iscore_method in all_ts['score_method'].unique():
            iRows = all_ts[(all_ts['patient']==ipatient) & (all_ts['cov_choice']==icov_choice) & (all_ts['score_method']==iscore_method)]
            if len(iRows)>1:
                time_1 = all_ts[(all_ts['patient']==ipatient[:2]) & (all_ts['cov_choice']==icov_choice) & (all_ts['score_method']==iscore_method)]['time'].max()
                iRows.loc[:, 'patient']=iRows.loc[:, 'patient'].values[0][:2] 
                iRows.loc[:, 'time']=iRows.loc[:, 'time'] + time_1 
                all_ts = all_ts.append(iRows)
                all_ts = all_ts.drop(iRows.index)
            elif len(iRow)==1:
                print('****problem only one row*********')

# Brier skill score per patient
lBSS_all_rand, lpatients_all, lcov_choice_all = [], [], []
for icov_choice in lcov_choice:
    print(icov_choice)
    lpatients =  np.unique(all_ts[all_ts['dataset']==dataset]['patient'])
    for ipatient, patient in enumerate(lpatients):
        ts = all_ts[(all_ts['patient']==patient) & (all_ts['cov_choice']==icov_choice) & (all_ts['score_method']=='prop_time_corr_1') & (all_ts['dataset']==dataset)]
        y = ts['test_Sz'].values.astype('int')
        pk = ts['preds_test'].values
        pk[pk>1] = 1
        BS = np.mean((pk-y)**2)
        BSrand = []
        for j_ran in range(1000):
            BSrand.append(np.mean((np.random.choice(pk, pk.shape[0])-y)**2))
        BSrand = np.mean(BSrand)
        BSSrand = 1-(BS/BSrand)
        lBSS_all_rand.append(BSSrand)
        lpatients_all.append(patient)
        lcov_choice_all.append(icov_choice)

dfBSSrand = pd.DataFrame.from_items([('patient', lpatients_all), 
                                 ('BSS', lBSS_all_rand), 
                                 ('cov_choice', lcov_choice_all), 
                                 ('dataset', dataset)])

if sampling=='day':
    g = pp.RainCloud(x="cov_choice", y="BSS", data=dfBSSrand, color='g', order=lcov_choice, width_box=0.2, ax=ax, alpha=0.55, dodge=True)
    ax.set_xticklabels(['Recent\nseizures',  'Weekly\ndistribution', 'Recent\nIEA', 'Multidien\nphases', 'Multivariate'])
elif sampling=='hour':
    lcolor = [plt.get_cmap('Set2')(0)]*8
    g = pp.RainCloud(x="cov_choice", y="BSS", data=dfBSSrand, order=lcov_choice, color='g', width_viol=5, alpha=1, palette=lcolor)
    ax.set_xticklabels(['Recent\nseizures', 'Recent\nIEA', 'Circadian\nphases', 'Multidien\nphases', 'Multivariate'])
ax.set_xlabel('')
if sampling=='day':
    ax.set_ylim([-.15, .65])
    ax.set_xlim([-.6, 3.5])
elif sampling=='hour':
    ax.set_ylim([-.03, .12])
    ax.set_xlim([-1, 4.5])    
plt.ylabel('BSS')
sns.despine(offset=0, trim=True);
plt.show()