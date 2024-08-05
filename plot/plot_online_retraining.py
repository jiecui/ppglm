# plot figure 1 of the paper

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import os
import feather
import numpy as np
import pdb
import scipy.io as sio
import ptitprince as pp

sns.set_context('paper')
sns.set(font_scale=1.1)
sns.set_style("ticks")

rootResults = 'YourPath'
dataset = 'UCSF_patients' 
path_res = os.path.join(rootResults, 'forecasting_over_days', dataset)

## Days
cov_type = 'acausal4'
cov_choice = 'phase_multi'
sampling = "day"

all_retrain_day = feather.read_dataframe(os.path.join(path_res, 'all_patients_online_retraining_day_' + cov_type + '.feather'))
all_retrain_day['dataset'] = 'dataset'
retrain_day = all_retrain_day[all_retrain_day['retraining_period'].isin([-1, 1, 10, 50])]
retrain_day = retrain_day.loc[retrain_day['score_method']=='prop_time_corr_1']
retrain_day = retrain_day.loc[retrain_day['cov_choice']==cov_choice]

# merging patients with two blocks
for index, iRow in retrain_day[retrain_day['patient'].str.contains('2')].iterrows():
    jRow = retrain_day[(retrain_day['patient']==iRow['patient'][:2]) & (retrain_day['cov_choice']==iRow['cov_choice']) & (retrain_day['score_method']==iRow['score_method']) & (retrain_day['retraining_period']==iRow['retraining_period'])  & (retrain_day['dataset']==dataset)]
    if len(jRow)>0:
        retrain_day.at[jRow.index.values[0], 'score'] = (jRow['score'] + iRow['score'])/2
        retrain_day = retrain_day.drop([iRow.name])


fig = plt.figure()
ax = fig.add_subplot(111)
pp.RainCloud(x="retraining_period", y="score", data=retrain_day, color='g', width_box=0.2, ax=ax, alpha=0.55, dodge=True)
ax.set_xlabel('Frequency of retraining')
ax.set_ylabel('AUC score')
ax.set_ylim([0.4, 1.01])
ax.set_xlim([-0.6, 3.5])
sns.despine(offset=5, trim=True);
ax.set_xticklabels(['every seizure', 'every day', 'every 10 days', 'every 50 days'])
plt.show()

## Hours
cov_type = 'acausal4'
cov_choice = 'multivariate'
sampling = 'hours'

all_retrain_hour = feather.read_dataframe(os.path.join(path_res, 'all_patients_online_retraining_hour_' + cov_type + '.feather'))
all_retrain_hour['dataset'] = dataset
retrain_hour = all_retrain_hour[all_retrain_hour['retraining_period'].isin([-1, 1, 12, 24, 24*10])]
retrain_hour = retrain_hour.loc[retrain_hour['score_method']=='prop_time_corr_1']
retrain_hour = retrain_hour.loc[retrain_hour['cov_choice']==cov_choice]

# merging patients with two blocks
for index, iRow in retrain_hour[retrain_hour['patient'].str.contains('2')].iterrows():
    jRow = retrain_hour[(retrain_hour['patient']==iRow['patient'][:2]) & (retrain_hour['cov_choice']==iRow['cov_choice']) & (retrain_hour['score_method']==iRow['score_method']) & (retrain_hour['retraining_period']==iRow['retraining_period'])]
    if len(jRow)>0:
        retrain_hour.at[jRow.index.values[0], 'score'] = (jRow['score'] + iRow['score'])/2
        retrain_hour = retrain_hour.drop([iRow.name])

lcolor = [plt.get_cmap('Set2')(0)]*8
fig = plt.figure()
ax = fig.add_subplot(111)
pp.RainCloud(x='retraining_period', y='score_value', data=retrain_hour, ax=ax, palette=lcolor)
ax.set_xlabel('Frequency of retraining')
ax.set_ylabel('AUC score')
ax.set_ylim([0.4, 1.01])
sns.despine(offset=0, trim=True);
ax.set_xticklabels(['every seizure', 'every hour', 'every\n12 hours', 'every\n24 hours', 'every\n10 days'])
plt.show()