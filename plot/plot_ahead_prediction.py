# Author: Timothée Proix
# License: GPL-3.0-only

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import os
import feather
import numpy as np
import scipy.io as sio
import scipy.stats as ss
import ptitprince as pp
sns.set_context('paper')
sns.set(font_scale=1.1)
sns.set_style("ticks")

rootResults = 'YourPath'
dataset = 'UCSF_patients' 
path_res = os.path.join(rootResults, 'forecasting_over_days', dataset)

cov_type = 'acausal4'
cov_choice = 'phase_multi'
sampling = 'day'

all_ahead_day = feather.read_dataframe(os.path.join(path_res, 'all_patients_ahead_perf_day_' + cov_type + '.feather'))
all_ahead_day['dataset'] = dataset
ahead_day = all_ahead_day.loc[all_ahead_day['score_method']=='prop_time_corr_1']
ahead_day = ahead_day.loc[ahead_day['cov_choice']==cov_choice]
ahead_day = ahead_day[ahead_day['n_ahead'].isin(np.arange(15))]

# merging patients with two blocks
for index, iRow in ahead_day[ahead_day['patient'].str.contains('2')].iterrows():
    jRow = ahead_day[(ahead_day['patient']==iRow['patient'][:2]) & (ahead_day['cov_choice']==iRow['cov_choice']) & (ahead_day['score_method']==iRow['score_method']) & (ahead_day['n_ahead']==iRow['n_ahead']) & (ahead_day['dataset']==dataset)]
    if len(jRow)>0:
        ahead_day.at[jRow.index.values[0], 'score'] = (jRow['score'] + iRow['score'])/2
        ahead_day = ahead_day.drop([iRow.name])

fig = plt.figure()
ax = fig.add_subplot(111)
pp.RainCloud(x="n_ahead", y="score", data=ahead_day, color='g', width_box=0.2, ax=ax, alpha=0.55, dodge=True)
ax.plot([-0.6, 6], [0.5, 0.5], '--k', zorder=1)
ax.set_xlabel('Forecasting horizon [number of days in advance]')
ax.set_ylabel('AUC')
ax.set_xticklabels(np.arange(2, 7))
ax.set_xlim([.4, 6.3])
ax.set_ylim([0.3, 1.01])
sns.despine(offset=0, trim=True);
plt.show()

## Hours
cov_type = 'acausal4'
cov_choice = 'multivariate'
sampling = 'hour'

all_ahead_hour = feather.read_dataframe(os.path.join(path_res, 'all_patients_ahead_perf_hour_' + cov_type + '.feather'))
all_ahead_hour[dataset] = dataset
ahead_hour = all_ahead_hour.loc[all_ahead_hour['score_method']=='prop_time_corr_1']
ahead_hour = ahead_hour.loc[ahead_hour['cov_choice']==cov_choice]
ahead_hour = ahead_hour[ahead_hour['n_ahead'].isin(np.arange(1, 8))]

# merging patients with two blocks
for index, iRow in ahead_hour[ahead_hour['patient'].str.contains('2')].iterrows():
    jRow = ahead_hour[(ahead_hour['patient']==iRow['patient'][:2]) & (ahead_hour['cov_choice']==iRow['cov_choice']) & (ahead_hour['score_method']==iRow['score_method']) & (ahead_hour['n_ahead']==iRow['n_ahead'])]
    if len(jRow)>0:
        ahead_hour.at[jRow.index.values[0], 'score'] = (jRow['score'] + iRow['score'])/2
        ahead_hour = ahead_hour.drop([iRow.name])

lcolor = [plt.get_cmap('Set2')(0)]*8
fig = plt.figure()
ax = fig.add_subplot(111)
pp.RainCloud(x="n_ahead", y="score_value", data=ahead_hour, width_box=0.2, ax=ax, alpha=1, palette=lcolor)
ax.plot([-0.6, 7], [0.5, 0.5], '--k', zorder=1)
ax.set_xlabel('Forecasting horizon [number of hours in advance]')
ax.set_ylabel('AUC')
ax.set_xticklabels(np.arange(2, 16, 2))
ax.set_xlim([-.6, 8.3])
ax.set_ylim([0.3, 1.01])
sns.despine(offset=0, trim=True);
plt.show()