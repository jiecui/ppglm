# Author: Timothée Proix
# License: GPL-3.0-only

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import os
import feather
import numpy as np
import pdb
import scipy.io as sio
sns.set_context('paper')
sns.set(font_scale=1.1)
sns.set_style("ticks")

rootResults = 'YourPath'
dataset = 'UCSF_patients' 
path_res = os.path.join(rootResults, 'forecasting_over_days', dataset)

sampling = 'day' # day, hour
cov_type = 'acausal4'
cov_choice = 'phase_multi'


min_train_day = feather.read_dataframe(os.path.join(path_res, 'all_patients_minimal_training_day_' + cov_type + '.feather'))
min_train_day['dataset'] = dataset
min_train_day = min_train_day.loc[min_train_day['score_method']=='prop_time_corr_1']
min_train_day = min_train_day.loc[min_train_day['cov_choice']==cov_choice]
min_train_day = min_train_day.loc[min_train_day['training_period']<=250]

# merging patients with two blocks
for index, iRow in min_train_day[min_train_day['patient'].str.contains('2')].iterrows():
    jRow = min_train_day[(min_train_day['patient']==iRow['patient'][:2]) & (min_train_day['cov_choice']==iRow['cov_choice']) & (min_train_day['score_method']==iRow['score_method']) & (min_train_day['training_period']==iRow['training_period'])]
    if len(jRow)>0:
        min_train_day.at[jRow.index.values[0], 'score'] = (jRow['score'] + iRow['score'])/2
        min_train_day = min_train_day.drop([iRow.name])


ax = sns.lineplot(x='training_period', y='score', 
                  hue='dataset', data=min_train_day,  linewidth=1., 
                  sizes=10, dashes=False,
                  palette='Set2', legend=False)
ax.set_xlabel('Length of data history for training [days]')
ax.set_ylabel('AUC')
ticks = list(range(0, 251, 50))
ax.set_xticks(ticks)
ticklabs = list(range(0, 251, 50))
ax.set_xticklabels(ticklabs)
ax.set_xlim([0, 251])
ax.set_ylim([0.2, 1])
sns.despine(offset=0, trim=True);
plt.show()

all_min_train_hour = feather.read_dataframe(os.path.join(path_res, 'all_patients_minimal_training_hour_' + cov_type + '.feather'))
min_train_hour = all_min_train_hour.loc[all_min_train_hour['score_method']=='prop_time_corr_1']
min_train_hour = min_train_hour.loc[min_train_hour['training_period']<=250*24]
lpatients = np.unique(min_train_hour['patient'])

ltraining = []
for patient in lpatients:
    mtd = min_train_hour[min_train_hour['patient']==patient]['score'].values
    mtdmax = np.max(mtd)
    lmtd = np.nonzero(mtd>mtdmax*.95)[0]
    idx95 = lmtd[np.nonzero(np.diff(lmtd)==1)[0][0]]
    ltraining.append(min_train_hour[min_train_hour['patient']==patient].iloc[idx95]['training_period'])

print(np.mean(ltraining)/24.)

fig = plt.figure(figsize = (7, 5))
palette = ['k']*len(lpatients)
palette = sns.color_palette("Greys", len(lpatients))
lcolor = [plt.get_cmap('Set2')(0)]
ax = sns.lineplot(x='training_period', y='score', 
                  data=min_train_hour, linewidth=2., 
                  dashes=False,
                  palette=lcolor, legend=False, ci="sd")
ax.set_xlabel('Length of data history for training [days]')
ax.set_ylabel('AUC')
ax.set_xlim([0, 251*24])
ax.set_ylim([0.2, 1])
ticks = list(range(0, 251*24, 50*24))
ax.set_xticks(ticks)
ticklabs = list(range(0, 251, 50))
ax.set_xticklabels(ticklabs)
sns.despine(offset=0, trim=True);
plt.show()

