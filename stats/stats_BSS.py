# Author: Timothée Proix
# License: GPL-3.0-only

import numpy as np
import seaborn as sns
import pandas as pd
import os
import feather
from scipy import stats 
import pdb

rootResults = 'YourPath'
dataset = 'UCSF_patients' 
path_res = os.path.join(rootResults, 'forecasting_over_days', dataset)

sampling = 'day' # day, hour
cov_type = 'acausal4'
surr = 0

if sampling=='day':
    lcov_choice = ['only_auto', 'week_hist', 'only_cross', 'phase_multi']
    lcov_choice_nca = ['week_hist', 'phase_multi']
elif sampling=='hour':
    lcov_choice = ['only_auto', 'circ_hist', 'only_cross', 'phase_circ', 'phase_multi', 'multivariate']
    lcov_choice_nca = ['circ_hist', 'phase_circ', 'phase_multi', 'multivariate']

try:
    all_scores_ca = feather.read_dataframe(os.path.join(path_res, 'all_patients_optimum_test_' + sampling + '_' + cov_type + '.feather'))
    all_ts_ca = feather.read_dataframe(os.path.join(path_res, 'all_patients_ts_test_' + sampling + '_' + cov_type + '.feather'))
    all_scores_nca = feather.read_dataframe(os.path.join(path_res, 'all_patients_no_cross_auto_optimum_test_' + sampling + '_' + cov_type + '.feather'))
    all_ts_nca = feather.read_dataframe(os.path.join(path_res, 'all_patients_no_cross_auto_ts_test_' + sampling + '_' + cov_type + '.feather'))
    if surr:
        all_scores_surr_ca = feather.read_dataframe(os.path.join(path_res, 'all_patients_surrogate_optimum_test_' + sampling + '_' + cov_type + '.feather'))
        all_scores_surr_nca = feather.read_dataframe(os.path.join(path_res, 'all_patients_no_cross_auto_surrogate_optimum_test_' + sampling + '_' + cov_type + '.feather'))
except FileNotFoundError as e:
    print(e)
    print("file not found")

all_scores = pd.concat([all_scores_ca, all_scores_nca])
all_ts = pd.concat([all_ts_ca, all_ts_nca])
all_scores['dataset'] = dataset
all_ts['dataset'] = dataset
if surr:
    all_scores_surr = pd.concat([all_scores_surr_ca, all_scores_surr_nca])
    all_scores_surr['dataset'] = dataset

# removing NA
all_scores = all_scores[all_scores['score_method']=='prop_time_corr_1']
all_scores = all_scores.dropna()
all_scores = all_scores.reset_index(drop=True)
all_ts = all_ts[all_ts['score_method']=='prop_time_corr_1']
all_ts = all_ts.dropna()
all_ts = all_ts.reset_index(drop=True)
if surr:
    all_scores_surr = all_scores_surr[all_scores_surr['score_method']=='prop_time_corr_1']
    all_scores_surr = all_scores_surr.dropna()
    all_scores_surr = all_scores_surr.reset_index(drop=True)

# merging patients with two blocks
for index, iRow in all_scores[all_scores['patient'].str.contains('2')].iterrows():
    jRow = all_scores[(all_scores['patient']==iRow['patient'][:2]) &
                      (all_scores['cov_choice']==iRow['cov_choice']) & 
                      (all_scores['score_method']==iRow['score_method']) & 
                      (all_scores['dataset']==dataset)]
    if len(jRow)>0:
        all_scores.at[jRow.index.values[0], 'score'] = (jRow['score'] + iRow['score'])/2
        all_scores = all_scores.drop([iRow.name])


for ipatient in all_ts[all_ts['patient'].str.contains('2')]['patient'].unique():
    print(ipatient)
    for icov_choice in all_ts['cov_choice'].unique():
        for iscore_method in all_ts['score_method'].unique():
            iRows = all_ts[(all_ts['patient']==ipatient) & 
                           (all_ts['cov_choice']==icov_choice) & 
                           (all_ts['score_method']==iscore_method) & 
                           (all_ts['dataset']==dataset)]
            if len(iRows)>1:
                time_1 = all_ts[(all_ts['patient']==ipatient[:2]) & 
                                (all_ts['cov_choice']==icov_choice) & 
                                (all_ts['score_method']==iscore_method)]['time'].max()
                iRows.loc[:, 'patient']=iRows.loc[:, 'patient'].values[0][:2] 
                iRows.loc[:, 'time']=iRows.loc[:, 'time'] + time_1 
                all_ts = all_ts.append(iRows)
                all_ts = all_ts.drop(iRows.index)
            elif len(iRow)==1:
                print('****problem only one row*********')

if surr:
    for index, iRow in all_scores_surr[all_scores_surr['patient'].str.contains('2')].iterrows():
        jRow = all_scores_surr[(all_scores_surr['patient']==iRow['patient'][:2]) & 
                               (all_scores_surr['cov_choice']==iRow['cov_choice']) & 
                               (all_scores_surr['score_method']==iRow['score_method']) & 
                               (all_scores_surr['i_surr']==iRow['i_surr']) & 
                               (all_scores['dataset']==dataset)]
        if len(jRow)>0:
            all_scores_surr.at[jRow.index.values[0], 'score'] = (jRow['score'] + iRow['score'])/2
            all_scores_surr = all_scores_surr.drop([iRow.name])


# Brier skill score per patient
lBSS_all_clim, lBSS_all_rand, lpatients_all, lcov_choice_all, ldataset_all, lsig_all = [], [], [], [], [], []
for icov_choice in lcov_choice:
    print(icov_choice)
    lpatients =  np.unique(all_ts[all_ts['dataset']==dataset]['patient'])
    pvals = []
    for ipatient, patient in enumerate(lpatients):
        # select the time series of interest
        ts = all_ts[(all_ts['patient']==patient) & 
                    (all_ts['cov_choice']==icov_choice) & 
                    (all_ts['score_method']=='prop_time_corr_1') & 
                    (all_ts['dataset']==dataset)]
        # select the test seizures
        y = ts['test_Sz'].values.astype('int')
        # get the predicted probabilities
        pk = ts['preds_test'].values
        # threshold to 1
        pk[pk>1] = 1
        # Brier score
        BS = np.mean((pk-y)**2)
        # Brier score after rando; shuffling
        BSrand = []
        for j_ran in range(1000):
            BSrand.append(np.mean((np.random.choice(pk, pk.shape[0])-y)**2))
        BSrand = np.mean(BSrand)
        BSSrand = 1-(BS/BSrand)
        lBSS_all_rand.append(BSSrand)
        lpatients_all.append(patient)
        lcov_choice_all.append(icov_choice)
        ldataset_all.append(dataset)

        # find AUC pvals to find subjects with IoC
        score_true = all_scores[(all_scores['patient']==patient) & 
                                (all_scores['cov_choice']==icov_choice) & 
                                (all_scores['dataset']==dataset)]
        x_true = score_true['score'].values
        score_surr = all_scores_surr[(all_scores_surr['cov_choice']==icov_choice) & 
                                     (all_scores_surr['dataset']==dataset)]
        score_surr = score_surr[score_surr['score_method']=='prop_time_corr_1']
        score_surr = score_surr[score_surr['patient']==patient]
        y_surr = score_surr['score'].values
        if y_surr.shape[0]<200:
            print('for patient ' + patient + ' there are ' + str(y_surr.shape[0]) + ' surrogates')
        pval = (np.sum(y_surr>x_true) + 1)/(y_surr.shape[0])
        pvals.append(pval)

        # FDR correction        
        q = 0.05
        ps = np.sort(np.array(pvals))
        V = ps.shape[0]
        I = np.arange(1, V+1)
        cVID = 1.
        cVN = np.sum(1/I)

        if np.nonzero(ps<=I/V*q/cVID)[0].shape!=(0,):
            pID = ps[np.max(np.nonzero(ps<=I/V*q/cVID))]
            lsig_all.extend((np.array(pvals)<pID+0.0001).tolist())
        else:
            lsig_all.extend([False]*len(pvals))
            print('not significant anywhere')

dfBSSrand = pd.DataFrame.from_items([('patient', lpatients_all), 
                                 ('BSS', lBSS_all_rand), 
                                 ('cov_choice', lcov_choice_all), 
                                 ('dataset', ldataset_all),
                                 ('significant', lsig_all)])

for cov_choice in lcov_choice:
    curr_dat = dfBSSrand[(dfBSSrand['cov_choice']==cov_choice) & (dfBSSrand['dataset']==dataset)]
    print(cov_choice + ' ' + dataset + ' all : ' + str(curr_dat['BSS'].median()) + ' ' + str(np.sum(curr_dat['BSS']>0)/curr_dat['BSS'].shape[0]))
    print(cov_choice + ' ' + dataset + ' all : ' + str(curr_dat['BSS'].quantile(.75)) + ' ' + str(np.sum(curr_dat['BSS']>0)/curr_dat['BSS'].shape[0]))
    print(cov_choice + ' ' + dataset + ' all : ' + str(curr_dat['BSS'].quantile(.25)) + ' ' + str(np.sum(curr_dat['BSS']>0)/curr_dat['BSS'].shape[0]))
    print(cov_choice + ' ' + dataset + ' significant: ' + str(curr_dat[curr_dat['significant']]['BSS'].median()))
    print(cov_choice + ' ' + dataset + ' significant: ' + str(curr_dat[curr_dat['significant']]['BSS'].quantile(.75)))
    print(cov_choice + ' ' + dataset + ' significant: ' + str(curr_dat[curr_dat['significant']]['BSS'].quantile(.25)))
