# Author: Timothée Proix
# License: GPL-3.0-only

import numpy as np
import seaborn as sns
import pandas as pd
import os
import feather
import pdb
import scipy.special as ss


# Snyder and Mormann stats
def FB(k, n, p):
    res = 0
    for j in np.arange(k+1):
        res += fB(j, n, p)
    return res


def fB(k, n, p):
    # print(p)
    # print(k)
    # print(n)
    return ss.comb(n, k) * p**k * (1-p)**(n-k)

rootResults = 'YourPath'
dataset = 'UCSF_patients' 
path_res = os.path.join(rootResults, 'forecasting_over_days', dataset)

sampling = 'day' # day, hour
cov_type = 'acausal4'
output = 2 # 1: median AUC,  p-values surrogates, 3: p-values for optimal threshold

if sampling=='day':
    lcov_choice = ['only_auto', 'only_cross', 'week_hist', 'phase_multi']
    lcov_choice_nca = ['week_hist', 'phase_multi']
elif sampling=='hour':
    lcov_choice = ['only_auto', 'circ_hist', 'only_cross', 'phase_circ', 'phase_multi', 'multivariate']
    lcov_choice_nca = ['circ_hist', 'phase_circ', 'phase_multi', 'multivariate']

try:
    all_scores_ca = feather.read_dataframe(os.path.join(path_res, 'all_patients_optimum_test_' + sampling + '_' + cov_type + '.feather'))
    all_scores_nca = feather.read_dataframe(os.path.join(path_res, 'all_patients_no_cross_auto_optimum_test_' + sampling + '_' + cov_type + '.feather'))
    all_curves_day_ca = feather.read_dataframe(os.path.join(path_res, 'all_patients_curve_test_'+ sampling+'_' + cov_type + '.feather'))
    all_curves_day_nca = feather.read_dataframe(os.path.join(path_res, 'all_patients_no_cross_auto_curve_test_'+ sampling+'_' + cov_type + '.feather'))
    all_ts_ca = feather.read_dataframe(os.path.join(path_res, 'all_patients_no_cross_auto_ts_test_' + sampling + '_' + cov_type + '.feather'))
    all_ts_nca = feather.read_dataframe(os.path.join(path_res, 'all_patients_ts_test_' + sampling + '_' + cov_type + '.feather'))
    if output==2:
        all_scores_surr_ca = feather.read_dataframe(os.path.join(path_res, 'all_patients_surrogate_optimum_test_' + sampling + '_' + cov_type + '.feather'))
        all_scores_surr_nca = feather.read_dataframe(os.path.join(path_res, 'all_patients_no_cross_auto_surrogate_optimum_test_' + sampling + '_' + cov_type + '.feather'))
except FileNotFoundError as e:
    print(e)
    print("file not found")

score_selection_ca = all_scores_ca.loc[all_scores_ca['cov_choice'].isin(['only_auto', 'only_cross'])].copy()
if output==2:
    score_selection_surr_ca = all_scores_surr_ca.loc[all_scores_surr_ca['cov_choice'].isin(['only_auto', 'only_cross'])].copy()
if sampling=='day':
    score_selection_nca = all_scores_nca.loc[all_scores_nca['cov_choice'].isin(lcov_choice_nca)].copy()
    if output==2:
        score_selection_surr_nca = all_scores_surr_nca.loc[all_scores_surr_nca['cov_choice'].isin(lcov_choice_nca)].copy()
elif sampling=='hour':
    score_selection_nca = all_scores_nca.loc[all_scores_nca['cov_choice'].isin(lcov_choice_nca)].copy()
    if output==2:
        score_selection_surr_nca = all_scores_surr_nca.loc[all_scores_surr_nca['cov_choice'].isin(lcov_choice_nca)].copy()

score_selection = pd.concat([score_selection_ca, score_selection_nca])
score_selection['dataset'] = dataset
all_curves_day = pd.concat([all_curves_day_nca, all_curves_day_ca])
all_curves_day['dataset'] = dataset
all_ts = pd.concat([all_ts_ca, all_ts_nca])
all_ts['dataset'] = dataset
if output==2:
    score_selection_surr = pd.concat([score_selection_surr_ca, score_selection_surr_nca])
    score_selection_surr['dataset'] = dataset

score_selection = score_selection.loc[score_selection['score_method']=='prop_time_corr_1']
all_curves_day = all_curves_day.loc[all_curves_day['score_method']=='prop_time_corr_1']
all_ts = all_ts.loc[all_ts['score_method']=='prop_time_corr_1']
if output==2:
    score_selection_surr = score_selection_surr.loc[score_selection_surr['score_method']=='prop_time_corr_1']

# removing NA and reset index to guarantee uniqueness
score_selection = score_selection.dropna()
score_selection = score_selection.reset_index(drop=True)
all_curves_day = all_curves_day.dropna()
all_curves_day = all_curves_day.reset_index(drop=True)
all_ts = all_ts.dropna()
all_ts = all_ts.reset_index(drop=True)
if output==2:
    score_selection_surr = score_selection_surr.dropna()
    score_selection_surr = score_selection_surr.reset_index(drop=True)

# merging patients with two blocks only for score AUC
for index, iRow in score_selection[score_selection['patient'].str.contains('2')].iterrows():
    jRow = score_selection[(score_selection['patient']==iRow['patient'][:2]) & 
                           (score_selection['cov_choice']==iRow['cov_choice']) & 
                           (score_selection['score_method']==iRow['score_method']) & 
                           (score_selection['dataset']==dataset)]
    if len(jRow)>0:
        score_selection.at[jRow.index.values[0], 'score'] = (jRow['score'] + iRow['score'])/2
        score_selection = score_selection.drop([iRow.name])

if output==2:
    for index, iRow in score_selection_surr[score_selection_surr['patient'].str.contains('2')].iterrows():
        jRow = score_selection_surr[(score_selection_surr['patient']==iRow['patient'][:2]) & 
                                    (score_selection_surr['cov_choice']==iRow['cov_choice']) & 
                                    (score_selection_surr['score_method']==iRow['score_method']) & 
                                    (score_selection_surr['i_surr']==iRow['i_surr']) & 
                                    (score_selection['dataset']==dataset)]
        if len(jRow)>0:
            score_selection_surr.at[jRow.index.values[0], 'score'] = (jRow['score'] + iRow['score'])/2
            score_selection_surr = score_selection_surr.drop([iRow.name])

## Find median AUC
if output==1:
    for cov_choice in lcov_choice:
        med = score_selection[(score_selection['cov_choice']==cov_choice) & 
                              (score_selection['dataset']==dataset)]['score'].median()
        iqr1 = score_selection[(score_selection['cov_choice']==cov_choice) & 
                               (score_selection['dataset']==dataset)]['score'].quantile(0.75)
        iqr2 = score_selection[(score_selection['cov_choice']==cov_choice) & 
                               (score_selection['dataset']==dataset)]['score'].quantile(0.25)
        print('Median ' + dataset + ' ' + cov_choice + ':' + str(med))
        print('Quartile ' + dataset + ' ' + cov_choice + ':' + str(iqr1))
        print('Quartile ' + dataset + ' ' + cov_choice + ':' + str(iqr2))

# surrogate stats
if output==2:
    for cov_choice in lcov_choice:
        print(cov_choice)

        lpatients =  np.unique(score_selection[score_selection['dataset']==dataset]['patient'])

        pvals = []
        for ipatient, patient in enumerate(lpatients):
            # print(str(ipatient) + ' ' + patient)
            score_true = score_selection[(score_selection['patient']==patient) & (score_selection['cov_choice']==cov_choice) & (score_selection['dataset']==dataset)]
            x_true = score_true['score'].values
            score_surr = score_selection_surr[(score_selection_surr['cov_choice']==cov_choice) & (score_selection_surr['dataset']==dataset)]
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

        print('For a threshold of: ' + str(q))
        if np.nonzero(ps<=I/V*q/cVID)[0].shape!=(0,):
            pID = ps[np.max(np.nonzero(ps<=I/V*q/cVID))]
            print(np.sum(np.array(pvals)<pID+0.0001))
            print(len(pvals))
            print(np.sum(np.array(pvals)<pID+0.0001)/(len(pvals)))
            AUC_curr = score_selection[(score_selection['cov_choice']==cov_choice) & (score_selection['dataset']==dataset)]['score'].values
            print('median AUC of significant: ' + str(np.median(AUC_curr[np.array(pvals)<pID+0.0001])))
            print('quantile AUC of significant: ' + str(np.quantile(AUC_curr[np.array(pvals)<pID+0.0001], 0.75)))
            print('quantile AUC of significant: ' + str(np.quantile(AUC_curr[np.array(pvals)<pID+0.0001], 0.25)))
            score_selection.loc[(score_selection['cov_choice']==cov_choice) & (score_selection['dataset']==dataset), 'IOC'] = np.array(pvals)<pID+0.0001
        else:
            score_selection.loc[(score_selection['cov_choice']==cov_choice) & (score_selection['dataset']==dataset), 'IOC'] = False
            print('not significant anywhere')

## print sensitivity and ptw for some thresholds and p-values for optimal threshold
if output==3:
    for cov_choice in lcov_choice:
        curves_day = all_curves_day.loc[all_curves_day['cov_choice']==cov_choice]
        ts = all_ts.loc[all_ts['cov_choice']==cov_choice]
        lpatients_all, lsens_ptw_all, ldataset_all, lval_all, lthreshold, ln_pred_all, lN_tot_all, lSnc_all = [], [], [], [], [], [], [], []
        lpatients =  np.unique(curves_day[curves_day['dataset']==dataset]['patient'])
        for ipatient, patient in enumerate(lpatients):
            # print(patient)
            # print(ipatient +1)
            curves_day_pat = curves_day.loc[(curves_day['patient']==patient) & 
                                            (curves_day['dataset']==dataset)]
            ts_IEA = ts[(ts['patient']==patient) & (ts['dataset']==dataset)]['preds_test'].values
            ts_sort =  np.concatenate([[0], np.sort(ts_IEA), [1]])
            ts_sort = np.sort(ts_IEA)
            ts_sz = ts[(ts['patient']==patient) & (ts['dataset']==dataset)]['test_Sz'].values.astype(bool)
            c_day = curves_day_pat[['x','y']].values
            idx_min_day = np.argmin(np.sum(([0,1] - c_day)**2, 1))
            if idx_min_day==0:
                idx_min_day = 1
            elif idx_min_day==ts_sort.shape[0] + 1:
                idx_min_day = ts_sort.shape[0]
            sens_day = c_day[idx_min_day, 1] 
            ptw_day = c_day[idx_min_day, 0] 
            th_opt = ts_sort[ts_sort.shape[0] - idx_min_day] + 0.0000001
            n_pred_opt = np.sum(ts_sz[ts_IEA>=th_opt])
            N_tot_opt = np.sum(ts_sz)
            Snc_opt = len(ts_sz[ts_IEA>=th_opt])/len(ts_sz)
            lpatients_all.extend([patient]*2)
            lsens_ptw_all.extend(['sens', 'ptw'])
            lval_all.extend([sens_day, ptw_day])
            ldataset_all.extend([dataset]*2)
            lthreshold.extend([th_opt, th_opt])
            ln_pred_all.extend([n_pred_opt, n_pred_opt])
            lN_tot_all.extend([N_tot_opt, N_tot_opt])
            lSnc_all.extend([Snc_opt, Snc_opt])

        dfVals = pd.DataFrame.from_items([('patient', lpatients_all), 
                                         ('sens_ptw', lsens_ptw_all), 
                                         ('val', lval_all), 
                                         ('dataset', ldataset_all), 
                                         ('threshold', lthreshold),
                                         ('n_pred', ln_pred_all), 
                                         ('N_tot', lN_tot_all),
                                         ('Snc', lSnc_all)])

        # merge patients with two blocks
        for index, iRow in dfVals[dfVals['patient'].str.contains('2')].iterrows():
            # print(index)
            if len(iRow)>0:
                jRow = dfVals[(dfVals['patient']==iRow['patient'][:2]) & 
                              (dfVals['sens_ptw']==iRow['sens_ptw']) & 
                              (dfVals['dataset']==iRow['dataset'])]
                dfVals.at[jRow.index.values[0], 'val'] = (jRow['val'] + iRow['val'])/2
                dfVals.at[jRow.index.values[0], 'threshold'] = (jRow['threshold'] + iRow['threshold'])/2
                dfVals.at[jRow.index.values[0], 'n_pred'] = (jRow['n_pred'] + iRow['n_pred'])/2
                dfVals.at[jRow.index.values[0], 'N_tot'] = (jRow['N_tot'] + iRow['N_tot'])/2
                dfVals.at[jRow.index.values[0], 'Snc'] = (jRow['Snc'] + iRow['Snc'])/2
                dfVals = dfVals.drop([iRow.name])



        lpatients =  np.unique(dfVals[dfVals['dataset']==dataset]['patient'])
        p_all = []
        for patient in lpatients:
            for (n, N, Snc) in dfVals[(dfVals['patient']==patient) & 
                                      (dfVals['dataset']==dataset) & 
                                      (dfVals['sens_ptw']=='sens')][['n_pred', 'N_tot', 'Snc']].itertuples(index=False):
                if n/N>=Snc:
                    # print((n, N, Snc))
                    p = 1 - FB(n-1, N, Snc) + FB(np.floor(2*N*Snc-n), N, Snc)
                else:
                    p = 1 - FB(np.ceil(2*N*Snc - n), N, Snc) + FB(n, N, Snc)
            p_all.append(p)

        q = 0.05
        ps = np.sort(np.array(p_all))
        V = ps.shape[0]
        I = np.arange(1, V+1)
        cVID = 1.
        cVN = np.sum(1/I)

        print('For a threshold of: ' + str(q) + ' with all p-values')
        print('*** cov choice: ' + cov_choice + ', data type: ' + dataset + ' ***')
        if np.nonzero(ps<=I/V*q/cVID)[0].shape!=(0,):
            pID = ps[np.max(np.nonzero(ps<=I/V*q/cVID))]
            print(np.nonzero(np.array(p_all)<pID+0.0001)[0].shape[0]/len(p_all))
            print(np.nonzero(np.array(p_all)<pID+0.0001)[0].shape[0])
        else:
            print('not significant anywhere')
