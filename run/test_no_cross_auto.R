# Author: Timothée Proix
# License: GPL-3.0-only

library(R.utils)
library(feather)
sourceDirectory('../src')
sourceDirectory('../utils')

rootData = 'YourPath'
rootProcess = 'YourPath'
rootResults = 'YourPath'

dataset = 'UCSF_patients'

paths = list()
paths['path.data'] = paste0(rootData, '/forecasting_over_days/', dataset, '/patients/')
paths['path.process'] = paste0(rootProcess, '/forecasting_over_days/', dataset, '/')
paths['path.res'] = paste0(rootResults, '/forecasting_over_days/', dataset, '/')

patient.list.temp = read.csv(paste0(paths[['path.data']], 'patients.csv'), sep=",", header=F)
patient.list = as.character(patient.list.temp$V1)

params = list()
params['cov.type'] = "acausal4"
params['sampling'] = "day" # "day", "hour""
score_list = list('prop_time_corr_1')

if (params[['sampling']]=='hour'){
  params['cross.history.IEA.eval'] = 11
  cov.choice.list = c('circ_hist', 'phase_circ', 'phase_multi', 'phase_circ_multi', 'multivariate')
} else if (params[['sampling']]=='day'){
  params['cross.history.IEA.eval'] = 6
  cov.choice.list = c('week_hist', 'phase_multi')
}

if (file.exists(paste0(paths[['path.res']], "all_patients_no_cross_auto_ts_test_", params[['sampling']], "_", params[['cov.type']], ".feather"))){

  # load existing dataframe
  optimal.history.test = read_feather(paste0(paths[['path.res']], "all_patients_no_cross_auto_optimum_test_", params[['sampling']], "_", params[['cov.type']], ".feather"))
  ts.test = read_feather(paste0(paths[['path.res']], "all_patients_no_cross_auto_ts_test_", params[['sampling']], "_", params[['cov.type']], ".feather"))
  coeffs = read_feather(paste0(paths[['path.res']], "all_patients_no_cross_auto_coeffs_", params[['sampling']], "_", params[['cov.type']], ".feather"))
  curve.auc = read_feather(paste0(paths[['path.res']], "all_patients_no_cross_auto_curve_test_", params[['sampling']], "_", params[['cov.type']], ".feather"))

} else {
  # create data frame 
  optimal.history.test <- data.frame(
    patient = c(NA),
    score_method = c(NA), 
    score = c(NA),
    auto_history_Sz_opt = c(NA),
    cross_history_IEA_opt = c(NA),
    cov_choice = c(NA)
  )
  
  curve.auc<- data.frame(
    patient = c(NA),
    score_method = c(NA), 
    y = c(NA),
    x = c(NA),
    cov_choice = c(NA)
  )
  
  ts.test <- data.frame(
    patient = c(NA),
    score_method = c(NA), 
    preds_test = c(NA),
    test_Sz = c(NA),
    cov_test_IEA = c(NA),
    time = c(NA),
    cov_choice = c(NA)
  )

   coeffs <- data.frame(
    patient = c(NA),
    score_method = c(NA), 
    coeffs = c(NA),
    cov_choice = c(NA)
  )
}

env = environment()
for (cov.choice in cov.choice.list){
  params['cov.choice'] = cov.choice
  for (patient in patient.list){
    print('###')

    print(patient)
    print(cov.choice)

    if (dim(subset(ts.test, patient==get('patient', env) & score_method=='prop_time_corr_1' & cov_choice==get('cov.choice', env)))[1]==0){

      # run the test model for each score
      params['auto.history.Sz'] = 0
      params['cross.history.IEA'] = 0
     
      # error catching, if you want to run in safe mode
      possibleError = tryCatch({
        list(preds.test, test.Sz, cov.test.IEA, Sz.with.IEA.fit) := predictionTest(paths, patient, params)
      }, error = function(cond){
        browser()
        message("Prediction on test did not work")
      })
      if(inherits(possibleError, "error")) next

      for (score in score_list){
        print(score)
        score.test = scoreAuc(test.Sz, preds.test$pred, score)    

        #for figure box plots
        optimal.history.test.add = data.frame(
          patient = patient,
          score_method = score,
          score = score.test,
          auto_history_Sz_opt = 0,
          cross_history_IEA_opt = 0,
          cov_choice = params[['cov.choice']]
        )
        optimal.history.test = rbind(optimal.history.test, optimal.history.test.add)

        # for figure time series
        ts.test.add <- data.frame(
          patient = patient,
          score_method = score,
          preds_test = preds.test$pred,
          test_Sz = test.Sz,
          cov_test_IEA = cov.test.IEA[,1],
          time = c(1: length(test.Sz)),
          cov_choice = cov.choice
        )
        ts.test = rbind(ts.test, ts.test.add)


        # for figure curve of prop time
        list(y, x) := scoreCurve(test.Sz, preds.test$pred, score)
        curve.auc.add <- data.frame(
          patient = patient,
          score_method = score, 
          y = y,
          x = x,
          cov_choice = cov.choice
        )
        curve.auc = rbind(curve.auc, curve.auc.add)


        coeffs.add <- data.frame(
          patient = patient,
          score_method = score, 
          coeffs = Sz.with.IEA.fit$coefficients,
          cov_choice = cov.choice
        )
        coeffs = rbind(coeffs, coeffs.add)

        # save dataframe for later access in python
        write_feather(optimal.history.test, paste0(paths[['path.res']], "all_patients_no_cross_auto_optimum_test_", params[['sampling']], "_", params[['cov.type']], ".feather"))
        write_feather(ts.test, paste0(paths[['path.res']], "all_patients_no_cross_auto_ts_test_", params[['sampling']], "_", params[['cov.type']], ".feather"))
        write_feather(coeffs, paste0(paths[['path.res']], "all_patients_no_cross_auto_coeffs_", params[['sampling']], "_", params[['cov.type']], ".feather"))
        write_feather(curve.auc, paste0(paths[['path.res']], "all_patients_no_cross_auto_curve_test_", params[['sampling']], "_", params[['cov.type']], ".feather"))
      }
    } else {
      print(paste('already computed', cov.choice, patient))
    }
  }
}
