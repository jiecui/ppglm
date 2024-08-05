# Author: Timothée Proix
# License: GPL-3.0-only

library(R.utils)
library(feather)
library(reshape2)
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

cov.choice.list = c('only_cross', 'only_auto')

params = list()
params['cov.type'] = "acausal4"
params['sampling'] = "day" # "day", "hour"
params['surrogates'] = 0
score_list = list('prop_time_corr_1')

if (params['sampling']=='hour'){
  params['cross.history.IEA.eval'] = 11
} else if (params['sampling']=='day'){
  params['cross.history.IEA.eval'] = 6
}

if (file.exists(paste0(paths[['path.res']], "all_patients_ts_test_", params[['sampling']], "_", params[['cov.type']], ".feather"))){

  # load existing dataframe
  scores.validation = read_feather(paste0(paths[['path.res']], "all_patients_score_validation_", params[['sampling']], "_", params[['cov.type']], ".feather"))
  optimal.history.validation = read_feather(paste0(paths[['path.res']], "all_patients_optimum_validation_", params[['sampling']], "_", params[['cov.type']], ".feather"))
  optimal.history.test = read_feather(paste0(paths[['path.res']], "all_patients_optimum_test_", params[['sampling']], "_", params[['cov.type']], ".feather"))
  ts.val = read_feather(paste0(paths[['path.res']], "all_patients_ts_val_", params[['sampling']], "_", params[['cov.type']], ".feather"))
  ts.test = read_feather(paste0(paths[['path.res']], "all_patients_ts_test_", params[['sampling']], "_", params[['cov.type']], ".feather"))
  curve.auc = read_feather(paste0(paths[['path.res']], "all_patients_curve_test_", params[['sampling']], "_", params[['cov.type']], ".feather"))

} else {
  # create data frame
  scores.validation <- data.frame(
    patient = c(NA),
    score_method = c(NA), 
    score.Var1 = c(NA), 
    score.Var2 = c(NA), 
    score.value = c(NA),
    cov_choice = c(NA)
  )
  
  optimal.history.validation <- data.frame(
    patient = c(NA),
    score_method = c(NA), 
    score = c(NA),
    auto_history_Sz_opt = c(NA),
    cross_history_IEA_opt = c(NA),
    cov_choice = c(NA)
  )
  
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

  ts.val <- data.frame(
    patient = c(NA),
    score_method = c(NA), 
    preds_val = c(NA),
    val_Sz = c(NA),
    time = c(NA),
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
}

env = environment()
for (cov.choice in cov.choice.list){
  params['cov.choice'] = cov.choice
  print(params[['cov.choice']])
  for (patient in patient.list){
    ncfname = paste0(paths[['path.res']], '/PSE/', patient, "_predictions_", params[['cov.choice']], "_", params[['sampling']], "_", params[['cov.type']], ".nc")
    if ((params[['cov.choice']]=='only_auto') || (params[['cov.choice']]=='only_cross')){
      ncfname = paste0(paths[['path.res']], '/PSE/', patient, "_predictions_only_cross_auto_", params[['sampling']], "_", params[['cov.type']], ".nc")
    }

    if (file.exists(ncfname)){
      print('###')
      ncin = nc_open(ncfname)

      # dim of predictions.labels: (2, auto.history, cross.history, folds, ts)
      predictions.labels = ncvar_get(ncin, "predictions.labels")

      for (score in score_list){
        print(score)
        if (dim(subset(ts.test, patient==get('patient', env) & score_method==get('score', env) & cov_choice==get('cov.choice', env)))[1]==0){ 
          m.score = matrix(data=0, nrow=dim(predictions.labels)[2], ncol=dim(predictions.labels)[3])
          for (i.auto.history.Sz in c(1:dim(predictions.labels)[2])){
            for (j.cross.history.IEA in c(1:dim(predictions.labels)[3])){
              nb.fold = 5
              for (k.fold in c(1:nb.fold)){
                if (!is.na(sum(predictions.labels[1, i.auto.history.Sz, j.cross.history.IEA, k.fold,]))){
                  if ((sum(predictions.labels[2, i.auto.history.Sz, j.cross.history.IEA, k.fold,])!=0) && (sum(predictions.labels[1, i.auto.history.Sz, j.cross.history.IEA, k.fold,])!=0)){
                    m.score[i.auto.history.Sz, j.cross.history.IEA] = m.score[i.auto.history.Sz, j.cross.history.IEA] + scoreAuc(predictions.labels[2, i.auto.history.Sz, j.cross.history.IEA, k.fold,], predictions.labels[1, i.auto.history.Sz, j.cross.history.IEA, k.fold,], score)
                  } else if (i.auto.history.Sz==1 && j.cross.history.IEA==1){
                    nb.fold = nb.fold - 1
                  } else {
                    print("no seizures in fold")
                    print(patient)
                    print(i.auto.history.Sz)
                    print(j.cross.history.IEA)
                    print(k.fold)
                    nb.fold = nb.fold - 1
                  }
                } else {
                  print("Fold couldn't be predicted")
                  print(patient)
                  print(i.auto.history.Sz)
                  print(j.cross.history.IEA)
                  print(k.fold)
                }
              }
              if (nb.fold>0){
                m.score[i.auto.history.Sz, j.cross.history.IEA] = 1/nb.fold * m.score[i.auto.history.Sz, j.cross.history.IEA]
              } else {
                m.score[i.auto.history.Sz, j.cross.history.IEA] = 0.5
              }
            }
          }
          # save the matrices for the figures
          scores.validation.add <- data.frame(
            patient = patient,
            score_method = score, 
            score = melt(m.score),
            cov_choice = params[['cov.choice']],
            stringsAsFactors=FALSE
          )

          scores.validation = rbind(scores.validation, scores.validation.add)

          # find optimal
          list(score.validation, auto.history.Sz.max, cross.history.IEA.max) := findMax(m.score, params)
          # add to the dataframe
          optimal.history.validation.add = data.frame(
            patient = patient,
            score_method = score, 
            score = score.validation,
            auto_history_Sz_opt = auto.history.Sz.max,
            cross_history_IEA_opt = cross.history.IEA.max,
            cov_choice = params[['cov.choice']]
          )

          optimal.history.validation = rbind(optimal.history.validation, optimal.history.validation.add)

          # run the test model for each score
          params['auto.history.Sz'] = auto.history.Sz.max - 1
          params['cross.history.IEA'] = cross.history.IEA.max - 1 
          #list(preds.test, test.Sz, cov.test.IEA) := predictionTest(paths, patient, params)
          # error catching, if you want to run in safe mode
          possibleError = tryCatch({
            list(preds.test, test.Sz, cov.test.IEA) := predictionTest(paths, patient, params)
          }, error = function(cond){
            browser()
            message("Prediction on test did not work")
          })
          if(inherits(possibleError, "error")) next

          score.test = scoreAuc(test.Sz, preds.test$pred, score)    

          #for figure box plots
          optimal.history.test.add = data.frame(
            patient = patient,
            score_method = score,
            score = score.test,
            auto_history_Sz_opt = auto.history.Sz.max,
            cross_history_IEA_opt = cross.history.IEA.max,
            cov_choice = params[['cov.choice']]
          )
          optimal.history.test = rbind(optimal.history.test, optimal.history.test.add)

          # for figure time series
          ts.val.add <- data.frame(
            patient = patient,
            score_method = score,
            preds_val = as.vector(predictions.labels[1, auto.history.Sz.max, cross.history.IEA.max,,]),
            val_Sz = as.vector(predictions.labels[2, auto.history.Sz.max, cross.history.IEA.max,,]),
            time = c(1: length(as.vector(predictions.labels[1, auto.history.Sz.max, cross.history.IEA.max,,]))),
            cov_choice = cov.choice
          )
          ts.val = rbind(ts.val, ts.val.add)

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

          # save dataframe for later access in python
          write_feather(scores.validation, paste0(paths[['path.res']], "all_patients_score_validation_", params[['sampling']], "_", params[['cov.type']], ".feather"))
          write_feather(optimal.history.validation, paste0(paths[['path.res']], "all_patients_optimum_validation_", params[['sampling']], "_", params[['cov.type']], ".feather"))
          write_feather(optimal.history.test, paste0(paths[['path.res']], "all_patients_optimum_test_", params[['sampling']], "_", params[['cov.type']], ".feather"))
          write_feather(ts.test, paste0(paths[['path.res']], "all_patients_ts_test_", params[['sampling']], "_", params[['cov.type']], ".feather"))
          write_feather(ts.val, paste0(paths[['path.res']], "all_patients_ts_val_", params[['sampling']], "_", params[['cov.type']], ".feather"))
          write_feather(curve.auc, paste0(paths[['path.res']], "all_patients_curve_test_", params[['sampling']], "_", params[['cov.type']], ".feather"))
        } else {
          print(paste('already computed', cov.choice, patient, score))
        }
      }
    } else { print(ncfname)}
  }
}
