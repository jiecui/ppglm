# Author: Timothée Proix
# License: GPL-3.0-only

library(ncdf4)
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

params = list()
patient.list.temp = read.csv(paste0(paths[['path.data']], 'patients.csv'), sep=",", header=F)
patient.list = as.character(patient.list.temp$V1)

cov.choice.list = c('only_cross', 'only_auto')

params['cov.type'] = "acausal4"
params['sampling'] = "day" # "day", "hour"
params['surrogates'] = 1

if (params['sampling']=='hour'){
  params['cross.history.IEA.eval'] = 11
} else if (params['sampling']=='day'){
  params['cross.history.IEA.eval'] = 6
}

score_list = list('prop_time_corr_1')

env = environment()
for (i.surr in c(1:200)){
  print('i.surr')
  print(i.surr)
  for (patient in patient.list){

    if (file.exists(paste0(paths[['path.res']], "/surrogates/patients_", patient, "_surrogate_optimum_test_", params[['sampling']], "_", params[['cov.type']], ".feather"))){
      # load existing dataframe
      optimal.history.test = read_feather(paste0(paths[['path.res']], "/surrogates/patients_", patient, "_surrogate_optimum_test_", params[['sampling']], "_", params[['cov.type']], ".feather"))
    } else {
      # create data frame 
      optimal.history.test <- data.frame(
        patient = c(NA),
        score_method = c(NA), 
        score = c(NA),
        auto_history_Sz_opt = c(NA),
        cross_history_IEA_opt = c(NA),
        cov_choice = c(NA),
        i_surr = c(NA)
      )
    }

    for (cov.choice in cov.choice.list){
      params['cov.choice'] = cov.choice
      print(params[['cov.choice']])

      if (dim(subset(optimal.history.test, patient==get('patient', env) & score_method=='prop_time_corr_1' & cov_choice==get('cov.choice', env) & i_surr==get('i.surr', env)))[1]==0){
        # update the path for surrogate number i
        paths['path.process'] = paste0(rootProcess, '/forecasting_over_days/', dataset, '/surr_', toString(i.surr), '/') 
        paths['path.res'] = paste0(rootResults, '/forecasting_over_days/', dataset, '/surrogates', '/surr_', toString(i.surr), '/') 

        ncfname = paste0(paths[['path.res']], '/', patient, "_predictions_only_cross_auto_", params[['sampling']], "_", params[['cov.type']], ".nc")

        ncin = nc_open(ncfname)

        # dim of predictions.labels: (2, auto.history, cross.history, folds, ts)
        predictions.labels = ncvar_get(ncin, "predictions.labels")
        # if (params[['cov.choice']]=='only_auto'){
        #   dim(predictions.labels) <- c(dim(predictions.labels)[1], dim(predictions.labels)[2], 1, dim(predictions.labels)[3], dim(predictions.labels)[4])
        # } else if (params[['cov.choice']]=='only_cross') {
        #   dim(predictions.labels) <- c(dim(predictions.labels)[1], 1, dim(predictions.labels)[2], dim(predictions.labels)[3], dim(predictions.labels)[4])     
        # }

        for (score in score_list){ #
          print(score)

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

          # find optimal
          list(score.validation, auto.history.Sz.max, cross.history.IEA.max) := findMax(m.score, params)
            
          # run the test model for each score
          params['auto.history.Sz'] = auto.history.Sz.max - 1
          params['cross.history.IEA'] = cross.history.IEA.max - 1

          # run the test model for each score
          #list(preds.test, test.Sz, cov.test.IEA) := predictionTest(paths, patient, params)
          # error catching, if you want to run in safe mode

          possibleError = tryCatch({
            list(preds.test, test.Sz, cov.test.IEA) := predictionTest(paths, patient, params)
          }, error = function(cond){
            message("Prediction on test did not work")
          })
          if(inherits(possibleError, "error")) next

          score.test = scoreAuc(test.Sz, preds.test$pred, score)    
          #for figure box plots
          optimal.history.test.add = data.frame(
            patient = patient,
            score_method = score,
            score = score.test,
            auto_history_Sz_opt = params[['auto.history.Sz']],
            cross_history_IEA_opt = params[['cross.history.IEA']],
            cov_choice = params[['cov.choice']],
            i_surr = i.surr
          )
          optimal.history.test = rbind(optimal.history.test, optimal.history.test.add)

        }
      } else {
        print(paste('already computed', cov.choice, patient, i.surr))
      }
    }
    # save dataframe for later access in python
    write_feather(optimal.history.test, paste0(rootResults, '/forecasting_over_days/', dataset, "/surrogates/patients_", patient, "_surrogate_optimum_test_", params[['sampling']], "_", params[['cov.type']], ".feather"))
  }
}