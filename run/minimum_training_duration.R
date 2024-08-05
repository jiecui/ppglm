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
params['sampling'] = "day" # "day", "hour"
score_list = list('prop_time_corr_1')


if (params['sampling']=='hour'){
  params['cross.history.IEA.eval'] = 11
  cov.choice.list = c('multivariate')
} else if (params['sampling']=='day'){
  params['cross.history.IEA.eval'] = 6
  cov.choice.list = c('phase_multi')
}

ts.test = read_feather(paste0(paths[['path.res']], "all_patients_no_cross_auto_ts_test_", params[['sampling']], "_", params[['cov.type']], ".feather"))

if (file.exists(paste0(paths[['path.res']], "all_patients_minimal_training_", params[['sampling']], "_", params[['cov.type']], ".feather"))){
  # load existing dataframe
  minimal.training = read_feather(paste0(paths[['path.res']], "all_patients_minimal_training_", params[['sampling']], "_", params[['cov.type']], ".feather"))
  ts.minimal = read_feather(paste0(paths[['path.res']], "all_patients_minimal_training_ts_", params[['sampling']], "_", params[['cov.type']], ".feather"))
} else {
  # create data frame
  minimal.training <- data.frame(
    patient = c(NA),
    score_method = c(NA), 
    training_period = c(NA),
    length_training = c(NA),  
    score = c(NA),
    cov_choice = c(NA)
  )

  ts.minimal <- data.frame(
    patient = c(NA),
    score_method = c(NA), 
    training_period = c(NA),
    length_training = c(NA),  
    preds_test = c(NA),
    test_Sz = c(NA),
    time = c(NA),
    cov_choice = c(NA)
  )

}

env = environment()
for (cov.choice in cov.choice.list){
  params['cov.choice'] = cov.choice
  print(params[['cov.choice']])
  for (patient in patient.list){

    print('###')
   
    for (score in score_list){
      print(score)
      if (dim(subset(minimal.training, patient==get('patient', env) & score_method==get('score', env) & cov_choice==get('cov.choice', env)))[1]==0){ 

        ts.test.idx = subset(ts.test, patient==get('patient', env) & score_method==get('score', env) & cov_choice==get('cov.choice', env))
        test.Sz = ts.test.idx[['test_Sz']]
        length.training = floor((length(test.Sz) - 3)*(6/4)) # maximum length of the training data, 60%

        if (params['sampling']=='hour'){
          max.training.length = 480*24
          sample.rate = 24
        } else if (params['sampling']=='day'){
          max.training.length = 480
          sample.rate = 1
        }

        for (cut.training in seq((params[['cross.history.IEA.eval']] + 1), min(max.training.length, (length.training-1)), sample.rate)){
          print(cut.training)

          # run the test model for different length of the training period
          params['auto.history.Sz'] = 0
          params['cross.history.IEA'] = 0
          params['cut.training'] = cut.training

          # error catching, if you want to run in safe mode
          possibleError = tryCatch(
            list(preds.test, test.Sz, cov.test.IEA) := predictionTest(paths, patient, params),
            error = function(e) e)

          if(inherits(possibleError, "error")){
            print('Prediction on test did not work')
            next
          }

          score.test = scoreAuc(test.Sz, preds.test$pred, score)    
          print(score.test)
          # store the score
          minimal.training.add = data.frame(
            patient = patient,
            score_method = score,
            score = score.test,
            training_period = cut.training,
            length_training = length.training,
            cov_choice = params[['cov.choice']]
          )
          minimal.training = rbind(minimal.training, minimal.training.add)

          # and the ts test
          ts.minimal.add <- data.frame(
            patient = patient,
            score_method = score,
            training_period = cut.training,
            length_training = length.training,
            preds_test = preds.test$pred,
            test_Sz = test.Sz,
            time = c(1: length(test.Sz)),
            cov_choice =  params[['cov.choice']]
          )
          ts.minimal = rbind(ts.minimal, ts.minimal.add)

          # save dataframe for later access in python
          write_feather(minimal.training, paste0(paths[['path.res']], "all_patients_minimal_training_", params[['sampling']], "_", params[['cov.type']], ".feather"))
          write_feather(ts.minimal, paste0(paths[['path.res']], "all_patients_minimal_training_ts_", params[['sampling']], "_", params[['cov.type']], ".feather"))
        }
      } else {
        print(paste('already computed', cov.choice, patient, score))
      }
    }
  }
}
