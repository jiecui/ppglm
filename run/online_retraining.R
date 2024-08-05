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

cov.choice.list = c('phase_multi')

params = list()
params['cov.type'] = "acausal4"
params['sampling'] = "day" # "day", "hour"
score_list = list('prop_time_corr_1')

if (params['sampling']=='hour'){
  params['cross.history.IEA.eval'] = 11
  retraining.period.list = c(-1, 1, 12, 24, 24*10) # -1 means every seizure
  cut.training = 60*8*24
  cov.choice.list = c('multivariate')
} else if (params['sampling']=='day'){
  params['cross.history.IEA.eval'] = 6
  retraining.period.list = c(-1, 1, 10, 50) # -1 means every seizure
  cut.training = 60*8 
  cov.choice.list = c('phase_multi')
}

ts.test = read_feather(paste0(paths[['path.res']], "all_patients_no_cross_auto_ts_test_", params[['sampling']], "_", params[['cov.type']], ".feather"))

if (file.exists(paste0(paths[['path.res']], "all_patients_online_retraining_", params[['sampling']], "_", params[['cov.type']], ".feather"))){
  # load existing dataframe
  online.retraining = read_feather(paste0(paths[['path.res']], "all_patients_online_retraining_", params[['sampling']], "_", params[['cov.type']], ".feather"))
  ts.retraining = read_feather(paste0(paths[['path.res']], "all_patients_online_retraining_ts_", params[['sampling']], "_", params[['cov.type']], ".feather"))

} else {
  # create data frame
  online.retraining <- data.frame(
    patient = c(NA),
    score_method = c(NA), 
    retraining_period = c(NA),
    score = c(NA),
    cov_choice = c(NA)
  )

  ts.retraining <- data.frame(
    patient = c(NA),
    score_method = c(NA), 
    retraining_period = c(NA),
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
      if (dim(subset(online.retraining, patient==get('patient', env) & score_method==get('score', env) & cov_choice==get('cov.choice', env)))[1]==0){ 

        ts.test.idx = subset(ts.test, patient==get('patient', env) & score_method==get('score', env) & cov_choice==get('cov.choice', env))
        test.Sz.ref = ts.test.idx[['test_Sz']]
        length.training = min(floor((length(test.Sz.ref) - 3)*(6/4)), 60*8) # maximum length of the training data, 60% or 60*8

        for (retraining.period in retraining.period.list){
          print(paste('retraining period', retraining.period))
          preds.test.retrain = c()
          test.Sz.retrain = c()

          if (retraining.period>0){
            seq.retraining =  seq(from=0, to = length(test.Sz.ref), by = retraining.period)
          } else if (retraining.period==-1){
            seq.retraining = c(0, which(test.Sz.ref==1))
          }
          if (sum(seq.retraining==length(test.Sz.ref))==0){
            seq.retraining = c(seq.retraining, length(test.Sz.ref))
          }

          i.cursor = 1;
          for (cursor.retraining in seq.retraining[1:(length(seq.retraining)-1)]){
            print(paste('cursor.retraining', cursor.retraining))

            # run the test model for different length of the training period
            params['auto.history.Sz'] = 0
            params['cross.history.IEA'] = 0
            params['length.train'] = length.training + cursor.retraining
            params['length.test'] = diff(seq.retraining)[i.cursor]
            i.cursor = i.cursor + 1

            # error catching, if you want to run in safe mode
            possibleError = tryCatch(
              list(preds.test, test.Sz, cov.test.IEA) := predictionTest(paths, patient, params),
              error = function(e) e)

            if(inherits(possibleError, "error")){
              # the last period might be out of bound per the divisors of length(test.Sz.ref)
              print('Prediction on test did not work')
              next
            }
            preds.test.retrain[(cursor.retraining+1):(cursor.retraining + params[['length.test']])] = preds.test$pred
            test.Sz.retrain[(cursor.retraining+1):(cursor.retraining + params[['length.test']])] = test.Sz
          }

          score.test = scoreAuc(test.Sz.retrain, preds.test.retrain, score)

          # store the score
          online.retraining.add = data.frame(
            patient = patient,
            score_method = score,
            score = score.test,
            retraining_period = retraining.period,
            cov_choice = params[['cov.choice']]
          )
          online.retraining = rbind(online.retraining, online.retraining.add)

          # and the ts test
          ts.retraining.add <- data.frame(
            patient = patient,
            score_method = score,
            retraining_period = retraining.period,
            preds_test = preds.test.retrain,
            test_Sz = test.Sz.retrain,
            time = c(1: length(test.Sz.retrain)),
            cov_choice = params[['cov.choice']]
          )
          ts.retraining = rbind(ts.retraining, ts.retraining.add)

          # save dataframe for later access in python
          write_feather(online.retraining, paste0(paths[['path.res']], "all_patients_online_retraining_", params[['sampling']], "_", params[['cov.type']], ".feather"))
          write_feather(ts.retraining, paste0(paths[['path.res']], "all_patients_online_retraining_ts_", params[['sampling']], "_", params[['cov.type']], ".feather"))
          
        }
      } else {
        print(paste('already computed', cov.choice, patient, score))
      }
    }

  }
}
