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

optimal.history.test = read_feather(paste0(paths[['path.res']], "all_patients_no_cross_auto_optimum_test_", params[['sampling']], "_", params[['cov.type']], ".feather"))

if (file.exists(paste0(paths[['path.res']], "all_patients_ahead_perf_", params[['sampling']], "_", params[['cov.type']], ".feather"))){
  # load existing dataframe
  ahead.perf = read_feather(paste0(paths[['path.res']], "all_patients_ahead_perf_", params[['sampling']], "_", params[['cov.type']],  ".feather"))
  ts.ahead = read_feather(paste0(paths[['path.res']], "all_patients_ahead_ts_", params[['sampling']], "_", params[['cov.type']], ".feather"))

} else {
  # create data frame
  ahead.perf <- data.frame(
    patient = c(NA),
    score_method = c(NA), 
    n_ahead = c(NA),
    score = c(NA),
    cov_choice = c(NA)
  )

  ts.ahead <- data.frame(
    patient = c(NA),
    score_method = c(NA), 
    n_ahead = c(NA),
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

      if (dim(subset(ahead.perf, patient==get('patient', env) & score_method==get('score', env) & cov_choice==get('cov.choice', env)))[1]==0){ 

        params['auto.history.Sz'] = 0
        params['cross.history.IEA'] = 0

        # One step ahead prediction
        possibleError = tryCatch(
          list(preds.test.ref, test.Sz.ref, cov.test.IEA.ref) := predictionTest(paths, patient, params),
          error = function(e) e)

        if(inherits(possibleError, "error")){
          print('Prediction on test did not work')
        }
        
        score.test.ref = scoreAuc(test.Sz.ref, preds.test.ref$pred, score)
         
        # check the score is as expected
        score.test.check = subset(optimal.history.test, patient==get('patient', env) & score_method==get('score', env) & cov_choice==get('cov.choice', env))
        if (score.test.ref!=score.test.check$score){
          browser()
        }

        # run the test model for different number of steps ahead
        params.ahead = params
        for (n.ahead in c(1:10)){
          params.ahead['n.ahead'] = n.ahead # number of points ahead for prediction
          preds.test.ahead = c()
          test.Sz.ahead = c()
          for (length.test in c(1:(length(test.Sz.ref) - n.ahead + 1))){
            print(n.ahead)
            print(length.test)
            params.ahead['length.test'] = length.test
            # error catching, if you want to run in safe mode
            possibleError = tryCatch(
              list(preds.test, test.Sz, cov.test.IEA) := predictionTest(paths, patient, params.ahead),
              error = function(e) e)

            if(inherits(possibleError, "error")){
              print('Prediction on test did not work')
              next
            }
             
            preds.test.ahead[length.test] = preds.test$pred[length(preds.test$pred)]
            # we need to use the reference because the test.Sz stops at n.ahead=1
            test.Sz.ahead[length.test] = test.Sz.ref[length.test - 1 + n.ahead]
          }
          score.test = scoreAuc(test.Sz.ahead, preds.test.ahead, score)    

          # store the score
          ahead.perf.add = data.frame(
            patient = patient,
            score_method = score,
            n_ahead = params.ahead[['n.ahead']],
            score = score.test,
            cov_choice = params.ahead[['cov.choice']]
          )
          ahead.perf = rbind(ahead.perf, ahead.perf.add)

          # and the ts test
          ts.ahead.add <- data.frame(
            patient = patient,
            score_method = score,
            n_ahead = params.ahead[['n.ahead']],
            preds_test = preds.test.ahead,
            test_Sz = test.Sz.ahead,
            time = c(1: length(test.Sz.ahead)),
            cov_choice =  params.ahead[['cov.choice']]
          )
          ts.ahead = rbind(ts.ahead, ts.ahead.add)

        }
        # save dataframe for later access in python
        write_feather(ahead.perf, paste0(paths[['path.res']], "all_patients_ahead_perf_", params[['sampling']], "_", params[['cov.type']], ".feather"))
        write_feather(ts.ahead, paste0(paths[['path.res']], "all_patients_ahead_ts_", params[['sampling']], "_", params[['cov.type']], ".feather"))

      } else {
        print(paste('already computed', cov.choice, patient, score))
      }
    }
  }
}