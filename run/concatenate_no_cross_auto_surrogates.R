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

if (params[['sampling']]=='hour'){
  params['cross.history.IEA.eval'] = 11
  cov.choice.list = c('circ_hist', 'phase_circ', 'phase_multi', 'phase_circ_multi', 'multivariate')
} else if (params[['sampling']]=='day'){
  params['cross.history.IEA.eval'] = 6
  cov.choice.list = c('week_hist', 'phase_multi')
}

optimal.history.test <- data.frame(
  patient = c(NA),
  score_method = c(NA), 
  score = c(NA),
  auto_history_Sz_opt = c(NA),
  cross_history_IEA_opt = c(NA),
  cov_choice = c(NA),
  i_surr = c(NA)
)


for (patient in patient.list){
  print(patient)
  fname = paste0(paths[['path.res']], "/surrogates/", "patients_", patient, "_no_cross_auto_surrogate_optimum_test_", params[['sampling']], "_", params[['cov.type']], ".feather")
  if (file.exists(fname)){
        optimal.history.test.add = read_feather(fname)
        optimal.history.test = rbind(optimal.history.test, optimal.history.test.add)
  } else {
    print(paste('missing file : ', patient, cov.choice))
  }
}

write_feather(optimal.history.test, paste0(paths[['path.res']], "all_patients_no_cross_auto_surrogate_optimum_test_", params[['sampling']], "_", params[['cov.type']], ".feather"))
