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
paths['path.res'] = paste0(rootResults, '/forecasting_over_days/', dataset, '/surrogates/')

# hyper parameter fitting on validation dataset
params = list()
patient.list.temp = read.csv(paste0(paths[['path.data']], 'patients.csv'), sep=",", header=F)
patient.list = as.character(patient.list.temp$V1)
params['cov.type'] = 'acausal4'
params['cv.option'] = 'block'
params['cov.choice'] = c('only_cross_auto')
params['sampling'] = "day" # "day", "hour"
params['cross.history.IEA.eval'] = 6
params['surrogate'] = 1
if (params[['cov.choice']]=='only_cross'){
	params['auto.history.Sz.max'] = 0
	params['cross.history.IEA.max'] = 5
} else if (params[['cov.choice']]=='only_auto'){
	params['auto.history.Sz.max'] = 5
	params['cross.history.IEA.max'] = 0
} else {
	params['auto.history.Sz.max'] = 5
	params['cross.history.IEA.max'] = 5
}

for (i.surr in c(1:200)){
  print('i.surr')
  print(i.surr)
  paths['path.process'] = paste0(rootProcess, '/forecasting_over_days/', dataset, '/surr_', toString(i.surr), '/') 
  paths['path.res'] = paste0(rootResults, '/forecasting_over_days/', dataset, '/surrogates', '/surr_', toString(i.surr), '/') 
  dir.create(paths[['path.res']])
  for (patient in patient.list){
    if (file.exists(paste0(paths[['path.res']], patient, "_predictions_", params[['cov.choice']], "_", params[['sampling']], "_", params[['cov.type']], ".nc"))){
      print('already processed')
    } else {
      pseHistory(paths, patient, params)
    }
  }
}