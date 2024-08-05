# Author: Timothée Proix
# License: GPL-3.0-only

readData = function(paths, patient, params){
  raw.data = read.csv(paste0(paths[['path.process']], 'features/', 'features_', params[['sampling']], '_', patient, '_', params[['cov.type']], ".csv"))
  return(raw.data)
}


selectCovMulti = function(dat, params){
  # select the covariates and shift them in time to keep causality
  cov.choice = params[['cov.choice']]
  if (cov.choice=='phase_circ') {
    cov.multi.IEA = dat[ , grepl('CircadianCos', names(dat)) | grepl('CircadianSin', names(dat)), drop=FALSE]
  } else if (cov.choice=='phase_multi') {
    cov.multi.IEA = dat[ , grepl('MultidienCos', names(dat)) | grepl('MultidienSin', names(dat)), drop=FALSE]
  } else if (cov.choice=='phase_circ_multi') {
    cov.multi.IEA = dat[ , grepl('Cos', names(dat)) | grepl('Sin', names(dat)), drop=FALSE]
  } else if ((cov.choice=='only_cross_auto') || (cov.choice=='only_auto') || (cov.choice=='only_cross')) {
    cov.multi.IEA = NA
  } else if (cov.choice=='circ_hist') {
    cov.multi.IEA = dat[ , grepl('Daily', names(dat)), drop=FALSE]
  } else if (cov.choice=='week_hist') {
    cov.multi.IEA = dat[ , grepl('Weekly', names(dat)), drop=FALSE]
  } else if (cov.choice=='multivariate'){
    cov.multi.IEA = dat[ , grepl('Cos', names(dat)) | grepl('Sin', names(dat)) | grepl('Daily', names(dat)), drop=FALSE]
  } else {
    stop("cov.choice not valid")
  }
  return(as.matrix(cov.multi.IEA))
}


selectDat = function(dat, params){
  # To select the rest of the data
  t = dat['Time']
  dat.Sz = dat['Seizures']
  if (params[['surrogates']]==1 && params[['cov.choice']]=='only_auto'){
    dat.Sz = dat['SeizuresShuffled']
  }
  dat.IEA = dat[ , grepl('IEA', names(dat)), drop=FALSE]
    
  return(list(as.matrix(t), as.matrix(dat.Sz), as.matrix(dat.IEA)))
}