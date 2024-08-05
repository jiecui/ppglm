# Author: Timothée Proix
# License: GPL-3.0-only

library(tscount)

predictionTest = function(paths, patient, params){
 
  auto.history.Sz = params[['auto.history.Sz']]
  cross.history.IEA = params[['cross.history.IEA']]
  cross.history.IEA.eval = params[['cross.history.IEA.eval']]
  if (cross.history.IEA.eval<=cross.history.IEA){
      stop("cross.history.IEA.eval should be at least (cross.history.IEA + 1)")
  }
  

  # get features
  raw.data = readData(paths, patient, params)
  cov.multi.IEA = selectCovMulti(raw.data, params)
  list(t, dat.Sz, dat.IEA) := selectDat(raw.data, params)

  # binarize seizure times series
  threshold = 1  # for patients with low number of seizures
  category.dat.Sz = dat.Sz
  category.dat.Sz[category.dat.Sz<threshold] = 0
  category.dat.Sz[category.dat.Sz>0] = 1
  dat.Sz = category.dat.Sz

  ## Construct and shift covariate by one time point to ensure causality
  # for cov.multi.IEA we shift forward and repeat first time point
  cov.multi.IEA = rbind(cov.multi.IEA[1,,drop=FALSE], cov.multi.IEA[1:(nrow(cov.multi.IEA) - 1),,drop=FALSE])

  # for the IEA cross-history, build a matrix
  if (cross.history.IEA>0){
      cross.history.IEA.list = c(1:cross.history.IEA)
  } else { # just to avoid a special case for the return values of the function
      cross.history.IEA.list = c(1)
  }
  
  cov.IEA = dat.IEA[1:(dim(dat.IEA)[1]-1),,drop=FALSE]
  if (length(cross.history.IEA.list)>1){
      for (ilag in c(2:length(cross.history.IEA.list))){
          cov.IEA = cbind(cov.IEA[(1+cross.history.IEA.list[ilag]-cross.history.IEA.list[ilag-1]):dim(cov.IEA)[1],],
                          dat.IEA[1:(dim(dat.IEA)[1]-cross.history.IEA.list[ilag]),])
          }
  } else {
      ilag=1
  }
  cov.IEA.cut = cov.IEA[(cross.history.IEA.eval-length(cross.history.IEA.list)):(length(dat.Sz)-length(cross.history.IEA.list)),,drop=FALSE]
  dat.Sz.cut = dat.Sz[cross.history.IEA.eval:length(dat.Sz)]
  
  # length of train, validation, and test data
  # we remove cross.history.eval points in length.train to ease future comparison
  if (params[['sampling']]=="day"){
      length.both = min(floor(0.60*length(dat.Sz)), 48*8) - (cross.history.IEA.eval + 1)
      length.train = floor(0.80*length.both) 
      length.validation = length.both - length.train
  } else { 
      length.both = min(floor(0.60*length(dat.Sz)), 48*8*24) - (cross.history.IEA.eval + 1)
      length.train = floor(0.80*length.both)
      length.validation = length.both - length.train
  }

  # if we want a specific training length
  if (length(params[['length.train']])!=0) {
    length.train = params[['length.train']]
    length.validation = 0 # we do not use validation for the test
  }

  # if we don't want to use the full training set
  cut.training = params[['cut.training']]
  if (length(cut.training)==0) {cut.training = length.train + length.validation - 1}

  cov.train.IEA = cov.IEA.cut[(length.train + length.validation - cut.training):(length.train+length.validation),, drop=FALSE]
  train.Sz = dat.Sz.cut[(length.train + length.validation - cut.training):(length.train+length.validation)]

  # if we want a specific testing length
  length.test = length.train + length.validation + params[['length.test']]
  if (length(params[['length.test']])==0){length.test = length(dat.Sz.cut)}

  cov.test.IEA = cov.IEA.cut[(length.train+length.validation+1):length.test,, drop=FALSE]
  test.Sz = dat.Sz.cut[(length.train+length.validation+1):length.test]
  if (!(is.na(cov.multi.IEA[1]))){
      cov.multi.IEA.cut = cov.multi.IEA[cross.history.IEA.eval:length(dat.Sz),,drop=FALSE]
      cov.multi.train.IEA = cov.multi.IEA.cut[(length.train + length.validation - cut.training):(length.train+length.validation),,drop=FALSE]
      cov.train.IEA = cbind(cov.train.IEA, cov.multi.train.IEA)
      cov.multi.test.IEA = cov.multi.IEA.cut[(length.train+length.validation+1):length.test,,drop=FALSE]
      cov.test.IEA = cbind(cov.test.IEA, cov.multi.test.IEA)
  }

  ## prediction
  # if we want predictions more than one step ahead
  pred.n.ahead = length((length.train+length.validation+1):length.test) - 1 + params[['n.ahead']]
  if (length(params[['n.ahead']])==0){pred.n.ahead = length((length.train+length.validation+1):length.test)}

  if (auto.history.Sz==0 && cross.history.IEA==0){
    if (!(is.na(cov.multi.IEA[1]))){
        # fit on train data
        Sz.with.IEA.fit = tsglm(train.Sz, xreg=cov.multi.train.IEA, link="log", distr="poisson")

        # evaluate auc score on validation data
        preds.test = predict(Sz.with.IEA.fit, n.ahead=pred.n.ahead, 
                                   newobs=test.Sz, newxreg=cov.multi.test.IEA, 
                                   method="conddistr", level=0)
    } else {
      stop("cross.history.IEA and auto.history.Sz are 0 without any covariate")
    }
  } else if (auto.history.Sz==0 && cross.history.IEA>0){
      # fit on train data
      Sz.with.IEA.fit = tsglm(train.Sz, xreg=cov.train.IEA, link="log", distr="poisson")

      # evaluate auc score on test data
      preds.test = predict(Sz.with.IEA.fit, n.ahead=pred.n.ahead, newxreg=cov.test.IEA, level=0)

  } else if (auto.history.Sz>0 && cross.history.IEA==0){

      if (!(is.na(cov.multi.IEA[1]))){
          # fit on train data
          Sz.with.IEA.fit = tsglm(train.Sz, model=list(past_obs=c(1:auto.history.Sz)),
                                  xreg=cov.multi.train.IEA, link="log", distr="poisson")

         # evaluate auc score on validation data
          preds.test = predict(Sz.with.IEA.fit, n.ahead=pred.n.ahead, 
                                     newobs=test.Sz, newxreg=cov.multi.test.IEA, 
                                     method="conddistr", level=0)
      } else {
          # fit on train data
          Sz.with.IEA.fit = tsglm(train.Sz, model=list(past_obs=c(1:auto.history.Sz)),
                                  link="log", distr="poisson")

          # evaluate auc score on validation data
          preds.test = predict(Sz.with.IEA.fit, n.ahead=pred.n.ahead, 
                                     newobs=test.Sz, method="conddistr", level=0)
      }

  } else if (auto.history.Sz>0 && cross.history.IEA>0){
      Sz.with.IEA.fit = tsglm(train.Sz, model=list(past_obs=c(1:auto.history.Sz)), xreg=cov.train.IEA, 
                                                         link="log", distr="poisson")
      
      preds.test = predict(Sz.with.IEA.fit, n.ahead=pred.n.ahead, newobs=test.Sz, newxreg=cov.test.IEA, level=0)

  }

  return(list(preds.test, test.Sz, cov.test.IEA, Sz.with.IEA.fit))
}