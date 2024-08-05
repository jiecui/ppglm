# Author: Timothée Proix
# License: GPL-3.0-only

library(tscount)
library(ncdf4)


trainSzCv = function(train.Sz, cov.train, cov.validation, length.validation){
  # fit on train daFautota
  Sz.with.IEA.fit = tsglm(train.Sz, xreg=cov.train, link="log", distr="poisson")
  # predict on validation data 
  # issue when data are too short, trycatch them
  preds.validation = tryCatch({
    predict(Sz.with.IEA.fit, n.ahead=length.validation, 
                             newxreg=cov.validation, type="shortest")
    }, error = function(error_message){
          return(NA)
  })
  return(preds.validation)
}


#' Cross-validation
predictionValidationCv = function(paths, patient, params){

  auto.history.Sz = params[['auto.history.Sz']]
  cross.history.IEA = params[['cross.history.IEA']]
  cross.history.IEA.eval = params[['cross.history.IEA.eval']]
  if (cross.history.IEA.eval<=cross.history.IEA){
      stop("cross.history.IEA.eval should be at least (cross.history.IEA + 1)")
  }
  print(auto.history.Sz)
  print(cross.history.IEA)
  print(patient)

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

  # set the auto_history lists
  if (auto.history.Sz>0){
      auto.past.obs = c(1:auto.history.Sz)
  } else {
      auto.past.obs = c(1)
  }

  # for the IEA cross-history, build a matrix
  if (cross.history.IEA>0){
      cross.history.IEA.list = c(1:cross.history.IEA)
  } else { # just to avoid a special case for the return values of the function
      cross.history.IEA.list = c(1)
  }

  # build the Sz auto-history covariate
  cov.Sz = matrix(dat.Sz[1:(length(dat.Sz)-1)])
  if (length(auto.past.obs)>1){
    for (ilag in c(2:length(auto.past.obs))){
        cov.Sz = cbind(cov.Sz[(1+auto.past.obs[ilag]-auto.past.obs[ilag-1]):dim(cov.Sz)[1],],
                       dat.Sz[1:(length(dat.Sz)-auto.past.obs[ilag])])
    }
  } else {
      ilag=1
  }

  # build the IEA cross-history covariate
  cov.IEA = dat.IEA[1:(dim(dat.IEA)[1]-1),,drop=FALSE]
  if (length(cross.history.IEA.list)>1){
      for (ilag in c(2:length(cross.history.IEA.list))){
          cov.IEA = cbind(cov.IEA[(1+cross.history.IEA.list[ilag]-cross.history.IEA.list[ilag-1]):dim(cov.IEA)[1],],
                          dat.IEA[1:(dim(dat.IEA)[1]-cross.history.IEA.list[ilag]),])
          }
  } else {
      ilag=1
  }
  cov.Sz.cut = cov.Sz[(cross.history.IEA.eval-length(auto.past.obs)):(length(dat.Sz)-length(auto.past.obs)),,drop=FALSE]
  cov.IEA.cut = cov.IEA[(cross.history.IEA.eval-length(cross.history.IEA.list)):(length(dat.Sz)-length(cross.history.IEA.list)),,drop=FALSE]
  dat.Sz.cut = dat.Sz[cross.history.IEA.eval:length(dat.Sz)]

  # length of train, validation, and test data
  # we remove cross.history.eval points in length.train to ease future comparison
  if (params[['sampling']]=="day"){
      length.both = min(floor(0.60*length(dat.Sz)), 48*8) - (cross.history.IEA.eval - 1)
      length.train = floor(0.80*length.both)
      length.validation = length.both - length.train
  } else { 
      length.both = min(floor(0.60*length(dat.Sz)), 48*8*24) - (cross.history.IEA.eval - 1)
      length.train = floor(0.80*length.both)
      length.validation = length.both - length.train
  }

  # cross-validation loop
  m.validation.Sz = matrix(, nrow=5, ncol=length.validation) 
  m.preds.validation = matrix(, nrow=5, ncol=length.validation)

  for (i.cv in c(1:5)){

    # pairewise shuffling of covariate and Sz auto-history
    if (params[['cv.option']]=='block') {
      #pos.block.val = sample(c(1:(length.train-1)), size=(1))
      pos.block.val = (i.cv-1)*length.validation

      if (i.cv<5){
        pos = c(seq_len(pos.block.val), 
                pos.block.val + length.validation + seq_len(length.train - pos.block.val),
                c((pos.block.val+1):(pos.block.val + length.validation)))
      } else if (i.cv==5){
        pos = c(1:length.both)        
      }
    }

    cov.Sz.cut.cv = cov.Sz.cut[1:(length.train+length.validation),, 
                                 drop=FALSE][pos,,drop=FALSE] 
    cov.IEA.cut.cv = cov.IEA.cut[1:(length.train+length.validation),,
                                   drop=FALSE][pos,,drop=FALSE]
    dat.Sz.cut.cv = dat.Sz.cut[1:(length.train+length.validation)][pos]
     
    cov.train.Sz = cov.Sz.cut.cv[1:length.train,, drop=FALSE]
    cov.train.IEA = cov.IEA.cut.cv[1:length.train,, drop=FALSE]
    train.Sz = dat.Sz.cut.cv[1:length.train]
     
    cov.validation.Sz = cov.Sz.cut.cv[(length.train+1):(length.train+length.validation),, drop=FALSE]
    cov.validation.IEA = cov.IEA.cut.cv[(length.train+1):(length.train+length.validation),, drop=FALSE]
    validation.Sz = dat.Sz.cut.cv[(length.train+1):(length.train+length.validation)]
    m.validation.Sz[i.cv, ] = validation.Sz
      
    if (!(is.na(cov.multi.IEA[1]))){
      cov.multi.IEA.cut = cov.multi.IEA[params[['cross.history.IEA.eval']]:length(dat.Sz),,drop=FALSE]
      # pairewise shuffling of mulidian covariate
      cov.multi.IEA.cut.cv = cov.multi.IEA.cut[1:(length.train+length.validation),,drop=FALSE][pos,,drop=FALSE]
      cov.multi.train.IEA = cov.multi.IEA.cut.cv[1:length.train,,drop=FALSE]
      cov.train.IEA = cbind(cov.train.IEA, cov.multi.train.IEA)
      cov.multi.validation.IEA = cov.multi.IEA.cut.cv[(length.train+1):(length.train+length.validation),,drop=FALSE]
      cov.validation.IEA = cbind(cov.validation.IEA, cov.multi.validation.IEA)
    }

    if (auto.history.Sz==0 && cross.history.IEA==0){
      if (!(is.na(cov.multi.IEA[1]))){
        preds.validation = trainSzCv(train.Sz, cov.multi.train.IEA, cov.multi.validation.IEA, length.validation)
      } else {
        # no actual prediction here
        preds.validation = list()
        preds.validation['pred'] = numeric(length(validation.Sz))
      }

    } else if (auto.history.Sz==0 && cross.history.IEA>0){

      if (!(is.na(cov.multi.IEA[1]))){
        cov.train = cbind(cov.train.IEA, cov.multi.train.IEA)
        cov.validation = cbind(cov.validation.IEA, cov.multi.validation.IEA)
        preds.validation = trainSzCv(train.Sz, cov.train, cov.validation, length.validation)
      } else {
        preds.validation = trainSzCv(train.Sz, cov.train.IEA, cov.validation.IEA, length.validation)
      }
         
    } else if (auto.history.Sz>0 && cross.history.IEA==0){

      if (!(is.na(cov.multi.IEA[1]))){
        cov.train = cbind(cov.train.Sz, cov.multi.train.IEA)
        cov.validation = cbind(cov.validation.Sz, cov.multi.validation.IEA)
        preds.validation = trainSzCv(train.Sz, cov.train, cov.validation, length.validation)
      } else {
        preds.validation = trainSzCv(train.Sz, cov.train.Sz, cov.validation.Sz, length.validation)
                             }
    } else if (auto.history.Sz>0 && cross.history.IEA>0){
      if (!(is.na(cov.multi.IEA[1]))){
        cov.train = cbind(cov.train.Sz, cov.train.IEA, cov.multi.train.IEA)
        cov.validation = cbind(cov.validation.Sz, cov.validation.IEA, cov.multi.validation.IEA)
        preds.validation = trainSzCv(train.Sz, cov.train, cov.validation, length.validation)
      } else {
        cov.train = cbind(cov.train.Sz, cov.train.IEA)
        cov.validation = cbind(cov.validation.Sz, cov.validation.IEA)
        preds.validation = trainSzCv(train.Sz, cov.train, cov.validation, length.validation)
      }
    }

    if (is.na(preds.validation)){
        print('bad cv')
    } else { 
      m.preds.validation[i.cv,] = preds.validation$pred
    }

  }

  return(list(m.preds.validation, m.validation.Sz))
}

pseHistory = function(paths, patient, params){

  for (auto.history.Sz in c(0:params[['auto.history.Sz.max']])){
    for (cross.history.IEA in c(0:params[['cross.history.IEA.max']])){
      params['auto.history.Sz'] = auto.history.Sz
      params['cross.history.IEA'] = cross.history.IEA
      list(m.preds.validation, m.validation.Sz) := predictionValidationCv(paths, patient, params)
      if (!exists("predictions.labels")){
        predictions.labels = array(data=NA, dim=c(2, params[['auto.history.Sz.max']]+1, params[['cross.history.IEA.max']]+1, dim(m.preds.validation)[1], dim(m.validation.Sz)[2]))
      }
      predictions.labels[1, auto.history.Sz+1, cross.history.IEA+1,,] = m.preds.validation
      predictions.labels[2, auto.history.Sz+1,cross.history.IEA+1,,] = m.validation.Sz
    }
  }

  if (!exists('no_save', params) || params[["no_save"]]==0){
    # save in hdf5/ncdf4 format
    nc.predictions.labels.dim = ncdim_def("labels", "na", c(1:2))
    nc.auto.history.Sz.dim = ncdim_def("autohistory", "na", c(1:dim(predictions.labels)[2]))
    nc.cross.history.Sz.dim = ncdim_def("crosshistory", "na", c(1:dim(predictions.labels)[3]))
    nc.fold.dim = ncdim_def("fold", "na", c(1:dim(predictions.labels)[4]))
    nc.ts.dim = ncdim_def("ts", "na", c(1:dim(predictions.labels)[5]))
    nc.predictions.labels = ncvar_def("predictions.labels", "na", list(nc.predictions.labels.dim, nc.auto.history.Sz.dim, nc.cross.history.Sz.dim, nc.fold.dim,  nc.ts.dim))
    ncfname = paste0(paths[['path.res']], patient, "_predictions_", params[['cov.choice']], "_", params[['sampling']], "_", params[['cov.type']], ".nc")
    print(paste0(paths[['path.res']], patient, "_predictions_", params[['cov.choice']], "_", params[['sampling']], "_", params[['cov.type']], ".nc"))
    ncout = nc_create(ncfname,list(nc.predictions.labels),force_v4=TRUE)
    ncvar_put(ncout, nc.predictions.labels, predictions.labels)
    nc_close(ncout)
  }

  return(predictions.labels)
}

