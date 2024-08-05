# Author: Timothée Proix
# License: GPL-3.0-only

library(pracma)

scoreAuc = function(labels, predictions, score_method){
  if (score_method=='prop_time_corr_1'){
    prop.time.auc = propTimeWarningAuc(labels, predictions, 1)
  }
}


scoreCurve = function(labels, predictions, score_method){
  if (score_method=='prop_time_corr_1'){
    list(y, x) := propTimeWarningCurve(labels, predictions, 1)
  }
  return(list(y, x))
}


propTimeWarning = function(labels, predictions, threshold, correction=0){
	# proportion of time in warning
  # (FP + TP)/(TP+FP+TN+FN)
  predictions[predictions<=threshold] = 0
  predictions[predictions>threshold] = 1

  if (correction==1){
    TP = sum((labels + predictions)==2)
    ptw = (sum(predictions) - TP)/(length(predictions) - TP)
  }
  return(ptw)
}


sensitivity = function(labels, predictions, threshold){
  # TPR, recall
  predictions[predictions<=threshold] = 0
  predictions[predictions>threshold] = 1
  TP = sum((labels+predictions)==2)
  FN = sum((labels-predictions)==1)
  TP/(TP+FN)
}


propTimeWarningAuc = function(labels, predictions, correction=0){
	# auc of sensitivity as function of proportion of time in warning
  list(sens, ptw) := propTimeWarningCurve(labels, predictions, correction, FALSE)
  tiw = trapz(ptw, sens)
  return(tiw)
}


propTimeWarningCurve = function(labels, predictions, correction=0, plot.show=FALSE){
  # auc of sensitivity as function of proportion of time in warning
  labels[labels>0]=1
  labels[is.na(labels)] = 0
  ptw = c(0)
  sens = c(0)
  nsteps = length(predictions)
  y = sort(predictions)
  for (istep in c(1:(nsteps))){
    ptw[istep+1] = propTimeWarning(labels, predictions, y[istep], correction)
    sens[istep+1] = sensitivity(labels, predictions, y[istep])
  }
  if (correction==1){
    ptw[nsteps+2] = 1
    sens[nsteps+2] = 1
  }
  idx = order(sens, ptw)
  sens = sens[idx]
  ptw = ptw[idx]
  if (plot.show){
    X11()
    plot(ptw, sens, type='l')
  }
  return(list(sens, ptw))
}


findMax = function(m, params){
  m[1, 1] = 0 # Constant prediction is 0.5, set 0 to not choose this one
  all.inds = which(m == max(m, na.rm=TRUE), arr.ind=TRUE)
  if (params[['cov.choice']]=='only_auto'){
    all.inds = which(m[,1,drop=FALSE] == max(m[,1,drop=FALSE], na.rm=TRUE), arr.ind=TRUE)
  }
  if (params[['cov.choice']]=='only_cross'){
    all.inds = which(m[1,,drop=FALSE] == max(m[1,,drop=FALSE], na.rm=TRUE), arr.ind=TRUE)
  }
  inds = all.inds[dim(all.inds)[1],]
  auto.history.Sz.max = inds[1]
  cross.history.IEA.max = inds[2]
  # print(auto.history.Sz.max)
  # print(cross.history.IEA.max)
  return(list(m[auto.history.Sz.max, cross.history.IEA.max], auto.history.Sz.max, cross.history.IEA.max))
}