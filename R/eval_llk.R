eval_llk <- function(Env,Parameters){
  updateComponents(Env,Parameters)
  return(loglikelihood.oglmx(Env))
}

eval_llk_point <- function(params, threshparams,
                           X, Z, outcomeMat, w,
                           link,
                           whichparametersmean,
                           whichparametersscale,
                           whichparametersthresh,
                           sdmodel,
                           analhessian){


  if (sum(whichparametersmean)>0){
    beta <- params[whichparametersmean]
    XB <- as.vector(X%*%beta)
  }

  if (sum(whichparametersscale)>0) {
    delta <- params[whichparametersscale]
    gsdmodel<-D(sdmodel,"z")
    hsdmodel<-D(gsdmodel,"z")
    ZD<-as.vector(Z%*%delta)
  } else{
    ZD <- 0
    gsdmodel<-NULL; hsdmodel<-NULL
  }
  Std.Dev<-eval({z<-ZD;sdmodel})
  GStd.Dev<-eval({z<-ZD;gsdmodel})
  minStD<-.Machine$double.xmin^(1/2)
  maxStD<-.Machine$double.xmax
  Std.Dev[Std.Dev==0]<-minStD
  Std.Dev[Std.Dev==Inf]<-maxStD
  GStd.Dev[GStd.Dev==Inf]<-maxStD

  if (analhessian){
    HStd.Dev<-eval({z<-ZD;hsdmodel})
    HStd.Dev[HStd.Dev==Inf]<-maxStD
  }

  if (sum(whichparametersthresh)>0){
    alphas <- params[whichparametersthresh]
    ThresholdMatrix <- getThresholds(outcomeMat, alphas)
  } else{
    ThresholdMatrix <- getThresholds(outcomeMat, threshparams)
  }

  etas<-getEtas(ThresholdMatrix,XB,Std.Dev)
  probs<-Probability(etas[[1]],etas[[2]],link = link)
  logprobs<-suppressWarnings(log(probs))

  if (!is.null(w)){
    wll<-sum(w*logprobs)
    return(wll)
  }

  ll<-sum(logprobs)
  return(ll)

}
