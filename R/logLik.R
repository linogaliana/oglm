logLik.oglmx<-function(object, ...){
  value<-object$loglikelihood[1]
  attr(value,"df")<-length(object$coefficients)
  return(value)
}

# works
logLik.summary.oglmx<-function(object, ...){
  value<-object$loglikelihood[1]
  attr(value,"df")<-length(object$coefficients)
  return(value)
}

# works with output of new version, need to check with output of old version
.BaseLL<-function(object){
  data<-object$NoVarModData
  if (ncol(data)==2){
    BaseLL<-as.numeric(logLik(oglmx(Y~1,data=data,weights=data$weights)))
  } else {
    BaseLL<-as.numeric(logLik(oglmx(Y~1,data=data)))
  }
  return(BaseLL)
}
