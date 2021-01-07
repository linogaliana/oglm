# works
coef.oglmx<-function(object, ...){
  attr(object$coefficients,"coefftypes")<-NULL
  return(object$coefficients)
}

# works
coef.summary.oglmx<-function(object, ...){
  attr(object$coefficients,"coefftypes")<-NULL
  return(object$coefficients)
}
