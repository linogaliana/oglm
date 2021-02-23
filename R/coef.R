#' @export
coef.oglmx<-function(object, ...){
  attr(object$coefficients,"coefftypes")<-NULL
  return(object$coefficients)
}

#' @export
coef.summary.oglmx<-function(object, ...){
  attr(object$coefficients,"coefftypes")<-NULL
  return(object$coefficients)
}

#' @export
coef.oglmx.selection<-function(object, ...){
  return(object$estimate)
}

#' @export
coef.summary.selection<-function(object, ...){
  return(object$estimate)
}

