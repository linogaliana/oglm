#' Calculate Variance-Covariance Matrix for a Fitted Model Object
#' Returns the variance-covariance matrix of the main parameters of a
#'   fitted model object of class `oglm`.See \link[stats]{vcov} for
#'   more details
#'
#' @inheritParams stats::vcov
#' @inheritParams qr.solve
#' @inheritParams oglmx
#' @export

vcov.oglmx<-function(object,tol=1e-20,...){
  if (is.null(object$BHHHhessian)){
    vcov<-qr.solve(-object$hessian,tol=tol)
  } else {
    vcov<- qr.solve(object$hessian,tol=tol)%*%(object$BHHHhessian*(attr(object$loglikelihood,"No.Obs")/(attr(object$loglikelihood,"No.Obs")-1)))%*%qr.solve(object$hessian,tol=tol)
  }
  colnames(vcov)<-rownames(vcov)<-names(object$coefficients)
  return(vcov)
}

#' Calculate Variance-Covariance Matrix for a Fitted Model Object
#' Returns the variance-covariance matrix of the main parameters of a
#'   fitted model object of class `oglm`.See \link[stats]{vcov} for
#'   more details
#' @inheritParams stats::vcov
#' @export

vcov.oglmx.selection <-function(object, ...){
  return(object$vcov)
}
