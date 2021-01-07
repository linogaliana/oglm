# collection of functions used to produce regression output


# working
nobs.oglmx<-function(object, ...){
  return(attr(object$loglikelihood,"No.Obs"))
}



# working
vcov.oglmx<-function(object,tol=1e-20,...){
  if (is.null(object$BHHHhessian)){
    vcov<-qr.solve(-object$hessian,tol=tol)
  } else {
    vcov<- qr.solve(object$hessian,tol=tol)%*%(object$BHHHhessian*(attr(object$loglikelihood,"No.Obs")/(attr(object$loglikelihood,"No.Obs")-1)))%*%qr.solve(object$hessian,tol=tol)
  }
  colnames(vcov)<-rownames(vcov)<-names(object$coefficients)
  return(vcov)
}

# working
summary.oglmx<-function(object,tol=1e-20, ... ){
  stdEr.oglmx<-diag(vcov(object,tol=tol))^0.5
  t<-object$coefficients/stdEr.oglmx
  p <- 2*pnorm( -abs( t))
  results <- cbind("Estimate"=object$coefficients,
                   "Std. error"=stdEr.oglmx,
                   "t value"=t, "Pr(>|t|)"=p)
  betaresults<-results[attr(object$coefficients,"coefftypes")[[1]], ,drop=FALSE]
  deltaresults<-results[attr(object$coefficients,"coefftypes")[[2]], ,drop=FALSE]
  cutoffresults<-results[attr(object$coefficients,"coefftypes")[[3]], ,drop=FALSE]
  resultsSplit<-list(betaresults,deltaresults,cutoffresults)
  summary<-list(regtype=.regtype.oglmx(object),loglikelihood=object$loglikelihood,estimate=results,estimateDisplay=resultsSplit,no.iterations=object$no.iterations,McFaddensR2=McFaddensR2.oglmx(object),AIC=AIC(object),coefficients=object$coefficients)
  class(summary)<-"summary.oglmx"
  summary
}

# working
print.summary.oglmx<-function(x, ... ){
  cat(x$regtype,"\n")
  cat("Log-Likelihood:", x$loglikelihood, "\n")
  cat("No. Iterations:", x$no.iterations, "\n")
  cat("McFadden's R2:",x$McFaddensR2,"\n")
  cat("AIC:",x$AIC,"\n")
  if (nrow(x$estimateDisplay[[1]])>0 & nrow(x$estimateDisplay[[2]])==0 & nrow(x$estimateDisplay[[3]])==0){
    printCoefmat(x$estimateDisplay[[1]])
  } else if (nrow(x$estimateDisplay[[1]])>0){
    if (nrow(x$estimateDisplay[[2]])>0){
      cat("-----","Mean Equation","------\n")
    }
    printCoefmat(x$estimateDisplay[[1]],signif.legend=FALSE)
  }
  if (nrow(x$estimateDisplay[[2]])>0){
    if (nrow(x$estimateDisplay[[1]])>0){
      cat("-----","SD Equation","------\n")
    }
    if (nrow(x$estimateDisplay[[3]])>0){
      printCoefmat(x$estimateDisplay[[2]],signif.legend=FALSE)
    } else {
      printCoefmat(x$estimateDisplay[[2]])
    }
  }
  if (nrow(x$estimateDisplay[[3]])>0){
    cat("-----","Threshold Parameters","-----\n")
    printCoefmat(x$estimateDisplay[[3]])
  }
}

# works
McFaddensR2.oglmx<-function(object){
  value<-1-logLik(object)/.BaseLL(object)
  return(value)
}



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

print.margins.oglmx<-function(x, ... ){
  for (m in 1:length(x)){
    cat("Marginal Effects on Pr(Outcome==",names(x)[m],")","\n",sep="")
    if (m==length(x)){
      printCoefmat(x[[m]])
    } else {
      printCoefmat(x[[m]],signif.legend=FALSE)
      cat("------------------------------------","\n")
    }
  }
}
