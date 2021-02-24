#' @export


summary.oglmx<-function(object, ... ){
  stdEr.oglmx<-diag(object$vcov)^0.5
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

#' @export
print.summary.oglmx <- function(x, ... ){
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

#' @export

summary.oglmx.selection <- function(object, ... ){
  stdEr.oglmx<-diag(object$vcov)^0.5
  t<-object$coefAll/stdEr.oglmx
  p <- 2*pnorm( -abs( t))
  results <- cbind("Estimate"=object$coefAll,
                   "Std. error"=stdEr.oglmx,
                   "t value"=t, "Pr(>|t|)"=p)
  results_selection <- results[object$params$selection,]
  results_outcome <- results[object$params$outcome,]
  results_error <- results[object$params$error,]

  summary<-list(regtype=.regtype.oglmx(object),
                loglikelihood=object$loglikelihood,
                estimate=list("selection" = results_selection,
                              "outcome" = results_outcome,
                              "errors" = results_error),
                no.iterations=object$iterations,
                # McFaddensR2=McFaddensR2.oglmx(object),
                AIC=AIC(object),
                coefficients=object$estimate)
  class(summary) <- c("summary.oglmx.selection",
                      "summary.oglmx")
  summary
}


#' @export
print.summary.oglmx.selection <- function(x, ... ){
  cat(x$regtype,"\n")
  cat("Log-Likelihood:", x$loglikelihood, "\n")
  cat("Number of observations:", attr(x$loglikelihood, "No.Obs"))
  cat("No. Iterations:", x$no.iterations, "\n")
  cat("AIC:",x$AIC,"\n")

  cat("-----","Selection Model","------\n")
  printCoefmat(x$estimate$selection,signif.legend=FALSE, has.Pvalue = TRUE)
  cat("\n")
  cat("-----","Outcome Model","------\n")
  printCoefmat(x$estimate$outcome,signif.legend=FALSE, has.Pvalue = TRUE)
  cat("-----","Error terms parameters","------\n")
  printCoefmat(x$estimate$errors[c("log(sigma)","sigma","rho"),],signif.legend=TRUE, has.Pvalue = TRUE)
  cat("--------------------------------------------\n")
  invisible(x)
}
