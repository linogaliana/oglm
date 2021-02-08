oglmx.fit<-function(outcomeMatrix,X,Z,w,beta,delta,threshparam,link,start,sdmodel,
                    optmeth = c("NR", "BFGS", "BFGSR", "BHHH", "SANN", "CG", "NM"),
                    analhessian,robust,
                    start_method = c("default","search"),
                    search_iter = 10){

  optmeth <- match.arg(optmeth)
  start_method <- match.arg(start_method)

  # outcomeMatrix<-results.oprobhetOGLMXENV[[1]]
  # X<-results.oprobhetOGLMXENV[[2]]
  # Z<-results.oprobhetOGLMXENV[[3]]
  # w<-results.oprobhetOGLMXENV[[4]]
  # beta<-results.oprobhetOGLMXENV[[5]]
  # delta<-results.oprobhetOGLMXENV[[6]]
  # threshparam<-results.oprobhetOGLMXENV[[7]]
  # start<-results.oprobhetOGLMXENV[[8]]
  # optmeth<-results.oprobhetOGLMXENV[[9]]
  # analhessian<-results.oprobhetOGLMXENV[[10]]
  # sdmodel<-results.oprobhetOGLMXENV[[11]]
  # robust<-results.oprobhetOGLMXENV[[12]]
  # link<-results.oprobhetOGLMXENV[[13]]


  whichXest<-is.na(beta)
  whichZest<-is.na(delta)
  no.est.beta<-sum(whichXest)
  no.est.delta<-sum(whichZest)
  whichAlphaest<-is.na(threshparam)
  no.est.thresh<-sum(whichAlphaest)
  no.estparams<-no.est.beta+no.est.delta+no.est.thresh
  if (!is.null(start) & length(start)!=no.estparams){
    stop("Specified vector of start values for parameters is of incorrect length.")
  }

  whichparametersmean<-c(!logical(no.est.beta),logical(no.est.delta+no.est.thresh))
  whichparametersthresh<-c(logical(no.est.beta+no.est.delta),!logical(no.est.thresh))
  whichparametersscale<-c(logical(no.est.beta),!logical(no.est.delta),logical(no.est.thresh))

  if (sum(whichZest)>0){
    gsdmodel<-D(sdmodel,"z")
    hsdmodel<-D(gsdmodel,"z")
  } else {gsdmodel<-NULL; hsdmodel<-NULL}

  CalcEnv<-list(outcomeMat=outcomeMatrix,X=X,Z=Z,w=w,beta=beta,delta=delta,threshparam=threshparam,
                link=link,sdmodel=sdmodel,gsdmodel=gsdmodel,hsdmodel=hsdmodel,whichparametersmean=whichparametersmean,
                whichparametersthresh=whichparametersthresh,whichparametersscale=whichparametersscale,
                whichXest=whichXest,whichZest=whichZest,whichAlphaest=whichAlphaest,analhessian=analhessian)
  CalcEnv<-list2env(CalcEnv)

  parametertypes<-list(whichparametersmean,whichparametersscale,whichparametersthresh)

  if (is.null(start)){
    start_algo<-calcstartvalues(parametertypes,sdmodel,threshparam)
  }

  if ((start_method != "default") && is.null(start)){
    search_loglik <- function(start_algo){
      # In that case, we replace 0 start by random values
      start_algo[start_algo == 0] <- rnorm(sum(start_algo == 0),
                                           sd = mean(abs(start_algo[start_algo != 0]), na.rm = TRUE)) #sd is arbitrary
      return(list("init" = start_algo,
                  # "llk" = eval_llk(CalcEnv, start_algo)
                  "llk" = eval_llk_point(start_algo, threshparams = threshparam,
                                         X = X, Z = Z, outcomeMat = outcomeMatrix, w = w,
                                         link = link,
                                         whichparametersmean = whichparametersmean,
                                         whichparametersscale = whichparametersscale,
                                         whichparametersthresh = whichparametersthresh,
                                         sdmodel = sdmodel,
                                         analhessian = analhessian)
      ))


    }
    draws <- lapply(seq_len(search_iter), function(i){
      out <- search_loglik(start_algo)
      return(
        list('out' = data.frame(iter = i, llk = out$llk),
             "input" = out$init))
    })
    draws_out <- do.call(rbind, lapply(draws, function(x) x$out))
    draws_out <- draws_out[is.finite(draws_out$llk),]
    idx_max <- draws_out$iter[draws_out$llk == max(draws_out$llk)][1]
    start_algo <- draws[[idx_max]]$input
    # updateComponents(CalcEnv,start_algo)
    # ll<-loglikelihood.oglmx(inputenv)
  }



  maximum<-oglmx.maxlik(CalcEnv,start_algo, optmeth = optmeth)
  results<-collectmaxLikOutput(maximum)

  outcomenames<-colnames(outcomeMatrix)
  threshnames<-sapply(c(2:length(outcomenames)),function(x){paste("Threshold (",outcomenames[x-1],"->",outcomenames[x],")",sep="")})
  if (ncol(Z) == 1) colnames(Z) <- "ln(sigma)"
  names(results$coefficients)<-c(colnames(X)[whichXest],colnames(Z)[whichZest],threshnames[whichAlphaest])


  attr(results$coefficients,"coefftypes")<-list(whichparametersmean,whichparametersscale,whichparametersthresh)
  results$allparams<-list(beta=CalcEnv$beta,delta=CalcEnv$delta,threshparam=CalcEnv$threshparam)
  results$Est.Parameters<-list(beta=whichXest,delta=whichZest,alpha=whichAlphaest)
  #return(list(CalcEnv,results))
  if (robust){
    results$BHHHmatrix<-calcBHHHmatrix(CalcEnv)
  } else {results$BHHHmatrix<-NULL}

  class(results)<-"oglmx.fit"
  #results<-list(loglikelihood,link=link,no.iterations=maxLikRes$iterations,coefficients=coefficients,returnCode=maxLikRes$code,gradient=maxLikRes$gradient
  #              ,hessian=maxLikRes$hessian,BHHHhessian=BHHHmatrix,NOutcomes=no.outcomes,Outcomes=listoutcomes,sdmodel=sdmodel,allparams=allparameters,Est.Parameters=Est.Parameters)
  return(results)

}

collectmaxLikOutput<-function(x){
  output<-list()
  output$loglikelihood<-x$maximum
  output$coefficients<-x$estimate
  output$gradient<-x$gradient
  output$no.iterations<-x$iterations
  output$returnCode<-x$code
  output$hessian<-x$hessian
  return(output)
}

oglmx.maxlik<-function(inputenv,start, optmeth = c("NR", "BFGS", "BFGSR", "BHHH", "SANN", "CG", "NM")){
  optmeth <- match.arg(optmeth)
  inputfunc<-function(par){
    updateComponents(inputenv,par)
    ll<-loglikelihood.oglmx(inputenv)
    score<-score_oglmx(inputenv)
    attr(ll,"gradient")<-score
    if (inputenv$analhessian){
      hessian<-hessian_oglmx(inputenv)
      attr(ll,"hessian")<-hessian
    }
    return(ll)
  }
  output<-maxLik::maxLik(inputfunc,start=start,iterlim=300,finalHessian=TRUE,method=optmeth) # ,control=list(printLevel=4)
  return(output)
}


# oglmx.maxlik2<-function(inputenv,start, optmeth = c("NR", "BFGS", "BFGSR", "BHHH", "SANN", "CG", "NM")){
#   optmeth <- match.arg(optmeth)
#
#   ll <- maxLik::maxLik(eval_llk_point,start=start,iterlim=300,finalHessian=TRUE,method=optmeth,
#          threshparams = threshparams,
#          X = X, Z = Z, outcomeMat = outcomeMatrix, w = w,
#          link = link,
#          whichparametersmean = whichparametersmean,
#          whichparametersscale = whichparametersscale,
#          whichparametersthresh = whichparametersthresh,
#          sdmodel = sdmodel,
#          analhessian = analhessian)
#
#   inputfunc<-function(par){
#     updateComponents(inputenv,par)
#     ll<-loglikelihood.oglmx(inputenv)
#     score<-score_oglmx(inputenv)
#     attr(ll,"gradient")<-score
#     if (inputenv$analhessian){
#       hessian<-hessian_oglmx(inputenv)
#       attr(ll,"hessian")<-hessian
#     }
#     return(ll)
#   }
#   output<-maxLik(inputfunc,start=start,iterlim=300,finalHessian=TRUE,method=optmeth) # ,control=list(printLevel=4)
#   return(output)
# }
