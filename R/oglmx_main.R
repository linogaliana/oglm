#' Fit Ordered Generalized Linear Model.
#'
#' @description oglmx is used to estimate models for which
#' the outcome variable is discrete and the mean and/or
#' variance of the underlying latent variable can be
#' modelled as a linear combination of explanatory variables.
#' Standard models such as probit, logit, ordered probit
#' and ordered logit are included in the diverse
#' set of models estimated by the function.
#'
#' @param formulaMEAN an object of class \link[stats]{formula}:
#'  a symbolic description of the model used to explain the
#'  mean of the latent variable. The response variable
#'  should be a numeric vector or factor variable
#'  such that the numerical assignments for the
#'  levels of the factor have ordinal meaning.
#' @param formulaSD either `NULL` (homoskedastic model) or an object
#'   of class \link[stats]{formula}: a symbolic description of
#'   the model used to explain the variance of the latent variable.
#' @param selection Formula for Heckman selection model. If `NULL` (default),
#'   assuming no selection on observables
#'   (introduced in [#11](https://github.com/linogaliana/oglm/pull/11))
#' @param data a data frame containing the variables in the model.
#' @param start either `NULL` or a numeric vector specifying
#'   start values for each of the estimated parameters,
#'   passed to the maximisation routine.
#' @param link specifies a link function for the model to be
#'   estimated, accepted values are *"probit"*,
#'   *"logit"*, *"cauchit"*, *"loglog"* and *"cloglog"*
#' @param constantMEAN logical. Should an intercept be
#'   included in the model of the mean of the latent variable?
#'   Can be overwritten and set to `FALSE` using the
#'   `formulaMEAN` argument by writing `0 +`
#'   as the first element of the equation.
#' @param constantSD logical. Should an intercept be included
#'   in the model of the variance of the latent variable?
#'   Can be overwritten and set to `FALSE` using the `formulaSD`
#'   argument by writing `0 +` as the first element of the equation.
#' @param beta `NULL` or numeric vector. Used to prespecify
#'   elements of the parameter vector for the equation of the mean
#'   of the latent variable. Vector should be of length one or of
#'   length equal to the number of explanatory
#'   variables in the mean equation.
#'   If of length one the value is presumed to correspond to the
#'   constant if a constant is included or the first element of the
#'   parameter vector. If of length greater than one then `NA`
#'   should be entered for elements of the vector to be estimated.
#' @param delta `NULL` or numeric vector. Used to prespecify
#'   elements of the parameter vector for the equation of
#'   the variance of the latent variable. Vector should be of
#'   length one or of length equal to the number of explanatory
#'   variables in the variance equation.
#'   If of length one the value is presumed to
#'   correspond to the constant if a constant is included
#'   or the first element of the parameter vector. If
#'   of length greater than one then `NA` should be entered
#'   for elements of the vector to be estimated.
#' @param threshparam `NULL` or numeric vector. Used to prespecify
#'   the threshold parameters of the model. Vector should be
#'   of length equal to the number of outcomes minus one. `NA` should
#'   be entered for threshold parameters to be estimated by the model.
#' @param analhessian logical. Indicates whether the analytic
#'   Hessian should be calculated and used, default is
#'   `TRUE`, if set to `FALSE` a finite-difference
#'   approximation of the Hessian is used.
#' @param sdmodel object of mode *“expression”*. The expression
#'   defines the function that transforms the linear model
#'   for the standard deviation into the standard deviation.
#'   The expression should be written as a function of variable `z`.
#'   The default value is `expression(exp(z))`.
#' @param SameModelMEANSD logical. Indicates whether the matrix used to
#'   model the mean of the latent variable is identical to that used
#'   to model the variance. If `formulaSD=NULL` and `SameModelMEANSD=TRUE` a
#'   model with heteroskedasticity is estimated. If
#' `SameModelMEANSD=FALSE` and `formulaSD==formulaMEAN`
#'   value is overridden. Used to reduce memory requirements when models are identical.
#' @param savemodelframe logical. Indicates whether the
#'   model frame(s) should be saved for future use.
#'   Default is `FALSE`. Should be set to `TRUE` if
#'   intending to estimate Average Marginal Effects.
#' @param Force logical. If set to `FALSE` (the default),
#'   the function stops if the response
#'   variable has more than twenty categories.
#'   Should be changed to `TRUE` if a model with
#'   more than twenty categories is desired.
#' @param robust logical. If set to `TRUE` the outer product or
#'   BHHH estimate of the meat in the sandwich of the
#'   variance-covariance matrix is calculated. If calculated
#'   standard errors will be calculated using the sandwich estimator
#'   by default when calling `summary`.
#' @param optmeth specifies a method for the maximisation of the likelihood
#'   passed to [maxLik::maxLik()]. Default to *NR* (Newton-Raphson) when
#'   no selection is introduced. Forced to *BHHH* when selection on observables
#'   is introduced.
#' @param gradient Should we use analytical gradient (default)
#'   or numerical gradient ?
#'   Analytical gradient results in slower iterations but sometimes help
#'   the model to converge faster.
#' @param return_envir Logical indicating whether we want to stop early and
#'   return objects used to fit the model
#' @param tol Argument passed to [qr.solve], defines the tolerance
#' for detecting linear dependencies in the hessian matrix to be inverted.
#' @inheritParams stats::glm
#'
#' @return An object of class "\code{oglmx}" with the following components:
#' \describe{
#'     \item{loglikelihood}{log-likelihood for the estimated model.
#'       Includes as attributes the log-likelihood for the constant
#'       only model and the number of observations.}
#'     \item{link}{link function used in the estimated model.}
#'     \item{no.iterations}{number of iterations of maximisation algorithm.}
#'     \item{coefficients}{named vector of estimated parameters.}
#'     \item{returnCode}{code returned by
#'       the `maxLik` optimisation routine}
#'     \item{call}{the call used to generate the results.}
#'     \item{gradient}{numeric vector, the value of the
#'       gradient of the log-likelihood function at the obtained
#'       parameter vector. Should be approximately equal to zero.}
#'     \item{terms}{two element list. Each element is an object of
#'       type `terms` related to the mean and standard
#'       deviation equation respectively.
#'     }
#'     \item{formula}{two element list. Each element is an
#'       object of type [stats::formula()] related to
#'       the mean and standard deviation
#'       equation respectively.}
#'     \item{NoVarModData}{dataframe. Contains data required to
#'       estimate the no information model used in
#'       calculation of McFadden's R-squared measure.}
#'     \item{hessian}{hessian matrix of the log-likelihood
#'       function evaluated at the obtained parameter vector.}
#'     \item{BHHHhessian}{Either `NULL` if no weights were
#'       included and `robust = FALSE`, or the BHHH estimate.}
#'     \item{Hetero}{logical. If `TRUE` indicates that the
#'       estimated model includes a model for the variance
#'       of the error term, i.e. heteroskedasticity.}
#'     \item{NOutcomes}{the number of distinct outcomes
#'       in the response variable.}
#'     \item{Outcomes}{numeric vector of length equal to `NOutcomes`.
#'       Lists the values of the different outcomes.}
#'     \item{BothEq}{data.frame with either two or three
#'       columns. Lists the names of variables that are in both
#'       the mean and variance equations and their locations
#'       within their respective model frames. Information is required in
#'       the call of `margins.oglmx` to obtain correct marginal effects.}
#'     \item{allparams}{a list containing three numeric vectors,
#'       the vectors contain the parameters from the mean equation,
#'       the variance equation and the threshold parameters respectively.
#'       Includes the prespecified and estimated parameters together.}
#'     \item{varMeans}{a list containing two numeric vectors. The vectors
#'       list the mean values of the variables in the mean and variance
#'       equation respectively. Stored for use in a call of
#'       \code{margins.oglmx} to obtain marginal effects at means.}
#'     \item{varBinary}{a list containing two numeric vectors. The vectors
#'       indicate whether the variables in the mean and variance
#'       equations are binary indicators. Stored for use in a call
#'       of \code{margins.oglmx} to obtain marginal effects at means.}
#'     \item{Est.Parameters}{list containing three logical vectors.
#'       Indicates which parameters in the parameter vectors were estimated.}
#'     \item{modelframes}{If \code{savemodelframe} set to \code{FALSE}
#'       then returns \code{NULL}, otherwise returns a list with two
#'       elements, the model frames for the mean and variance equations.}
#' }
#'
#' @examples \dontrun{
#' # create random sample, three variables, two binary.
#' set.seed(242)
#' n<-250
#' x1<-sample(c(0,1),n,replace=TRUE,prob=c(0.75,0.25))
#' x2<-vector("numeric",n)
#' x2[x1==0]<-sample(c(0,1),n-sum(x1==1),replace=TRUE,prob=c(2/3,1/3))
#' z<-rnorm(n,0.5)
#' # create latent outcome variable
#' latenty<-0.5+1.5*x1-0.5*x2+0.5*z+rnorm(n,sd=exp(0.5*x1-0.5*x2))
#' # observed y has four possible values: -1,0,1,2
#' # threshold values are: -0.5, 0.5, 1.5.
#' y<-vector("numeric",n)
#' y[latenty< -0.5]<- -1
#' y[latenty>= -0.5 & latenty<0.5]<- 0
#' y[latenty>= 0.5 & latenty<1.5]<- 1
#' y[latenty>= 1.5]<- 2
#' dataset<-data.frame(y,x1,x2)
#' # estimate standard ordered probit
#' results.oprob<-oglmx(y ~ x1 + x2 + z, data=dataset,link="probit",constantMEAN=FALSE,
#'                      constantSD=FALSE,delta=0,threshparam=NULL)
#' coef(results.oprob) # extract estimated coefficients
#' summary(results.oprob)
#' # calculate marginal effects at means
#' margins.oglmx(results.oprob)
#' # estimate ordered probit with heteroskedasticity
#' results.oprobhet<-oglmx(y ~ x1 + x2 + z, ~ x1 + x2, data=dataset, link="probit",
#'                         constantMEAN=FALSE, constantSD=FALSE,threshparam=NULL)
#' summary(results.oprobhet)
#' library("lmtest")
#' # likelihood ratio test to compare model with and without heteroskedasticity.
#' lrtest(results.oprob,results.oprobhet)
#' # calculate marginal effects at means.
#' margins.oglmx(results.oprobhet)
#' # scale of parameter values is meaningless. Suppose instead two of the
#' # three threshold values were known, then can include constants in the
#' # mean and standard deviation equation and the scale is meaningful.
#' results.oprobhet1<-oglmx(y ~ x1 + x2 + z, ~ x1 + x2, data=dataset, link="probit",
#'                          constantMEAN=TRUE, constantSD=TRUE,threshparam=c(-0.5,0.5,NA))
#' summary(results.oprobhet1)
#' margins.oglmx(results.oprobhet1)
#' # marginal effects are identical to results.oprobithet, but using the true thresholds
#' # means the estimated parameters are on the same scale as underlying data.
#' # can choose any two of the threshold values and get broadly the same result.
#' results.oprobhet2<-oglmx(y ~ x1 + x2 + z, ~ x1 + x2, data=dataset, link="probit",
#'                          constantMEAN=TRUE, constantSD=TRUE,threshparam=c(-0.5,NA,1.5))
#' summary(results.oprobhet2)
#' margins.oglmx(results.oprobhet2)
#' # marginal effects are again identical. Parameter estimates do change.
#' }
#' @param start_method Should we use default intiial value or search for a better one ?
#' @param search_iter Number of values to look for when using *start_method = 'search'*
#'
#' @import maxLik
#' @import stats
#' @export
oglmx<-function(formulaMEAN, formulaSD=NULL,
                selection = NULL,
                data, start=NULL, weights=NULL, link="probit",
                constantMEAN=TRUE, constantSD=TRUE, beta=NULL, delta=NULL, threshparam=NULL,
                analhessian=TRUE, sdmodel=expression(exp(z)), SameModelMEANSD=FALSE, na.action,
                savemodelframe=TRUE, Force=FALSE, robust=FALSE,
                optmeth = c("NR", "BFGS", "BFGSR", "BHHH", "SANN", "CG", "NM"),
                gradient = c("analytical","numerical"),
                tol=1e-20,
                start_method = c("default","search"),
                search_iter = 10,
                return_envir = FALSE){

  optmeth <- match.arg(optmeth)
  start_method <- match.arg(start_method)
  gradient <- match.arg(gradient)


  cl<-match.call()
  oglmxoutput<-list()
  fitinput<-list()
  fitinput$analhessian<-analhessian
  fitinput$sdmodel<-sdmodel
  fitinput$robust<-robust
  fitinput$link<-link
  oglmxoutput$link<-link
  oglmxoutput$sdmodel<-sdmodel
  oglmxoutput$call<-cl
  # if (!constantMEAN){formulaMEAN<-update(formulaMEAN,~0+.)}

  if (!is.null(formulaMEAN)) formulaMEAN <- as.formula(formulaMEAN)
  if (!is.null(formulaSD)) formulaSD <- as.formula(formulaSD)
  if (!is.null(selection)) selection <- as.formula(selection)

  if ((!is.null(selection)) && (!is.null(formulaSD))){
    message("Heteroskedastic model not implemented with censored model. Forcing homoskedastic model")
    formulaSD <- NULL
  }


  if (!is.null(formulaSD)){
    #  if (!constantSD){formulaSD<-update(formulaSD,~0+.)}
    cl$formulaMEAN<-mergeformulas(formulaMEAN,formulaSD)
  } else if (SameModelMEANSD){
    formulaSD<-formulaMEAN
  }


  if (!is.null(selection)){
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("selection", "data", "subset"), names(mf),
               0)
    mfS <- mf[c(1, m)]
    mfS$drop.unused.levels <- TRUE
    mfS$na.action <- na.pass
    mfS[[1]] <- as.name("model.frame")
    names(mfS)[2] <- "formula"
    mfS <- eval(mfS, parent.frame())
    mtS <- terms(mfS)
    Z <- model.matrix(mtS, mfS)
    y_selection <- model.response(mfS)

    mo <- match(c("formulaMEAN", "data", "subset"), names(mf),
                0)
    mfO <- mf[c(1, mo)]
    if (is.null(threshparam)) {
      mfO$drop.unused.levels <- TRUE
    }
    mfO$na.action <- na.pass
    mfO[[1]] <- as.name("model.frame")
    names(mfO)[2] <- "formula"
    mfO <- eval(mfO, parent.frame())
    mtO <- attr(mfO, "terms")
    X <- model.matrix(mtO, mfO)
    Y <- as.factor(model.response(mfO))


    # names(cls)[match("selection",names(cls))]<-"formula"
    # ms<-match(c("selection","data","subset","weights","na.action","offset"),names(cl),0L)
    # mfs<-cls[c(1L,ms)]
    # mfs[["formula"]] <- as.formula(mfs[["formula"]])
    # mfs$drop.unused.levels <- TRUE
    # mfs[[1L]] <- quote(stats::model.frame)
    # # terrible temporary hack
    # selection_model <- model.frame(mfs, data = dat)
    # y_selection <- selection_model[, as.character(mfs[["formula"]][[2]])]
    # Z <- selection_model[, !(colnames(selection_model) %in% as.character(mfs[["formula"]][[2]]))]
    # if (!grepl("0 +", as.character("selection"))) Z <- as.matrix(
    #   cbind.data.frame("(Intercept)" = 1L, Z)
    # )

    mf <- mfO

  } else{
    names(cl)[match("formulaMEAN",names(cl))]<-"formula"
    m<-match(c("formula","data","subset","weights","na.action","offset"),names(cl),0L)
    mf<-cl[c(1L,m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    mf<-eval(mf,parent.frame())

    factorvars<-names(attr(attr(mf,"terms"),"dataClasses"))[attr(attr(mf,"terms"),"dataClasses")=="factor"]
    attr(factorvars,"levels")<-lapply(factorvars,function(x){levels(mf[[x]])})
    oglmxoutput$factorvars<-factorvars

    mt<- attr(mf,"terms")
    Y<-as.factor(model.response(mf,"numeric"))
    X<-model.matrix(formulaMEAN,mf)


  }
  factorvarsX<-names(attr(X,"contrasts"))

  outcomenames<-levels(Y)
  oglmxoutput$Outcomes<-outcomenames



  if (!constantMEAN){
    Xint<-match("(Intercept)",colnames(X),nomatch = 0L)
    if (Xint>0L){X<-X[,-Xint,drop=FALSE]}
  }
  #termsMEAN<-terms(formulaMEAN)
  weights<-as.vector(model.weights(mf))

  oglmxoutput$NoVarModData<-data.frame(cbind(Y,weights))

  No.Outcomes<-nlevels(Y)
  if (No.Outcomes > 20 & !Force){
    stop("More than 20 different values for outcome variable.\n If you are sure you wish to estimate this model rerun command with Force option set to TRUE.")
  }

  Y<-as.numeric(Y)
  outcomeMatrix<-1 * ((col(matrix(0, length(Y), No.Outcomes))) == Y)
  colnames(outcomeMatrix)<-outcomenames

  oglmxoutput$NOutcomes<-No.Outcomes


  if (!is.null(formulaSD)){
    Z<-model.matrix(formulaSD,mf)
    #termsSD<-terms(formulaSD)
    oglmxoutput$Hetero<-TRUE
    if (!constantSD){
      Zint<-match("(Intercept)",colnames(Z),nomatch = 0L)
      if (Zint>0L){Z<-Z[,-Zint,drop=FALSE]}
    }
    oglmxoutput$selection<-FALSE
  } else if (is.null(selection)) {
    Z<-matrix(rep(1,nrow(X)),ncol=1)
    oglmxoutput$Hetero<-FALSE
    oglmxoutput$selection<-FALSE
  } else{
    # Z <- model.matrix(as.formula(selection),mf)
    # Z constructed later on
    oglmxoutput$Hetero<-FALSE
    oglmxoutput$selection<-TRUE
  }

  oglmxoutput$formula<-list(meaneq=formulaMEAN,sdeq=formulaSD,selection=selection)

  # beta
  if (!is.null(beta) & length(beta)==1){
    beta<-c(beta,rep(NA,ncol(X)-1))
  } else if (!is.null(beta) & length(beta)>1){
    # check that the specified vector is of correct length
    if (length(beta)!=ncol(X)){stop("Specified beta vector of incorrect length.")}
  } else if (is.null(beta)){
    beta<-rep(NA,ncol(X))
  }
  # delta
  if (!is.null(delta) & length(delta)==1){
    delta<-c(delta,rep(NA,ncol(Z)-1))
  } else if (!is.null(delta) & length(delta)>1){
    # check that the specified vector is of correct length
    if (length(delta)!=ncol(Z)){stop("Specified delta vector of incorrect length.")}
  } else if (is.null(delta)){
    delta<-rep(NA,ncol(Z))
  }
  # threshparam
  if (!is.null(threshparam) & length(threshparam)==1){
    threshparam<-c(threshparam,rep(NA,No.Outcomes-2))
  } else if (!is.null(threshparam) & length(threshparam)>1){
    # check that the specified vector is of correct length
    if (is.null(selection) && length(threshparam)!=No.Outcomes-1){stop("Specified vector of threshold parameters of incorrect length.")}
  } else if (is.null(threshparam)){
    threshparam<-rep(NA,No.Outcomes-1)
  }

  if (savemodelframe){
    oglmxoutput$modelframes<-list(X=X,Z=Z, outcomeMatrix = outcomeMatrix)
  }

  if (oglmxoutput$Hetero){
    namesX<-colnames(X)[colnames(X)!="(Intercept)"]
    namesZ<-colnames(Z)[colnames(Z)!="(Intercept)"]
    meanandvarNAME<-namesX[namesX %in% namesZ]
    meanandvarLOC<-match(meanandvarNAME,colnames(X))
    meanandvarLOCZ<-match(meanandvarNAME,colnames(Z))
    oglmxoutput$BothEq<-data.frame(meanandvarNAME,meanandvarLOC,meanandvarLOCZ,stringsAsFactors = FALSE)
  } else {oglmxoutput$BothEq<-NULL}

  # collect variable means and check which variables are binary.
  XVarMeans<-apply(X,2,mean)
  XVarBinary<-apply(X,2,.checkbinary)

  ZVarMeans<-apply(Z,2,mean)
  ZVarBinary<-apply(Z,2,.checkbinary)

  oglmxoutput$varMeans<-list(XVarMeans,ZVarMeans)
  oglmxoutput$varBinary<-list(XVarBinary,ZVarBinary)



  if (is.null(selection)){
    FitInput<-append(list(outcomeMatrix=outcomeMatrix,X=X,Z=Z,w=weights,beta=beta,delta=delta,threshparam=threshparam,
                          start=start,optmeth=optmeth, start_method = start_method, search_iter = search_iter),fitinput)
  } else{
    FitInput<- list(y=Y,y_selection = y_selection, X=X,Z=Z,thresholds=threshparam,
                    start=start,optmeth=optmeth, start_method = start_method, search_iter = search_iter,
                    gradient = gradient)
  }
  if (return_envir) return(FitInput)

  if (is.null(selection)){
    results<-append(oglmxoutput,do.call("oglmx.fit",FitInput))
  } else{
    results<-append(oglmxoutput,do.call("oglmx.fit2", FitInput))
    results$loglikelihood <- results$maximum
  }

  attr(results$loglikelihood,"No.Obs") <- length(Y)


  class(results)<-"oglmx"


  if (!is.null(selection)){
    names(results$estimate) <- c(colnames(Z), colnames(X), "log(sigma)", "atanhRho")
    results$coefAll <- c(results$estimate, sigma = unname(exp(results$estimate["log(sigma)"])),
                         sigmaSq = unname(exp(2 * results$estimate["log(sigma)"])),
                         rho = unname(tanh(results$estimate["atanhRho"])))
    jac <- cbind(diag(length(results$estimate)), matrix(0, length(results$estimate),
                                                        3))
    rownames(jac) <- names(results$estimate)
    colnames(jac) <- c(names(results$estimate), "sigma", "sigmaSq",
                       "rho")
    jac["log(sigma)", "sigma"] <- exp(results$estimate["log(sigma)"])
    jac["log(sigma)", "sigmaSq"] <- 2 * exp(2 * results$estimate["log(sigma)"])
    jac["atanhRho", "rho"] <- 1 - (tanh(results$estimate["atanhRho"]))^2
    results$vcov <- t(jac) %*% vcov(results) %*% jac


    results$params <- list(
      "selection" = 1:ncol(Z),
      "outcome"   = seq(from = ncol(Z) + 1, to = ncol(Z) + ncol(X)),
      "error" = seq(from = ncol(Z) + ncol(X) + 1,
                    to = length(results$coefAll)
      )
    )
    class(results) <- c("oglmx.selection",class(results))

  } else{
    results$vcov <- vcov(results,tol=tol)
  }


  return(results)

  #return(list(Y,X,Z,outcomeMatrix,weights))
}




