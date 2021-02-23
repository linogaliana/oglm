#' Predictions from ordered discrete model
#'
#' @inheritParams stats::predict.glm
#' @param object An ordered discrete model computed using `oglmx` package
#' @param type Prediction considered. Can be \emph{'class'}
#'  (predicted label as maximum probability label, default);
#'  \emph{'probs'} (probabilities of each observation to belong to a particular class)
#'  or \emph{'latent'} (\eqn{Y^*} latent variable)
#' @return Depends of the value of prediction \code{type} \describe{
#'   \item{class}{Most likely label, i.e. \eqn{l = \arg \max_l p_j}}
#'   \item{probs}{Prediced probabilities for each label, i.e. \eqn{p_ij} matrix with
#'    $i$ observation index and $j$ class label}
#'   \item{latent}{Prediced value in latent space \eqn{y^*}}
#' }
#'
#'
#' @importFrom stats .checkMFClasses delete.response
#' @importFrom stats predict
#' @importFrom stats model.frame model.matrix napredict pcauchy
#' @importFrom stats plogis pnorm rcauchy rlogis rnorm terms
#' @export


predict.oglmx <- function(object, newdata = NULL, type = c("class", "probs","latent","xb"), ...){

  # CHECK IF oglmx OBJECT
  # --------------------------
  if (!inherits(object, "oglmx")) stop("not a \"oglmx\" object")

  # CHECK PREDICTION TYPE
  # --------------------------
  type <- match.arg(type)

  # FIT IF missing(newdata)
  # ---------------------------

  if (missing(newdata)){
    #Y <- object$fitted
    stop('fit method not yet implemented: use predict')
  }

  newdata <- as.data.frame(newdata)
  formula <- object$formula$`meaneq`


  # CREATE TERMS TO REPLICATE MASS::predict.polr BEHAVIOR
  # --------------------------------------------------

  # Transform formula in terms
  object$terms <- terms(as.formula(formula))

  # Keep only covariates
  Terms <- delete.response(object$terms)

  # Covariates matrix
  m <- model.frame(Terms, newdata, na.action = function(x) x,
                   xlev = object$factorvars)

  # Check factors are ok
  if (!is.null(cl <- attr(Terms, "dataClasses")))
    .checkMFClasses(cl, m)

  X <- model.matrix(Terms, m, contrasts = object$contrasts)

  xnocol <- match(names(object$coefficients),
                  colnames(X),
                  nomatch = 0L)
  if (length(xnocol)>0L){
    X2 <- X[, xnocol]
  }


  # Finalize covariates matrix by removing intercept
  xint <- match("(Intercept)", colnames(X), nomatch = 0L)
  if (xint > 0L){
    X <- X[, -xint, drop = FALSE]
  }

  # Parameters for logistic/normal/... distribution
  n <- nrow(X)
  q <- length(object$allparams$threshparam)


  # MAKE PREDICTION
  # -------------------------------

  coeff_list <- object$coefficients

  # REMOVE THRESHOLD COEFFICIENTS
  if (sum(grepl("Threshold", names(coeff_list)))>0){
    coeff <- coeff_list[-grep("Threshold", names(coeff_list))]
  }else{
    coeff <- coeff_list
  }

  # REMOVE STANDARD ERROR RELATED PARAMETER
  coeff <- coeff[attr(coeff_list, "coefftypes")[[1]]]
  coeff <- coeff[!is.na(names(coeff))]
  coeff <- coeff[names(coeff) != "ln(sigma)"]


  if (ncol(X2) != length(coeff)){
    X2 <- X2[,attr(coeff_list, "coefftypes")[[1]]]
  }
  # X2 <- X2[,attr(coeff_list, "coefftypes")[[1]]]


  # POSSIBILITY THAT SOME FIXED EFFECTS ARE NOT PRESENT IN X2
  # REMOVE THEM FROM coeff
  coeff2 <- coeff[names(coeff) %in% colnames(X2)]

  # x*\beta vector (nb: intercept column should be added back)
  eta <- drop(X2 %*% coeff2)

  if (type == "xb") return(eta)

  # if type == "latent", we simulate epsilon
  if (type == "latent"){

    epsilon_distribution <- switch(object$link, logit = rlogis, probit = rnorm,
                                   loglog = rgumbel, cloglog = rGumbel, cauchit = rcauchy)


    # sigma identified to 1 in case of non-user defined interval regression
    if (sum(grepl("Threshold", names(coeff_list)))>0){

      simulated_epsilon <- epsilon_distribution(length(eta))
      sigma_vector <- rep(1L, times = length(eta)) # A vérifier

    } else{

      sigma_vector <- sigma(object, newdata = newdata)

      if (object$link == "logit"){
        m <- 0
        sigma_vector <- sqrt(3)*sigma_vector/pi #Variance for logistic is s^2 * pi^2/3
        simulated_epsilon <- sapply(sigma_vector, function(s) epsilon_distribution(1L, location = 0,
                                                                                   scale = s)
        )
      } else if (object$link == "probit"){
        simulated_epsilon <- sapply(sigma_vector, function(s) epsilon_distribution(1L, mean = 0,
                                                                                   sd = s)
        )
      } else{
        # Other distributions for epsilon
        simulated_epsilon <- epsilon_distribution(length(eta))
      }

    }

    eta <- eta + simulated_epsilon

  }

  # Which distribution should be applied ?
  pfun <- switch(object$link, logit = plogis, probit = pnorm,
                 loglog = pgumbel, cloglog = pGumbel, cauchit = pcauchy)

  # Transform from latent space y = x*beta to probabilities
  cumpr <- matrix(pfun(matrix(object$allparams$threshparam, n, q, byrow = TRUE) -
                         eta), , q)

  # CREATE PREDICTION
  # ---------------------------------

  Y <- t(apply(cumpr, 1L, function(x) diff(c(0, x, 1))))
  dimnames(Y) <- list(rownames(X), object$Outcomes)


  if (missing(newdata) && !is.null(object$na.action)){
    Y <- napredict(object$na.action, Y)
  }


  if (type == "class"){
    factor(max.col(Y), levels = seq_along(object$Outcomes),
           labels = object$Outcomes)
  }else if (type == "probs") {
    drop(Y)
  } else{
    # When type = "latent"
    return(list("xb" = eta - simulated_epsilon, "y_latent_pred" = eta,
                "sigma_vector" = sigma_vector))
  }
}


#' @export
predict.oglmx.selection <- function(object, newdata = NULL,
                                    type = c("class", "probs","latent","xb", "E[y|X]", "P[y == 0|Z]", "E[y|X,y>0]", "E[y|X,Z]"),
                                    model = c("both", "outcome", "selection"),
                                    threshold_proba_selection = .5,
                                    ...){

  # CHECK IF oglmx OBJECT
  # --------------------------
  if (!inherits(object, "oglmx")) stop("not a \"oglmx\" object")


  # CHECK PREDICTION TYPE
  # --------------------------
  type <- match.arg(type)
  model <- match.arg(model)


  # FIT IF missing(newdata)
  # ---------------------------

  if (missing(newdata)){
    #Y <- object$fitted
    stop('fit method not yet implemented: use predict')
  }

  newdata <- as.data.frame(newdata)
  formula_outcome <- object$formula$`meaneq`
  formula_selection <- object$formula$selection


  # CREATE TERMS TO REPLICATE MASS::predict.polr BEHAVIOR
  # --------------------------------------------------

  # Transform formula in terms
  object$terms_outcome <- terms(as.formula(formula_outcome))
  object$terms_selection <- terms(as.formula(formula_selection))

  # Keep only covariates
  Terms_outcome <- delete.response(object$terms_outcome)
  Terms_selection <- delete.response(object$terms_selection)

  # Covariates matrix
  mo <- model.frame(Terms_outcome, newdata, na.action = function(x) x,
                   xlev = object$factorvars)
  ms <- model.frame(Terms_selection, newdata, na.action = function(x) x,
                    xlev = object$factorvars)

  # selection
  y_selection <- model.response(model.frame(object$terms_selection,
                             newdata, na.action = function(x) x))


  # Check factors are ok
  if (!is.null(cl <- attr(Terms_outcome, "dataClasses")))
    .checkMFClasses(cl, mo)
  if (!is.null(cl <- attr(Terms_selection, "dataClasses")))
    .checkMFClasses(cl, ms)

  X <- model.matrix(Terms_outcome, mo, contrasts = object$contrasts)
  Z <- model.matrix(Terms_selection, ms, contrasts = object$contrasts)

  # xnocol <- match(names(object$coefficients),
  #                 colnames(X),
  #                 nomatch = 0L)
  # if (length(xnocol)>0L){
  #   X2 <- X[, xnocol]
  # }
  #
  #
  # # Finalize covariates matrix by removing intercept
  # xint <- match("(Intercept)", colnames(X), nomatch = 0L)
  # if (xint > 0L){
  #   X <- X[, -xint, drop = FALSE]
  # }

  # Parameters for logistic/normal/... distribution
  n <- nrow(X)
  q <- length(object$allparams$threshparam)


  # MAKE PREDICTION
  # -------------------------------

  coeff_list <- object$coefAll

  # REMOVE THRESHOLD COEFFICIENTS
  if (sum(grepl("Threshold", names(coeff_list)))>0){
    coeff <- coeff_list[-grep("Threshold", names(coeff_list))]
  }else{
    coeff <- coeff_list
  }

  x_params <- object$params$outcome
  z_params <- object$params$selection

  beta <- coeff_list[x_params]
  gamma <- coeff_list[z_params]
  rho <- as.numeric(coeff_list['rho'])
  sigma <- as.numeric(coeff_list['sigma'])

  zhat <- drop(Z %*% gamma)
  xhat <- drop(X %*% beta)
  if (type == "E[y|X,y>0]") xhat[y_selection == 0] <- NA



  epsilon_distribution <- switch(object$link, logit = rlogis, probit = rnorm,
                                 loglog = rgumbel, cloglog = rGumbel, cauchit = rcauchy)
  pfun <- switch(object$link, logit = plogis, probit = pnorm,
                                 loglog = pgumbel, cloglog = pGumbel, cauchit = pcauchy)


  if ((type == "xb" && model == "outcome") || (type %in% c("E[y|X]", "E[y|X,y>0]"))) return(xhat)
  if ((type == "xb" && model == "selection")) return(zhat)
  if ((type == "xb" && model == "both")) return(list("E[y|X]" = xhat,
                                                     "gammaZ" = zhat))
  if ((type == "probs" && model == "selection") || (type == "P[y == 0|Z]")) return(pfun(zhat))

  if (type == "E[y|X,Z]") return(
    cbind("E[y|X,Z,y == 0]" = xhat - rho*sigma*inverse_mills_ratio(-zhat),
          "E[y|X,Z,y > 0]" = xhat + rho*sigma*inverse_mills_ratio(zhat)
    )
  )

  if ((type == "class") && (model != "outcome")){
    class_selection <- as.numeric(pfun(zhat)>threshold_proba_selection)
    if (model == "selection") return(class_selection)
  }


  # PROBS WHEN USING BOTH OR SELECTION ----------------
  # must account for the error terms correlation



  rho_matrix <- matrix(c(1, -rho, -rho, 1), nrow = 2)

}
