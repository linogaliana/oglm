nObs <- 300
betaS <- c( 1, 1, -1 )
betaO <- c( 10, 4 )
rho <- 0.4
sigma <- 5
# boundaries of the intervals
bound <- c(-Inf,5,15,Inf)
# set 'seed' of the pseudo random number generator
# in order to always generate the same pseudo random numbers
set.seed(123)
# generate variables x1 and x2
dat <- data.frame( x1 = rnorm( nObs ), x2 = rnorm( nObs ) )
# variance-covariance matrix of the two error terms
vcovMat <- matrix( c( 1, rho*sigma, rho*sigma, sigma^2 ), nrow = 2 )
# generate the two error terms
eps <- mvtnorm::rmvnorm( nObs, sigma = vcovMat )
dat$epsS <- eps[,1]
dat$epsO <- eps[,2]
# generate the selection variable
dat$yS <- with( dat, betaS[1] + betaS[2] * x1 + betaS[3] * x2 + epsS ) > 0
table( dat$yS )
# generate the unobserved/latent outcome variable
dat$yOu <- with( dat, betaO[1] + betaO[2] * x1 + epsO )
dat$yOu[ !dat$yS ] <- NA
# obtain the intervals of the outcome variable
dat$yO <- cut( dat$yOu, bound )
table( dat$yO )

selection_model <- sampleSelection::selection( yS ~ x1 + x2, yO ~ x1, data = dat, boundaries = bound,
                                               ys = TRUE, xs = TRUE, yo = TRUE, xo = TRUE)



# TERMS ---------------

dat$y <- as.numeric(dat$yS)
dat$y_outcome <- findInterval(dat$yOu, vec = bound)
dat$y_outcome[is.na(dat$y_outcome)] <- 1

fitInput <- oglm::oglmx(selection = "y ~ x1 + x2", yO ~ x1, data = dat,
                        threshparam = c(-Inf, 5, 15, Inf),
                        return_envir = TRUE)

testthat::test_that("Same matrices constructed by our package and sample selection", {
  testthat::expect_equal(
    data.frame(fitInput$X[fitInput$y_selection == 1,]),
    data.frame(na.omit(selection_model$xo))
  )
  testthat::expect_equal(
    data.frame(fitInput$Z),
    data.frame(na.omit(selection_model$xs))
  )
  testthat::expect_equal(
    as.numeric(fitInput$y[fitInput$y_selection == 1]),
    as.numeric(na.omit(selection_model$yo))
  )
  testthat::expect_equal(
    as.numeric(fitInput$y_selection),
    as.numeric(na.omit(selection_model$ys))
  )
})




# LOGLIKELIHOOD ----------------------------


testthat::test_that("Log-likelihood consistent with sampleSelection",{
  # Total
  testthat::expect_equal(
    sum(llk_selection(y = dat$y_outcome,
                      y_selection = dat$y,
                      beta = selection_model$estimate[4:5],
                      X = fitInput$X,
                      gamma = selection_model$estimate[1:3],
                      Z = fitInput$Z,
                      thresholds = bound,
                      rho = selection_model$estimate["atanhRho"],
                      sigma = selection_model$estimate['logSigma'])
    ),
    as.numeric(logLik(selection_model))
  )
  # Term by term
  testthat::expect_equal(
    llk_selection(y = dat$y_outcome,
                  y_selection = dat$y,
                  beta = selection_model$estimate[4:5],
                  X = fitInput$X,
                  gamma = selection_model$estimate[1:3],
                  Z = fitInput$Z,
                  thresholds = bound,
                  rho = selection_model$estimate["atanhRho"],
                  sigma = selection_model$estimate['logSigma']),
    as.numeric(selection_model$objectiveFn(selection_model$estimate))
  )


})




# LOG LIKELIHOOD GRADIENTS -----------------


testthat::test_that("Gradient consistent with sampleSelection",{
  # Total
  testthat::expect_equal(
    colSums(grad_llk_selection(y = dat$y_outcome,
                               y_selection = dat$y,
                               beta = selection_model$estimate[4:5],
                               X = fitInput$X,
                               gamma = selection_model$estimate[1:3],
                               Z = fitInput$Z,
                               thresholds = bound,
                               rho = selection_model$estimate["atanhRho"],
                               sigma = selection_model$estimate['logSigma'])
    ),
    as.numeric(selection_model$gradient))

  # Term by term
  testthat::expect_equal(
    grad_llk_selection(y = dat$y_outcome,
                       y_selection = dat$y,
                       beta = selection_model$estimate[4:5],
                       X = fitInput$X,
                       gamma = selection_model$estimate[1:3],
                       Z = fitInput$Z,
                       thresholds = bound,
                       rho = selection_model$estimate["atanhRho"],
                       sigma = selection_model$estimate['logSigma']),
    attr(selection_model$objectiveFn(selection_model$estimate), "gradient")
  )

})


theta <- selection_model$estimate

testthat::test_that("Derivates for loglikelihood evaluated around beta", {


  grad_numeric <- maxLik::numericGradient(
    function(theta){llk_selection_wrapper(theta, y = dat$y_outcome,
                                          y_selection = dat$y,
                                          X = fitInput$X, Z = fitInput$Z,
                                          thresholds = bound)},
    theta
  )


  grad_analytical <- grad_llk_selection_wrapper(theta, y = dat$y_outcome,
                                                y_selection = dat$y,
                                                X = fitInput$X, Z = fitInput$Z,
                                                thresholds = bound)

  testthat::expect_equal(
    grad_numeric,
    grad_analytical,
    check.attributes = FALSE
  )

})


})


y <- dat$yOu
y[is.na(y)] <- 0

X <- as.matrix(cbind(1, dat$x1))
Z <- as.matrix(cbind(1, dat[,c("x1", "x2")]))
beta <-  as.numeric(res$estimate[4:5])
gamma <- as.numeric(res$estimate[1:3])
sigma <- as.numeric(res$estimate[6])
rho <- res$estimate[7]
thresholds <- bound



# INTEGRATION IN OGLM MASTER FUNCTION -------------------


fitInput <- sampleSelection::selection(selection =  yS ~ x1 + x2,
                                       outcome = yO ~ x1, data = dat,
                                       boundaries = bound)
dat$y <- dat$yOu
dat$y[is.na(dat$y)] <- 0

dat$yO[is.na(dat$yO)] <- "(-Inf,5]"

oglm::oglmx("yO ~ x1 + x2",
            data=dat,link="probit",constantMEAN=FALSE,
            constantSD=FALSE,delta=0,
            threshparam = c(5, 15))




fitInput <- oglmx(selection =  "y ~ x1 + x2",
                  formulaMEAN ="yO ~ x1", data = dat,
                  threshparam = c(5, 15))







y <- dat$yOu
y[is.na(y)] <- 0

X <- as.matrix(cbind(1, dat$x1))
Z <- as.matrix(cbind(1, dat[,c("x1", "x2")]))
beta <-  as.numeric(res$estimate[4:5])
gamma <- as.numeric(res$estimate[1:3])
sigma <- as.numeric(res$estimate[6])
rho <- res$estimate[7]
thresholds <- bound
boundaries = bound








selection_model <- sampleSelection::selection( yS ~ x1 + x2, yO ~ x1, data = dat, boundaries = bound )

start <- toto( yS ~ x1 + x2, yO ~ x1, data = dat, boundaries = bound )



head(
  grad_llk_selection(y, selection_model$estimate[4:5], X,
                     selection_model$estimate[1:3], Z,
                     thresholds, selection_model$estimate["atanhRho"], selection_model$estimate['logSigma'])
)

head(
  selection_model$gradientObs
)

theta <- c(beta, gamma, rho, sigma)

opt <- maxLik::maxLik(
  llk_selection_wrapper,
  grad = grad_llk_selection_wrapper,
  hess = NULL,
  start = startVal,
  y = y, X = X, Z = Z, thresholds = thresholds
)


