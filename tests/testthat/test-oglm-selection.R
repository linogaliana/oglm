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

requireNamespace("sampleSelection", quietly = TRUE)

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


  testthat::expect_equal(
    llk_selection_wrapper(
      theta = selection_model$estimate,
      y = dat$y_outcome,
      y_selection = dat$y,
      X = fitInput$X,
      Z = fitInput$Z,
      thresholds = bound),
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




# OPTIMIZATION IS OK -------------------

opt_maxLik <- maxLik::maxLik(
  llk_selection_wrapper,
  grad = grad_llk_selection_wrapper,
  hess = NULL,
  method = "BHHH",
  start = selection_model$start,
  y = dat$y_outcome,
  y_selection = dat$y,
  X = as.matrix(fitInput$X), Z = as.matrix(fitInput$Z),
  thresholds = bound
)

testthat::test_that("maxLik optimization consistent with sampleSelection output", {

  testthat::expect_equal(
    opt_maxLik$maximum,
    selection_model$maximum
  )
  testthat::expect_equal(
    opt_maxLik$estimate,
    selection_model$estimate
  )
  testthat::expect_equal(
    opt_maxLik$gradient,
    selection_model$gradient
  )
  testthat::expect_equal(
    opt_maxLik$hessian,
    selection_model$hessian
  )
  testthat::expect_equal(
    opt_maxLik$code,
    selection_model$code
  )
  testthat::expect_equal(
    opt_maxLik$message,
    selection_model$message
  )
  testthat::expect_equal(
    opt_maxLik$iterations,
    selection_model$iterations
  )
  testthat::expect_equal(
    opt_maxLik$type,
    selection_model$type
  )
  testthat::expect_equal(
    opt_maxLik$gradientObs,
    selection_model$gradientObs
  )
  testthat::expect_equal(
    opt_maxLik$objectiveFn,
    llk_selection_wrapper
  )
}
)



# oglm.fit2 function ------------------------


# fit_output <- oglmx.fit2(y = dat$y_outcome,
#                          y_selection = dat$y,
#                          X = as.matrix(fitInput$X), Z = as.matrix(fitInput$Z),
#                          start = selection_model$start,
#                          thresholds = bound)
#
# opt_maxLik$start <- selection_model$start
#
# testthat::test_that("Wrapper does not change anything", {
#   testthat::expect_equal(
#     fit_output,
#     opt_maxLik
#   )
# })


# integration in oglm ------------------------


output_oglm <- oglm::oglmx(selection = "y ~ x1 + x2", yO ~ x1, data = dat,
                           threshparam = c(-Inf, 5, 15, Inf),
                           start = selection_model$start)


testthat::expect_s3_class(
  output_oglm, "oglmx.selection"
)
testthat::expect_s3_class(
  output_oglm, "oglmx"
)

testthat::test_that("Reported formulas correct",{
  testthat::expect_equal(class(output_oglm$formula),
                         "list")
  testthat::expect_equal(output_oglm$formula$meaneq,
                         as.formula("yO ~ x1"))
  testthat::expect_equal(output_oglm$formula$selection,
                         as.formula("y ~ x1 + x2"))
  testthat::expect_equal(output_oglm$formula$sdeq,
                         NULL)
})


testthat::test_that("oglm master function consistent with sampleSelection output", {

  testthat::expect_equal(
    output_oglm$maximum,
    selection_model$maximum
  )
  testthat::expect_equal(
    as.numeric(output_oglm$estimate),
    as.numeric(selection_model$estimate)#,
    #    tolerance = 1e-4
  )
  testthat::expect_equal(
    as.numeric(output_oglm$gradient),
    as.numeric(selection_model$gradient)#,
    #    tolerance = 1e-3
  )
  testthat::expect_equal(
    output_oglm$hessian,
    selection_model$hessian,
    check.attributes = FALSE#,
    #    tolerance = 1e-3
  )
  testthat::expect_equal(
    output_oglm$code,
    selection_model$code
  )
  testthat::expect_equal(
    output_oglm$message,
    selection_model$message
  )
  testthat::expect_equal(
    output_oglm$iterations,
    selection_model$iterations
  )
  testthat::expect_equal(
    output_oglm$type,
    selection_model$type
  )
  testthat::expect_equal(
    output_oglm$gradientObs,
    selection_model$gradientObs
  )
  testthat::expect_equal(
    output_oglm$objectiveFn,
    llk_selection_wrapper
  )
})


testthat::test_that("Variance-covariance matrix is the same", {
  testthat::expect_equal(
    output_oglm$vcov,
    selection_model$vcovAll,
    check.attributes = FALSE#,
#    tolerance = 1e-2
  )
})


# OPTIONS ---------

testthat::expect_error(
  oglm::oglmx(selection = "y ~ x1 + x2", yO ~ x1, data = dat,
              threshparam = c(-Inf, 5, 15, Inf),
              start = selection_model$start,
              gradient = "nonauthorized")
)

testthat::test_that("Default to analytical",{

  out1 <- oglm::oglmx(selection = "y ~ x1 + x2", yO ~ x1, data = dat,
                      threshparam = c(-Inf, 5, 15, Inf),
                      start = selection_model$start,
                      gradient = "analytical")
  out2 <- oglm::oglmx(selection = "y ~ x1 + x2", yO ~ x1, data = dat,
                      threshparam = c(-Inf, 5, 15, Inf),
                      start = selection_model$start)
  testthat::expect_equal(
    out1$estimate,
    out2$estimate
  )
})

testthat::test_that("Results when using numerical gradient", {
  out3 <- oglm::oglmx(selection = "y ~ x1 + x2", yO ~ x1, data = dat,
                      threshparam = c(-Inf, 5, 15, Inf),
                      start = selection_model$start,
                      gradient = "numerical")

  testthat::expect_equal(
    out2$estimate,
    out3$estimate
  )

})


# compare with stata ---------------------


# womensat <- haven::read_dta("https://www.stata-press.com/data/r16/womensat.dta")
#
# model_stata <- oglm::oglmx(selection = "work ~ education + age + factor(married) + factor(children)",
#             "satisfaction ~ education + age", data = womensat,
#             threshparam = c(-Inf, 1.728757, 2.64357,
#                             3.642911, Inf))



# microbenchmark::microbenchmark(
#   oglm::oglmx(selection = "y ~ x1 + x2", yO ~ x1, data = dat,
#               threshparam = c(-Inf, 5, 15, Inf),
#               start = selection_model$start),
#   sampleSelection::selection( yS ~ x1 + x2, yO ~ x1, data = dat, boundaries = bound,
#                               start = selection_model$start),
#   times = 10L
# )



# profvis::profvis({
#   x1 = oglm::oglmx(selection = "y ~ x1 + x2", yO ~ x1, data = dat,
#                              threshparam = c(-Inf, 5, 15, Inf),
#                              start = selection_model$start)})
# profvis::profvis({
#   x2 = sampleSelection::selection( yS ~ x1 + x2, yO ~ x1, data = dat, boundaries = bound,
#                               ys = TRUE, xs = TRUE, yo = TRUE, xo = TRUE)
# })
