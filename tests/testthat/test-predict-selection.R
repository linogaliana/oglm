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
oglm_model <- oglm::oglmx(selection = "y ~ x1 + x2", yO ~ x1, data = dat,
                           threshparam = c(-Inf, 5, 15, Inf),
                           start = selection_model$start)

library(sampleSelection)

testthat::test_that("linear predictors of the selection equation", {
  testthat::expect_equal(
    predict(selection_model, type = "link", part = "selection"),
    predict(oglm_model, model = "selection", type = "xb", newdata = dat)
  )
})
testthat::test_that("predicted probabilities of the selection equation", {
  testthat::expect_equal(
    as.numeric(predict(selection_model, type = "response", part = "selection")),
    as.numeric(predict(oglm_model, model = "selection", type = "probs", newdata = dat))
  )
  testthat::expect_equal(
    as.numeric(predict(selection_model, type = "response", part = "selection")),
    as.numeric(predict(oglm_model, type = "P[y == 0|Z]", newdata = dat))
  )
})
testthat::test_that("unconditional expectations", {
  testthat::expect_equal(
    predict(selection_model, type = "unconditional")[isTRUE(dat$yS)],
    predict(oglm_model, model = "outcome", type = "xb", newdata = dat)[isTRUE(dat$yS)]
  )
  testthat::expect_equal(
    predict(selection_model, type = "unconditional")[isTRUE(dat$yS)],
    predict(oglm_model, type = "E[y|X]", newdata = dat)[isTRUE(dat$yS)]
  )
  testthat::expect_equal(
    predict(selection_model, type = "unconditional"),
    predict(oglm_model, type = "E[y|X,y>0]", newdata = dat)
  )
})




