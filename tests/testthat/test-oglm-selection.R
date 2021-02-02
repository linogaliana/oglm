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

res <- sampleSelection::selection( yS ~ x1 + x2, yO ~ x1, data = dat, boundaries = bound )

mf <- sel(selection =  yS ~ x1 + x2, outcome = yO ~ x1, data = dat, boundaries = bound )




y <- dat$yOu
y[is.na(y)] <- 0

X <- as.matrix(cbind(1, dat$x1))
Z <- as.matrix(cbind(1, dat[,c("x1", "x2")]))
beta <-  as.numeric(res$estimate[4:5])
gamma <- as.numeric(res$estimate[1:3])
sigma <- as.numeric(res$estimate[6])
rho <- res$estimate[7]
thresholds <- bound





testthat::test_that("Log-likelihood consistent with sampleSelection",{
  testthat::expect_equal(
    sum(llk_selection(y, beta, X, gamma, Z,
                      thresholds, rho, sigma)
    ),
    as.numeric(logLik(res))
  )
})
