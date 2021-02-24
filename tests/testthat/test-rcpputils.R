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


dat$y <- as.numeric(dat$yS)
dat$y_outcome <- findInterval(dat$yOu, vec = bound)
dat$y_outcome[is.na(dat$y_outcome)] <- 1

fitInput <- oglm::oglmx(selection = "y ~ x1 + x2", yO ~ x1, data = dat,
                        threshparam = c(-Inf, 5, 15, Inf),
                        return_envir = TRUE)


thresholds = fitInput$thresholds
outcome_index = fitInput$y
y_selection <- fitInput$y_selection
idx_0 = y_selection == 0
xhat <- drop(fitInput$X%*%c(1,2))
zhat <- drop(fitInput$Z%*%c(1,2,-4))
rho_matrix <- matrix(c(1, -rho, -rho, 1), nrow = 2)
sigma = 2

testthat::expect_equal(
  oglm:::dff_dnorm(thresholds, outcome_index,
          xhat, zhat, sigma,
          rho_matrix),
  oglm:::dff_dnorm_cpp(thresholds, outcome_index,
          xhat, zhat, sigma,
          rho_matrix),
  check.attributes = FALSE
)
# microbenchmark::microbenchmark(
#   oglm:::dff_dnorm(thresholds, outcome_index,
#             xhat, zhat, sigma,
#             rho_matrix),
#   oglm:::dff_dnorm_cpp(thresholds, outcome_index,
#                 xhat, zhat, sigma,
#                 rho_matrix),
#   oglm:::dff_dnorm_cpp(thresholds, outcome_index,
#                 xhat, zhat, sigma,
#                 rho_matrix, ncores = 4)
# )
