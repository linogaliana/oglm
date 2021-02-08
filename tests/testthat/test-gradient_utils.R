

x1 <- rnorm(1)
x2 <- rnorm(1)
rho <- runif(1)

testthat::test_that("Wrapper for multivariate density", {
  testthat::expect_equal(
    dens_mvnorm(x1,x2,rho),
    exp(
      - (x1^2 - 2*rho*x1*x2 + x2^2)/(2*(1-rho^2))
    )/(2*pi*sqrt(1 - rho^2))
  )
})

testthat:::test_that("Numerical derivatives for multivariate normal ok",{

  grad_x1 <- maxLik::numericGradient(
    function(x){mvtnorm::pmvnorm(upper = c(x,x2) ,
                                 sigma = matrix(c(1, rho, rho, 1),
                                                nrow = 2))},
    x1
  )
  grad_x2 <- maxLik::numericGradient(
    function(y){mvtnorm::pmvnorm(upper = c(x1,y) ,
                                 sigma = matrix(c(1, rho, rho, 1),
                                                nrow = 2))},
    x2
  )
  grad_rho <- maxLik::numericGradient(
    function(r){mvtnorm::pmvnorm(upper = c(x1,x2) ,
                                 sigma = matrix(c(1, r, r, 1),
                                                nrow = 2))},
    rho
  )


  testthat::expect_equal(
    as.numeric(grad_x1),
    dPhidx1(x1,x2,rho)
  )
  testthat::expect_equal(
    as.numeric(grad_x2),
    dPhidx2(x1,x2,rho)
  )
  testthat::expect_equal(
    as.numeric(grad_rho),
    dPhidrho(x1,x2,rho)
  )
})


testthat::test_that("Derivates for multivariate normal evaluated around gamma", {

  X <- matrix(rnorm(30),10,3)
  Z <- matrix(rnorm(20),10,2)
  alpha <- 2
  beta <- c(0,-1,2)
  gamma <- c(2,-0.5)
  sigma <- 3

  grad_gamma <- function(i){
    maxLik::numericGradient(
      function(g){mvtnorm::pmvnorm(upper = c((alpha  - X %*% beta)[i,]/sigma,
                                             Z[i,] %*% g),
                                   sigma = matrix(c(1, rho, rho, 1),
                                                  nrow = 2))},
      gamma
    )
  }

  testthat::expect_equal(
    do.call(rbind, lapply(1:10, grad_gamma)),

    do.call(rbind, lapply(1:10, function(i){
      dPhidgamma(X[i,], Z[i,], -rho, sigma,
                 beta, gamma, alpha_m = alpha)
    }))
  )

})

X <- matrix(rnorm(30),10,3)
Z <- matrix(rnorm(20),10,2)
alpha <- 2
beta <- c(0,-1,2)
gamma <- c(2,-0.5)
sigma <- 3


testthat::test_that("Derivates for multivariate normal evaluated around beta", {



  grad_beta <- function(i){
    maxLik::numericGradient(
      function(b){mvtnorm::pmvnorm(upper = c((alpha  - X %*% b)[i,]/sigma,
                                             Z[i,] %*% gamma),
                                   sigma = matrix(c(1, rho, rho, 1),
                                                  nrow = 2))},
      beta
    )
  }

  testthat::expect_equal(
    do.call(rbind, lapply(1:10, grad_beta)),

    do.call(rbind, lapply(1:10, function(i){
      dPhidbeta(X[i,], Z[i,], -rho, sigma,
                 beta, gamma, alpha_m = alpha)
    }))
  )

})




