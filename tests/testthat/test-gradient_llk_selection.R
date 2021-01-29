

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


