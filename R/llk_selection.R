llk_selection <- function(y_selection, beta, X, gamma, Z,
                          y_outcome, thresholds, rho, sigma){

  rho_matrix <- matrix(c( 1, -rho, -rho, 1), nrow = 2 )

  llk_censor <- (1 - y_selection)*log(pnorm(-gamma %*% Z))
  llk_observed <- y_selection*log(
    mvtnorm::pmvnorm((thresholds[m + 1] - beta %*% X)/sigma,
                     (thresholds[m] - beta %*% X)/sigma,
                     sigma = rho_matrix)
  )

  llk <- llk_censor + llk_observed

  return(llk)
}


dens_mvnorm <- function(x1,x2,rho){

  rho_matrix <- matrix(c( 1, rho, rho, 1), nrow = 2 )
  dens <- mvtnorm::dmvnorm(c(x1, x2), sigma = rho_matrix)

  return(dens)
}

dPhidx2 <- function(x1,x2,rho){

  return(
    pnorm((x1 - rho*x2)/sqrt(1-rho^2))*dnorm(x2)
  )
}

dPhidx1 <- function(x1,x2,rho){

  return(
    pnorm((x2 - rho*x1)/sqrt(1-rho^2))*dnorm(x1)
  )
}


dPhidrho <- function(x1,x2,rho){
  return(
    dens_mvnorm(x1,x2,rho)
  )
}


dPhidgamma <- function(x, z, rho, sigma,
                                beta, gamma, alpha_m){

  dphi_eval <- dPhidx2((alpha_m - x %*% beta) /sigma, z %*% gamma, -rho)

  return(
    t(dphi_eval) %*% z
  )

}


dPhidbeta <- function(x, z, rho, sigma,
                       beta, gamma, alpha_m){


  dphi_eval <- dPhidx1((alpha_m - x %*% beta) /sigma, z %*% gamma, -rho)

  return(
    -t(dphi_eval) %*% x/sigma
  )

}
