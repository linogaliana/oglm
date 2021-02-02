llk_selection <- function(y, beta, X, gamma, Z,
                          thresholds, rho, sigma){

  rho <- tanh(rho)
  sigma <- exp(sigma)
  zhat <- drop(Z %*% gamma)
  rho_matrix <- matrix(c(1, -rho, -rho, 1), nrow = 2)
  xhat <- drop(X %*% beta)
  idx_0 <- (y == 0)

  llk <- rep(NA, nrow(X))
  llk[idx_0] <- pnorm(-zhat, log.p = TRUE)[idx_0]

  if (is.finite(max(thresholds))) thresholds <- c(thresholds, Inf)
  m_index <- findInterval(y, thresholds)

  llk_observed_i <- function(i){
    p1 <- mvtnorm::pmvnorm(
      upper = c(
        ((thresholds[m_index + 1] - xhat)/sigma)[i],
        zhat[i]
      ),
      sigma = rho_matrix)
    p2 <- mvtnorm::pmvnorm(
      upper = c(
        ((thresholds[m_index] - xhat)/sigma)[i],
        zhat[i]
      ),
      sigma = rho_matrix)

    return(p1 - p2)
  }

  llk[!idx_0] <- log(
    pmax(
      do.call(rbind,
              lapply(which(!idx_0),
                     llk_observed_i)
      ), .Machine$double.eps)
  )

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


dllkdbeta <- function(y, beta, X, gamma, Z,
                      thresholds, rho, sigma){

  zhat <- drop(Z %*% gamma)
  rho_matrix <- matrix(c(1, -rho, -rho, 1), nrow = 2)
  xhat <- drop(X %*% beta)
  m_index <- findInterval(y, thresholds)
  idx_0 <- which(y == 0)

  grad_llk <- rep(0, nrow(X))

  grad_llk[idx_0] <- inverse_mills_ratio(- zhat) %*% Z

  p2 <- wrap_pnorm((thresholds[m + 1] - X %*% beta)/sigma,
                   z %*% gamma,
                   rho)
  p2 <- p2 - wrap_pnorm((thresholds[m] - X %*% beta)/sigma,
                        z %*% gamma,
                        rho)
  p2 <- p2*dnorm(z %*% gamma)*z
  p2_denom <- mvtnorm::pmvnorm((thresholds[m + 1] - beta %*% X)/sigma,
                               (thresholds[m] - beta %*% X)/sigma,
                               sigma = rho_matrix)

  grad_llk[!idx_0] <- p2/p2_denom


  return(p2)
}

inverse_mills_ratio <- function(x){
  dnorm(x)/pnorm(x)
}


wrap_pnorm <- function(x1,x2,rho){
  return(
    pnorm((x1 + rho*x2)/sqrt(1 - rho^2))
  )
}

