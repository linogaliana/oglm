llk_selection <- function(y, beta, X, gamma, Z,
                          thresholds, rho, sigma){

  # tanh rho is estimated
  rho <- tanh(rho)

  # exp(sigma) estimated to get positive sigma
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


dllkdgamma <- function(y, beta, X, gamma, Z,
                      thresholds, rho, sigma){

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


  # tanh rho is estimated
  rho <- tanh(rho)

  # exp(sigma) estimated to get positive sigma
  sigma <- exp(sigma)

  zhat <- drop(Z %*% gamma)
  rho_matrix <- matrix(c(1, -rho, -rho, 1), nrow = 2)
  xhat <- drop(X %*% beta)
  m_index <- findInterval(y, thresholds)
  idx_0 <- (y == 0)

  grad_llk <- matrix(0, nrow(Z), length(gamma))

  # Gradient for observations where y == 0
  grad_llk[idx_0,] <- - (inverse_mills_ratio(- zhat) * Z)[idx_0,]

  p2 <- wrap_pnorm((thresholds[m_index + 1] - xhat)/sigma,
                   zhat,
                   rho)
  p2 <- p2 - wrap_pnorm((thresholds[m_index] - xhat)/sigma,
                        zhat,
                        rho)
  p2 <- p2*dnorm(zhat)*Z

  diff_pmvnorm <- pmax(
      do.call(rbind,
              lapply(which(!idx_0),
                     llk_observed_i)
      ), .Machine$double.eps)

  grad_llk[!idx_0,] <- p2[!idx_0,]/diff_pmvnorm


  return(grad_llk)
}


dllkdbeta <- function(y, beta, X, gamma, Z,
                       thresholds, rho, sigma){

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


  # tanh rho is estimated
  rho <- tanh(rho)

  # exp(sigma) estimated to get positive sigma
  sigma <- exp(sigma)

  zhat <- drop(Z %*% gamma)
  rho_matrix <- matrix(c(1, -rho, -rho, 1), nrow = 2)
  xhat <- drop(X %*% beta)
  m_index <- findInterval(y, thresholds)
  idx_0 <- (y == 0)

  grad_llk <- matrix(0, nrow(X), length(beta))

  # Gradient for observations where y == 0
  # stays 0 !

  p2_pnorm1 <- wrap_pnorm(zhat,
                   (thresholds[m_index + 1] - xhat)/sigma,
                   rho)
  p2_dnorm1 <- dnorm((thresholds[m_index + 1] - xhat)/sigma)
  p2_pnorm2 <- wrap_pnorm(zhat,
                        (thresholds[m_index] - xhat)/sigma,
                        rho)
  p2_dnorm2 <- dnorm((thresholds[m_index] - xhat)/sigma)
  p2 <- p2_pnorm1*p2_dnorm1 -  p2_pnorm2*p2_dnorm2
  p2 <- -p2*X/sigma

  diff_pmvnorm <- pmax(
    do.call(rbind,
            lapply(which(!idx_0),
                   llk_observed_i)
    ), .Machine$double.eps)

  grad_llk[!idx_0,] <- p2[!idx_0,]/as.numeric(diff_pmvnorm)


  return(grad_llk)
}

dllkdrho <- function(y, beta, X, gamma, Z,
                      thresholds, rho, sigma){

  dff_pnorm <- function(i){
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
  dff_dnorm <- function(i){
    p1 <- mvtnorm::dmvnorm(
      x = c(
        ((thresholds[m_index + 1] - xhat)/sigma)[i],
        zhat[i]
      ),
      sigma = rho_matrix)
    p2 <- mvtnorm::dmvnorm(
       c(
        ((thresholds[m_index] - xhat)/sigma)[i],
        zhat[i]
      ),
      sigma = rho_matrix)

    return(p1 - p2)
  }

  # tanh rho is estimated
  rho <- tanh(rho)

  # exp(sigma) estimated to get positive sigma
  sigma <- exp(sigma)

  zhat <- drop(Z %*% gamma)
  rho_matrix <- matrix(c(1, -rho, -rho, 1), nrow = 2)
  xhat <- drop(X %*% beta)
  m_index <- findInterval(y, thresholds)
  idx_0 <- (y == 0)


  # Gradient for observations where y == 0
  # stays 0 !

  diff_pmvnorm <- pmax(
    do.call(rbind,
            lapply(which(!idx_0),
                   dff_pnorm)
    ), .Machine$double.eps)
  diff_dmvnorm <- -do.call(rbind,
            lapply(which(!idx_0),
                   dff_dnorm)
    )*(1-rho^2)


  # because we estimate arctanh(rho)
  grad_llk <- as.numeric(diff_dmvnorm)/as.numeric(diff_pmvnorm)


  return(grad_llk)
}

inverse_mills_ratio <- function(x){
  dnorm(x)/pnorm(x)
}


wrap_pnorm <- function(x1,x2,rho){
  return(
    pnorm((x1 + rho*x2)/sqrt(1 - rho^2))
  )
}

