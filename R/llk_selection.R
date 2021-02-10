llk_selection <- function(y, y_selection, beta, X, gamma, Z,
                          thresholds, rho, sigma){

  # tanh rho is estimated
  rho <- tanh(rho)

  # exp(sigma) estimated to get positive sigma
  sigma <- exp(sigma)
  zhat <- drop(Z %*% gamma)
  rho_matrix <- matrix(c(1, -rho, -rho, 1), nrow = 2)
  xhat <- drop(X %*% beta)
  idx_0 <- (y_selection == 0)

  if (is.finite(max(thresholds))) thresholds <- c(thresholds, Inf)
  if (is.finite(min(thresholds))) thresholds <- c(-Inf, thresholds)
  # outcome_index <- findInterval(y, thresholds)
  outcome_index <- y


  llk <- rep(NA, nrow(X))
  llk[idx_0] <- pnorm(-zhat, log.p = TRUE)[idx_0]

  llk[!idx_0] <- log(
    pmax(
      do.call(rbind,
              lapply(which(!idx_0),
                     dff_pnorm,
                     thresholds = thresholds,
                     outcome_index = outcome_index,
                     xhat = xhat, zhat = zhat, sigma = sigma,
                     rho_matrix = rho_matrix)
      ), .Machine$double.eps)
  )

  return(llk)
}


grad_llk_selection_wrapper <- function(theta, y, y_selection,
                                       X, Z, thresholds){
  beta <- theta[1:ncol(X)]
  gamma <- theta[seq(from = ncol(X) + 1,
                     length.out = ncol(Z))]
  rho <- theta[ncol(X) + ncol(Z) + 1]
  sigma <- theta[ncol(X) + ncol(Z) + 2]
  return(
    grad_llk_selection(y, y_selection, beta, X, gamma, Z,
                       thresholds, rho, sigma)
  )
}

llk_selection_wrapper <- function(theta, y, y_selection,
                                  X, Z, thresholds){
  beta <- theta[1:ncol(X)]
  gamma <- theta[seq(from = ncol(X) + 1,
                     length.out = ncol(Z))]
  rho <- theta[ncol(X) + ncol(Z) + 1]
  sigma <- theta[ncol(X) + ncol(Z) + 2]
  return(
    llk_selection(y, y_selection, beta, X, gamma, Z,
                  thresholds, rho, sigma)
  )
}


grad_llk_selection <- function(y, y_selection,
                               beta, X, gamma, Z,
                               thresholds, rho, sigma){

  # tanh rho is estimated
  rho <- tanh(rho)

  # exp(sigma) estimated to get positive sigma
  sigma <- exp(sigma)

  zhat <- drop(Z %*% gamma)
  rho_matrix <- matrix(c(1, -rho, -rho, 1), nrow = 2)
  xhat <- drop(X %*% beta)
  # outcome_index <- findInterval(y, thresholds)
  outcome_index <- y
  idx_0 <- (y_selection == 0)


  # Vector of parameters:  (beta, gamma, rho, sigma)
  grad <- matrix(0, length(y), length(beta) + length(gamma) + 2)

  # Plug gradient wrt beta
  grad[,1:length(beta)] <- dllkdbeta(y, y_selection,
                                     beta, X, gamma, Z,
                                     thresholds, rho, sigma)

  # Plug gradient wrt gamma
  grad[,seq(from = length(beta)+1,
            length.out = length(gamma))] <- dllkdgamma(y, y_selection,
                                                       beta, X, gamma, Z,
                                                       thresholds, rho, sigma)

  # Plug gradient wrt rho
  grad[y_selection!=0,length(beta) +length(gamma) + 1] <- dllkdrho(y, y_selection,
                                                         beta, X, gamma, Z,
                                                         thresholds, rho, sigma)

  # Plug gradient wrt sigma
  grad[y_selection!=0,length(beta) +length(gamma) + 2] <- dllkdsigma(y, y_selection,
                                                           beta, X, gamma, Z,
                                                           thresholds, rho, sigma)

  return(grad)
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


dllkdgamma <- function(outcome_index, idx_0,
                       xhat, zhat,
                       thresholds, rho, sigma,
                       rho_matrix){


  grad_llk <- matrix(0, nrow(Z), length(gamma))

  # Gradient for observations where y == 0
  grad_llk[idx_0,] <- - (inverse_mills_ratio(- zhat) * Z)[idx_0,]

  p2 <- wrap_pnorm((thresholds[outcome_index + 1] - xhat)/sigma,
                   zhat,
                   rho)
  p2 <- p2 - wrap_pnorm((thresholds[outcome_index] - xhat)/sigma,
                        zhat,
                        rho)
  p2 <- p2*dnorm(zhat)*Z

  diff_pmvnorm <- pmax(
    do.call(rbind,
            lapply(which(!idx_0),
                   dff_pnorm,
                   thresholds = thresholds,
                   outcome_index = outcome_index,
                   xhat = xhat, zhat = zhat, sigma = sigma,
                   rho_matrix = rho_matrix)
    ), .Machine$double.eps)

  grad_llk[!idx_0,] <- p2[!idx_0,]/as.numeric(diff_pmvnorm)


  return(grad_llk)
}


dllkdbeta <- function(outcome_index, idx_0,
                      xhat, zhat,
                      thresholds, rho, sigma,
                      rho_matrix){

  grad_llk <- matrix(0, nrow(X), length(beta))

  # Gradient for observations where y == 0
  # stays 0 !

  p2_pnorm1 <- wrap_pnorm(zhat,
                          (thresholds[outcome_index + 1] - xhat)/sigma,
                          rho)
  p2_dnorm1 <- dnorm((thresholds[outcome_index + 1] - xhat)/sigma)
  p2_pnorm2 <- wrap_pnorm(zhat,
                          (thresholds[outcome_index] - xhat)/sigma,
                          rho)
  p2_dnorm2 <- dnorm((thresholds[outcome_index] - xhat)/sigma)
  p2 <- p2_pnorm1*p2_dnorm1 -  p2_pnorm2*p2_dnorm2
  p2 <- -p2*X/sigma

  diff_pmvnorm <- pmax(
    do.call(rbind,
            lapply(which(!idx_0),
                   dff_pnorm,
                   thresholds = thresholds,
                   outcome_index = outcome_index,
                   xhat = xhat, zhat = zhat, sigma = sigma,
                   rho_matrix = rho_matrix)
    ), .Machine$double.eps)

  grad_llk[!idx_0,] <- p2[!idx_0,]/as.numeric(diff_pmvnorm)


  return(grad_llk)
}

dllkdrho <- function(outcome_index, idx_0,
                     xhat, zhat,
                     thresholds, rho, sigma,
                     rho_matrix){

  # Gradient for observations where y == 0
  # stays 0 !

  diff_pmvnorm <- pmax(
    do.call(rbind,
            lapply(which(!idx_0),
                   dff_pnorm,
                   thresholds = thresholds,
                   outcome_index = outcome_index,
                   xhat = xhat, zhat = zhat, sigma = sigma,
                   rho_matrix = rho_matrix)
    ), .Machine$double.eps)
  diff_dmvnorm <- -do.call(rbind,
                           lapply(which(!idx_0),
                                  dff_dnorm)
  )*(1-rho^2)


  # because we estimate arctanh(rho)
  grad_llk <- as.numeric(diff_dmvnorm)/as.numeric(diff_pmvnorm)


  return(grad_llk)
}

dllkdsigma <- function(outcome_index, idx_0,
                       xhat, zhat,
                       thresholds, rho, sigma,
                       rho_matrix){


  diff_numerator <- function(i){
    p1 <- wrap_dpnorm(zhat[i], ((thresholds[outcome_index + 1] - xhat)/sigma)[i],
                      rho, sigma)
    p2 <- wrap_dpnorm(zhat[i], ((thresholds[outcome_index] - xhat)/sigma)[i],
                      rho, sigma)
    return(p1 - p2)
  }


  diff_pmvnorm <- pmax(
    do.call(rbind,
            lapply(which(!idx_0),
                   dff_pnorm,
                   thresholds = thresholds,
                   outcome_index = outcome_index,
                   xhat = xhat, zhat = zhat, sigma = sigma,
                   rho_matrix = rho_matrix)
    ), .Machine$double.eps)


  diff_dpnorm <- lapply(which(!idx_0),
                        diff_numerator)


  grad_llk <- as.numeric(diff_dpnorm)/as.numeric(diff_pmvnorm)

  # we need gradient with respect to log(sigma)
  grad_llk <- grad_llk*sigma

  return(grad_llk)
}

dff_pnorm <- function(i, thresholds, outcome_index,
                      xhat, zhat, sigma,
                      rho_matrix){

  p1 <- mvtnorm::pmvnorm(
    upper = c(
      ((thresholds[outcome_index + 1] - xhat)/sigma)[i],
      zhat[i]
    ),
    sigma = rho_matrix)
  p2 <- mvtnorm::pmvnorm(
    upper = c(
      ((thresholds[outcome_index] - xhat)/sigma)[i],
      zhat[i]
    ),
    sigma = rho_matrix)

  return(p1 - p2)

}

dff_dnorm <- function(i, thresholds, outcome_index,
                      xhat, zhat, sigma,
                      rho_matrix){

  p1 <- mvtnorm::dmvnorm(
    x = c(
      ((thresholds[outcome_index + 1] - xhat)/sigma)[i],
      zhat[i]
    ),
    sigma = rho_matrix)
  p2 <- mvtnorm::dmvnorm(
    c(
      ((thresholds[outcome_index] - xhat)/sigma)[i],
      zhat[i]
    ),
    sigma = rho_matrix)

  return(p1 - p2)
}




inverse_mills_ratio <- function(x){
  dnorm(x)/pnorm(x)
}


wrap_pnorm <- function(x1,x2,rho){
  return(
    pnorm((x1 + rho*x2)/sqrt(1 - rho^2))
  )
}


wrap_dpnorm <- function(x1,x2,rho,sigma){
  p <- -wrap_pnorm(x1,x2,rho)*dnorm(x2)*x2/sigma
  return(
    ifelse(is.infinite(x2), 0, p)
  )
}

