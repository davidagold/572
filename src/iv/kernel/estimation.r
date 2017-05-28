# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#########################################################################
# Dependencies

suppressWarnings(library(dplyr))
suppressWarnings(library(tidyr))
suppressWarnings(library(purrr))
suppressWarnings(library(glmnet))

#########################################################################
# Estimation

beta0hat <- function(x, Z) {
  # fit <- cv.glmnet(Z, x)
}

theta0_tsls. <- function(y, x, Z) {
  ZtZ <- t(Z) %*% Z
  # P_Z <- Z %*% solve(t(Z) %*% Z) %*% t(Z)
  # theta0_hat <- (t(x) %*% P_Z %*% y) / (t(x) %*% P_Z %*% x) %>%
  #   as.numeric
  
  num <- (t(x) %*% Z %*% solve(ZtZ, y)) %>% as.numeric
  denom <- (t(x) %*% Z %*% solve(ZtZ, x)) %>% as.numeric
  theta0_hat <- num / denom

  sigma0_hhat <- sd(y - x %*% theta0_hat)
  var_theta0_hat <- sigma0_hhat^2 / (t(x) %*% P_Z %*% x)
  
  list(theta0_hat = theta0_hat, sigma0_hhat = sigma0_hhat, 
       var_theta0_hat = var_theta0_hat)
}
