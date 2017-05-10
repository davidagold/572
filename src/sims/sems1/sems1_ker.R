# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(dplyr)
library(tidyr)
library(MASS)
library(mvtnorm)
library(glmnet)
library(flare)

# Approximation of Lambda defined in eq (3.10), p8
Lambdahat <- function(X, t = .05, m = 100){
  n <- nrow(X)
  g <- rnorm(n * m, 0, 1)
  (do.call(rbind, replicate(100, X, simplify = FALSE)) * g) %>% 
    { data.frame(i = rep(1:n, times=m), r = rep(1:m, each = n), .)} %>%
    gather(X_j, X_ij, contains("X")) %>%
    group_by(r, X_j) %>% 
    summarize(sum_X_ij = sum(X_ij)) %>%
    summarize(maxsum_X_ij = max(sum_X_ij)) %>%
    summarize(q = quantile(maxsum_X_ij, 1-t)) %>% 
    as.numeric
}

.lambda <- function(sigma0hat, t = .05, p, c = 1.1) { 2 * c * sigma0hat * qnorm(1-(t/(2*p))) }

# Compute sigma0hat using iterative procedure described in Algorithm 1, p35
.sigma0hat <- function(Y, X, t = .05, c = 1.1, psi = .01, K = 100, nu = 1e-2){
  p <- ncol(X) + 1
  sigma0hat_k0 <- 0
  sigma0hat_k1 <- psi * sd(Y)
  sigma0hats <- numeric(K)
  k <- 1
  while ( abs(sigma0hat_k1 - sigma0hat_k0) > nu & k < K ) {
    # lambda <- 2 * c * sigma0hat_k0 * Lambdahat(X, t = t)
    # lambda <- 2 * c * sigma0hat_k1 * qnorm(1-(t/(2*p)))
    lambda <- .lambda(sigma0hat_k1, t, p)
    fit <- glmnet(X, Y)
    sigma0hat_k0 <- sigma0hat_k1
    sigma0hat_k1 <- (predict(fit, X, s = lambda) - Y) %>% sd
    sigma0hats[k] <- sigma0hat_k1
    k <- k + 1
  }
  sigma0hat_k1
}

.Y <- function(n, X, beta, gs){
  cbind(rep(1, n), X) %*% beta + rnorm(n, 0, gs)
}

.X <- function(n, p, Sigma){
  Z <- rmvnorm(n, rep(0, p-1), Sigma)
}

sim <- function(t){
  n <- 100
  rho0 <- .5
  s <- 6
  p <- 500
  Sigma <- rho0^abs(outer(1:(p-1), 1:(p-1), "-"))
  beta0_S <- c(1, 1, 1/2, 1/3, 1/4, 1/5)
  s <- length(beta0_S)
  beta0 <- c(beta0_S, rep(0, p-s))
  # sigma0s <- c(1, .1)
  sigma0 <- .1
  n_estimators <- 2
  
  R <- p * n_estimators
  
  trial <- numeric(R)
  sigma0_ <- numeric(R)
  estimator <- numeric(R)
  j <- numeric(R)
  betahat_j <- numeric(R)
  beta0_j <- numeric(R)
  
  # Generate data
  X <- .X(n, p, Sigma)
  Y <- .Y(n, X, beta0, sigma0)
  
  # Lasso
  r <- 0
  lasso_fit <- glmnet(X, Y)
  sigma0hat <- .sigma0hat(Y, X)
  lambda <- .lambda(sigma0hat, t, p)
  
  trial[r+1:(r+p)] <- t
  sigma0_[r+1:(r+p)] <- sigma0
  estimator[r+1:(r+p)] <- "Lasso"
  j[r+1:(r+p)] <- 1:p
  betahat_j[r+1:(r+p)] <- predict(lasso_fit, type = "coefficients", s = lambda) %>% as.numeric
  beta0_j[r+1:(r+p)] <- beta0
    
  # CV Lasso
  r <- r + p
  cv_lasso_fit <- cv.glmnet(X, Y, alpha=1, nfolds=5)
  
  trial[r+1:(r+p)] <- t
  sigma0_[r+1:(r+p)] <- sigma0
  estimator[r+1:(r+p)] <- "CV_Lasso"
  j[r+1:(r+p)] <- 1:p
  betahat_j[r+1:(r+p)] <- as.numeric(coef(cv_lasso_fit, s = "lambda.min"))
  beta0_j[r+1:(r+p)] <- beta0
  
  data.frame(
    trial = trial,
    sigma0 = sigma0_,
    estimator = estimator,
    j = j,
    betahat_j = betahat_j,
    beta0_j = beta0_j
  )
}
