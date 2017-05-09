# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# library(dplyr)
library(MASS)
library(mvtnorm)
library(glmnet)

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
  n_estimators <- 1
  
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
    
  # CV Lasso
  r <- 0
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
