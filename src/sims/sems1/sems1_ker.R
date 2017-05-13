# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(dplyr)
library(tidyr)
library(purrr)
library(MASS)
library(mvtnorm)
library(glmnet)
library(flare)

# Approximation of Lambda defined in eq (3.10), p8
Lambdahat <- function(X, t = .05, m = 100){
  n <- nrow(X)
  g <- rnorm(n * m, 0, 1)
  map(1:m, ~ (X * rnorm(m, 0, 1)) %>% 
        apply(2, sum) %>% 
        max) %>%
    as.numeric %>%
    quantile(1 - t)
}

# .lambda <- function(sigma0hat, n, p, t = .05, c = 1.1) { 2 * sqrt(n) * c * sigma0hat * qnorm(1-(t/(2*p))) }
.lambda <- function(sigma0hat, X, t = 0.05, c = 1.1) {  2 * c * sigma0hat * Lambdahat(X, t = t) }

# Compute sigma0hat using iterative procedure described in Algorithm 1, p35
.sigma0hat <- function(Y, X, t = .05, c = 1.1, psi = .1, K = 100){
  n <- nrow(X)
  p <- ncol(X) + 1
  sigma0hat_k0 <- -Inf
  sigma0hat_k1 <- psi * sd(Y)
  nu <- sigma0hat_k1
  # sigma0hats <- numeric(K)
  k <- 1
  while ( abs(sigma0hat_k1 - sigma0hat_k0) > nu & k < K ) {
    lambda <- .lambda(sigma0hat_k1, X)
    # lambda <- .lambda(sigma0hat_k1, n, p, t = t)
    fit <- glmnet(X, Y)
    sigma0hat_k0 <- sigma0hat_k1
    sigma0hat_k1 <- (predict(fit, X, s = lambda/n, exact = TRUE, x = X, y = Y) - Y) %>% sd
    # sigma0hats[k] <- sigma0hat_k1
    k <- k + 1
    # print(sigma0hat_k1)
  }
  sigma0hat_k1
}

.Y <- function(n, X, beta0, gs){
  cbind(rep(1, n), X) %*% beta0 + rnorm(n, 0, gs)
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
  sigma0s <- c(1, .1)
  n_sigmas <- length(sigma0s)
  n_estimators <- 5
  
  R <- p * n_estimators * n_sigmas
  
  trial <- numeric(R)
  sigma0_ <- numeric(R)
  estimator <- numeric(R)
  j <- numeric(R)
  betahat_j <- numeric(R)
  beta0_j <- numeric(R)
  
  for ( i in 1:2 ) {
    # Generate data
    sigma0 <- sigma0s[i]
    X <- .X(n, p, Sigma)
    Y <- .Y(n, X, beta0, sigma0)
    
    # Lasso
    r <- p * n_estimators * (i - 1)
    
    lasso_fit <- glmnet(X, Y)
    # sigma0hat <- .sigma0hat(Y, X)
    sigma0hat <- sigma0
    # lambda <- .lambda(sigma0hat, n, p, t = t)
    lambda <- .lambda(sigma0hat, X)
    L_betahat <- predict(lasso_fit, type = "coefficients", 
                         s = lambda/n, exact = TRUE, x = X, y = Y) %>% 
      as.numeric
    
    trial[(r+1):(r+p)] <- t
    sigma0_[(r+1):(r+p)] <- sigma0
    estimator[(r+1):(r+p)] <- "Lasso"
    j[(r+1):(r+p)] <- 1:p
    betahat_j[(r+1):(r+p)] <- L_betahat
    beta0_j[(r+1):(r+p)] <- beta0

    # Post-Lasso
    r <- r + p
    Shat <- which(L_betahat != 0)
    shat <- length(Shat)
    if ( 1 %in% Shat ) {
      if ( shat > 1 ) {
        PL_fit <- lm(Y ~ 1 + X[, Shat[2:shat]-1])
      } else {
        PL_fit <- lm(Y ~ 1)
      }
    } else {
      PL_fit <- lm(Y ~ -1 + X[, Shat-1])
    }
    PL_betahat <- numeric(p)
    PL_betahat[Shat] <- coefficients(PL_fit)

    trial[(r+1):(r+p)] <- t
    sigma0_[(r+1):(r+p)] <- sigma0
    estimator[(r+1):(r+p)] <- "Post-Lasso"
    j[(r+1):(r+p)] <- 1:p
    betahat_j[(r+1):(r+p)] <- PL_betahat
    beta0_j[(r+1):(r+p)] <- beta0

    # CV Lasso
    r <- r + p
    cv_lasso_fit <- cv.glmnet(X, Y, alpha=1, nfolds=5)

    trial[(r+1):(r+p)] <- t
    sigma0_[(r+1):(r+p)] <- sigma0
    estimator[(r+1):(r+p)] <- "CV_Lasso"
    j[(r+1):(r+p)] <- 1:p
    betahat_j[(r+1):(r+p)] <- as.numeric(coef(cv_lasso_fit, s = "lambda.min"))
    beta0_j[(r+1):(r+p)] <- beta0
   
    # Iterated Lasso
    r <- r + p
    
    sigma0hat <- .sigma0hat(Y, X)
    lambda <- .lambda(sigma0hat, X)
    IL_betahat <- predict(lasso_fit, type = "coefficients", 
                         s = lambda/n, exact = TRUE, x = X, y = Y) %>% 
      as.numeric
    
    trial[(r+1):(r+p)] <- t
    sigma0_[(r+1):(r+p)] <- sigma0
    estimator[(r+1):(r+p)] <- "Iter-Lasso"
    j[(r+1):(r+p)] <- 1:p
    betahat_j[(r+1):(r+p)] <- IL_betahat
    beta0_j[(r+1):(r+p)] <- beta0
    
    # Post-Iterated-Lasso
    r <- r + p
    
    IL_Shat <- which(IL_betahat != 0)
    IL_shat <- length(IL_Shat)
    if ( 1 %in% IL_Shat ) {
      if ( IL_shat > 1 ) {
        PIL_fit <- lm(Y ~ 1 + X[, IL_Shat[2:IL_shat]-1])
      } else {
        PIL_fit <- lm(Y ~ 1)
      }
    } else {
      PIL_fit <- lm(Y ~ -1 + X[, IL_Shat-1])
    }
    PIL_betahat <- numeric(p)
    PIL_betahat[IL_Shat] <- coefficients(PIL_fit)
    # print(Shat)
    # print(PL_betahat)
    
    trial[(r+1):(r+p)] <- t
    sigma0_[(r+1):(r+p)] <- sigma0
    estimator[(r+1):(r+p)] <- "Post-Iter-Lasso"
    j[(r+1):(r+p)] <- 1:p
    betahat_j[(r+1):(r+p)] <- PIL_betahat
    beta0_j[(r+1):(r+p)] <- beta0
  }

  data.frame(
    trial = trial,
    sigma0 = sigma0_,
    estimator = estimator,
    j = j,
    betahat_j = betahat_j,
    beta0_j = beta0_j
  )
}

# d <- sim(1)
# d %>%
#   filter(is.na(betahat_j))
