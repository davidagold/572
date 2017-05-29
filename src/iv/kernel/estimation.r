# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#########################################################################
# Dependencies

suppressWarnings(library(dplyr))
suppressWarnings(library(tidyr))
suppressWarnings(library(purrr))
suppressWarnings(library(glmnet))

#########################################################################
# Estimation

.Lambda.hat <- function(Z, t = .05, m = 500){
  n <- nrow(Z)
  g <- rnorm(n * m, 0, 1)
  map(1:m, ~ (Z * rnorm(n, 0, 1)) %>% 
        apply(2, sum) %>% 
        max) %>%
    as.numeric %>%
    quantile(1 - t)
}

.lambda <- function(sigma0.hat, Z, t = 0.05, c = 1.1) {  2 * c * sigma0.hat * .Lambda.hat(Z, t = t) }

sigma0_v.hat_iter <- function(x, Z, t = .05, c = 1.1, psi = .1, K = 100){
  n <- nrow(Z)
  p <- ncol(Z)
  # sigma0v.hat_k0 <- psi * sd(Y)
  Lambda.hat <- .Lambda.hat(Z, t = t)
  sigma0v.hat_k0 <- Inf
  sigma0v.hat_k1 <- psi * sd(x)
  nu <- .2 * sd(x)
  # sigma0hats <- numeric(K)
  k <- 1
  fit <- glmnet(Z, x)
  while ( abs(sigma0v.hat_k1 - sigma0v.hat_k0) > nu & k < K ) {
    sigma0v.hat_k0 <- sigma0v.hat_k1
    lambda <- 2 * c * sigma0v.hat_k1 * Lambda.hat
    (predict(fit, Z, s = lambda/n, exact = TRUE, x = Z, y = x) - x)^2/n
    
    sigma0v.hat_k1 <- sum((predict(fit, Z, s = lambda/n, exact = TRUE, x = Z, y = x) - x)^2/n) %>% sqrt
    k <- k + 1
  }
  sigma0v.hat_k0
}

beta0hat <- function(x, Z) {
  # fit <- cv.glmnet(Z, x)
}

fit.IV <- function(y, x, Z) {
  n <- length(y)
  tZ <- t(Z)
  ZtZ <- tZ %*% Z
  tx <- t(x)
  
  num <- (tx %*% Z %*% solve(ZtZ, tZ %*% y)) %>% as.numeric
  denom <- (tx %*% Z %*% solve(ZtZ, tZ %*% x)) %>% as.numeric
  theta0.hat <- num / denom

  sigma0h.hat <- (sum((y - x %*% theta0.hat)^2)/(n-1)) %>% sqrt
  SE_theta0.hat <- (sigma0h.hat / sqrt(denom))
  
  list(theta0.hat = theta0.hat, sigma0h.hat = sigma0h.hat, 
       SE_theta0.hat = SE_theta0.hat)
}
