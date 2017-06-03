# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#########################################################################
# Dependencies

#########################################################################
# lambda specification

tau.hat_iter <- function(x, Z, t = .05, c = 1.1, psi = .1, K = 100){
  n <- nrow(Z)
  Lambda.hat <- Lambda.hat_(Z, t = t)
  tau.hat_k0 <- Inf
  tau.hat_k1 <- psi * sd(x)
  nu <- .2 * sd(x)
  k <- 1
  while ( abs(tau.hat_k1 - tau.hat_k0) > nu & k < K ) {
    tau.hat_k0 <- tau.hat_k1
    lambda <- 2 * c * tau.hat_k1 * Lambda.hat / n
    fit <- lasso(x, Z, lambda = lambda)
    tau.hat_k1 <- (sum((x - fit$beta0.hat - Z %*% fit$beta)^2)/n) %>% sqrt
    print(tau.hat_k1)
    k <- k + 1
  }
  tau.hat_k1
}

Lambda.hat_ <- function(Z, t = .05, m = 500){
  n <- nrow(Z)
  g <- rnorm(n * m, 0, 1)
  map(1:m, ~ (Z * rnorm(n, 0, 1)) %>% 
        apply(2, function(x) { sum(x) %>% abs }) %>% 
        max) %>%
    as.numeric %>%
    quantile(1 - t)
}

lambda_ <- function(x, Z, c = 1.1, m = 500) {  
  2*c * tau.hat_iter(x, Z) * Lambda.hat_(Z, m = m) / nrow(Z) 
}
