# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#########################################################################
# Dependencies

#########################################################################
# lambda specification

tau.hat_iter <- function(x, Z, t = .05, c = 1.1, psi = .1, K = 100){
  n <- nrow(Z); pz <- ncol(Z)
  # Lambda.hat <- Lambda.hat_(Z, t = t)
  tau.hat_k0 <- Inf
  id.maxcor <- map_dbl(1:pz, ~cor(x, Z[,.])) %>% which.max
  tau.hat_k1 <- psi * mean((x - predict(lm(x ~ Z[,id.maxcor])))^2)
  sprintf("tau.hat_k1 = %.6f", tau.hat_k1) %>% print
  sprintf("sd(x) = %.6f", mean((x - mean(x))^2)) %>% print
  nu <- .2 * sd(x)
  k <- 1
  while ( abs(tau.hat_k1 - tau.hat_k0) > nu & k < K ) {
    tau.hat_k0 <- tau.hat_k1
    # lambda <- 2 * c * tau.hat_k1 * Lambda.hat / n
    lambda <- c * tau.hat_k1 * sqrt(2*log(2*pz)/n)
    fit <- lasso(x, Z, lambda = lambda)
    tau.hat_k1 <- (sum((x - fit$beta0.hat - Z %*% fit$beta)^2)/n) %>% sqrt
    sprintf("tau.hat_k1 = %.6f", tau.hat_k1) %>% print
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

lambda_ <- function(x, Z, c = 1.1, m = 500, gamma = .05, tau.est = "sd.x") {  
  n <- nrow(Z); pz <- ncol(Z)
  # if ( tau.est = "sd.x" ) {
  #   tau.hat <- 
  # }
  tau.hat <- tau.hat_iter(x, Z)
  lambda.dat <- c * tau.hat * Lambda.hat_(Z, m = m) / n
  lambda.thr <- c * tau.hat * sqrt(2*log(2*pz)) / sqrt(n)
  list(lambda.dat = lambda.dat, lambda.thr = lambda.thr)
}
