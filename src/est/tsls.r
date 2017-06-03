# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#########################################################################
# Dependencies

#########################################################################
# 2-Stage Least Squares

tsls <- function(y, X, Z, ...) {
  n <- nrow(Z)
  Z <- cbind(ones(n), Z)
  X <- cbind(ones(n), X); Xt <- t(X); XtX <- Xt %*% X
  num <- Xt %*% P(Z, y) 
  denom <- Xt %*% P(Z, X)
  delta.hat <- solve(denom, num) %>%
    as.numeric
  # theta.hat <- (num/denom) %>% as.numeric
  
  sigma.hat2 <- sum((y - X %*% delta.hat)^2) / (n-1)
  sigma.hat <- sqrt(sigma.hat2)
  SE_delta.hat <- (sigma.hat2 * solve(denom)) %>% diag %>% sqrt
  
  list(theta.hat = delta.hat[2], sigma.hat = sigma.hat, 
       SE_theta.hat = SE_delta.hat[2], ...)
}