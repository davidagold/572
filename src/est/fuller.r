# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#########################################################################
# Dependencies

#########################################################################
# Fuller IV estimator

fuller <- function(y, X, Z, W = NULL, C = 1, ...) {
  n <- nrow(Z); pz <- ncol(Z)
  # X <- cbind(ones(n), x); Xt <- t(X); XtX <- Xt %*% X
    # I follow the development in [Hansen, Hausman, Newey 08, Section 2]
    # of the Fuller and related IV estimators

  # W are exogenous control variables. If empty, fill with an intercept
  if ( is.null(W) ) { W <- ones(n) }
  Wt <- t(W)
  
  # Needed for generality
  X <- as.matrix(X)
  Xt <- t(X)
  # Partial out controls for y, x, Z
  y <- y - P(W, y, Wt)
  X <- X - P(W, X, Wt)
  Z <- Z - P(W, Z, Wt)
  Zt <- t(Z)
  P.Z <- Z %*% solve(Zt %*% Z) %*% Zt
  
  O = cbind(y, X); Ot <- t(O)
  a.tilde <- solve(Ot %*% O, Ot %*% P(Z, O)) %>%
    eigen %>% { .$values } %>% min
  # a.hat <- ( a.tilde - (1 - a.tilde) * C/n ) / 
  #   ( 1 - (1 - a.tilde) * C/n )
  a.hat <- ( a.tilde - (1 - a.tilde) * C/(n-1-pz) ) / 
      ( 1 - (1 - a.tilde) * C/(n-1-pz) )
    
  XtPX <- Xt %*% P(Z, X)
  aXtX <- a.hat * Xt %*% X
  XtPy <- Xt %*% P(Z, y)
  aXty <- a.hat * Xt %*% y
  delta.hat <- solve(XtPX - aXtX, XtPy - aXty) # %>% as.numeric

  u.hat <- (y - X %*% delta.hat) %>% as.vector
  sigma2.hat <- as.numeric(t(u.hat) %*% u.hat) / (n-ncol(X))
  alpha.tilde <- as.numeric(t(u.hat) %*% P(Z, u.hat)) / as.numeric(t(u.hat) %*% u.hat)
  Upsilon.hat <- P(Z, X)
  X.tilde <- X - u.hat %*% (t(u.hat) %*% X) / as.numeric(t(u.hat) %*% u.hat)
  V.hat <- X.tilde - P(Z, X.tilde)
  kappa <- diag(P.Z)^2 %>% sum/n
  tau <- pz/n
  
  H.hat <- Xt %*% P(Z, X) - alpha.tilde * Xt %*% X
  SigmaB.hat <- sigma2.hat * (1 - alpha.tilde)^2 * t(X.tilde) %*% P(Z, X.tilde) +
    alpha.tilde^2 * t(X.tilde) %*% (X.tilde - P(Z, X.tilde))
  A.hat <- map(1:n, ~ { (P.Z[.,.] - tau) * Upsilon.hat[.,] %*% 
                        (map(1:n, ~ t(u.hat[.]^2 * V.hat[.,] / n)) 
                                    %>% reduce(sum)) 
                      }) %>% reduce(sum)
  B.hat <- map(1:n, ~ { (u.hat[.]^2 - sigma2.hat) * V.hat[.,] %*% t(V.hat[.,]) }) %>%
    reduce(sum) * pz * (kappa - tau) / (n*(1-2*tau+kappa*tau))
  
  Sigma.hat <- SigmaB.hat + A.hat + t(A.hat) + B.hat
  
  Var.hat <- solve(H.hat, Sigma.hat %*% solve(H.hat))
  # SE_delta.hat <- Var %>% diag %>% sqrt
  # SE_theta.hat <- SE_delta.hat[2]
  
  # In our case, we only have 1 endogenous x
  theta.hat <- delta.hat %>% as.numeric
  SE_theta.hat <- Var.hat %>% as.numeric %>% sqrt
  
  sigma.hat <- sqrt(sigma2.hat)
  
  list(theta.hat = theta.hat, sigma.hat = sigma.hat, 
       SE_theta.hat = SE_theta.hat, ...)
}
