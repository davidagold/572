# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#########################################################################
# Dependencies

#########################################################################
# Fuller IV estimator

fuller <- function(y, x, Z, C = 1, ...) {
  n <- nrow(Z)
  Z <- cbind(ones(n), Z)
  X <- cbind(ones(n), x); Xt <- t(X); XtX <- Xt %*% X
    # I follow the development in [Hansen, Hausman, Newey 08, Section 2]
    # of the Fuller and related IV estimators
  O = cbind(y, X); Ot <- t(O)
  a.tilde <- solve(Ot %*% O, Ot %*% P(Z, O)) %>%
    eigen %>% { .$values } %>% min
  a.hat <- ( a.tilde - (1 - a.tilde) * C/n ) / 
    ( 1 - (1 - a.tilde) * C/n )
  
  # xtPx <- t(x) %*% P(S.op, x)
  # axtx <- a.hat * t(x) %*% x
  # xtPy <- t(x) %*% P(S.op, y)
  # axty <- a.hat * t(x) %*% y
  # theta.hat <- solve(xtPx - axtx, xtPy - axty) %>%
  # as.numeric
  # delta.hat <- theta.hat
  XtPX <- Xt %*% P(Z, X)
  aXtX <- a.hat * Xt %*% X
  XtPy <- Xt %*% P(Z, y)
  aXty <- a.hat * Xt %*% y
  
  delta.hat <- solve(XtPX - aXtX, XtPy - aXty) %>%
    as.numeric
  theta.hat <- delta.hat[2]
  # print(delta.hat)
  
  u.hat <- y - X %*% delta.hat
  sigma.hat2 <- as.numeric(t(u.hat) %*% u.hat) / (n-ncol(X))
  alpha.tilde <- as.numeric(t(u.hat) %*% P(Z, u.hat)) / as.numeric(t(u.hat) %*% u.hat)
  X.tilde <- X - u.hat %*% (t(u.hat) %*% X) / as.numeric(t(u.hat) %*% u.hat)
  
  H.hat <- Xt %*% P(Z, X) - alpha.tilde * Xt %*% X
  Sigma.hat <- sigma.hat2 * (1 - alpha.tilde)^2 * t(X.tilde) %*% P(Z, X.tilde) +
    alpha.tilde^2 * t(X.tilde) %*% (X.tilde - P(Z, X.tilde))
  
  Var <- solve(H.hat, Sigma.hat %*% solve(H.hat))
  SE_delta.hat <- Var %>% diag %>% sqrt
  SE_theta.hat <- SE_delta.hat[2]
  sigma.hat <- sqrt(sigma.hat2)
  
  list(theta.hat = theta.hat, sigma.hat = sigma.hat, 
       SE_theta.hat = SE_theta.hat, ...)
}
