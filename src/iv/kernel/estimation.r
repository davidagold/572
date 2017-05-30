# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#########################################################################
# Dependencies
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

.lambda <- function(sigma.hat, Z, t = 0.05, c = 1.1) {  2 * c * sigma.hat * .Lambda.hat(Z, t = t) }

lambda.IL <- function(obs) {
  with(obs, {
    tau.hat <- tau.hat_iter(x, Z)
    lambda <- .lambda(sigma0.hat = tau.hat, Z = Z)
  })
  lambda
}

lambda.CV <- function(obs) "lambda.min"

tau.hat_iter <- function(obs, t = .05, c = 1.1, psi = .1, K = 100){
  y <- obs$y; x <- obs$x; Z <- obs$Z
  Lambda.hat <- .Lambda.hat(Z, t = t)
  tau.hat_k0 <- Inf
  tau.hat_k1 <- psi * sd(x)
  nu <- .2 * sd(x)
  k <- 1
  fit <- glmnet(Z, x)
  while ( abs(tau.hat_k1 - tau.hat_k0) > nu & k < K ) {
    tau.hat_k0 <- tau.hat_k1
    lambda <- 2 * c * tau.hat_k1 * Lambda.hat
    (predict(fit, Z, s = lambda/n, exact = TRUE, x = Z, y = x) - x)^2/n
    
    tau.hat_k1 <- sum((predict(fit, Z, s = lambda/n, exact = TRUE, x = Z, y = x) - x)^2/n) %>% sqrt
    k <- k + 1
  }
  tau.hat_k0
}

fit.Fuller <- function(obs, S.op, C = 1, ...) {
  with(obs, {
    # print(S.op)
    # I follow the development in [Hansen, Hausman, Newey 08, Section 2]
    # of the Fuller and related IV estimators
    print(C)
    O = cbind(y, X); Ot <- t(O)
    a.tilde <- solve(Ot %*% O, Ot %*% P(S.op, O)) %>%
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
    XtPX <- Xt %*% P(S.op, X)
    aXtX <- a.hat * Xt %*% X
    XtPy <- Xt %*% P(S.op, y)
    aXty <- a.hat * Xt %*% y

    delta.hat <- solve(XtPX - aXtX, XtPy - aXty) %>%
      as.numeric
    theta.hat <- delta.hat[2]
    print(delta.hat)
    
    u.hat <- y - X %*% delta.hat
    sigma.hat2 <- as.numeric(t(u.hat) %*% u.hat) / (n-ncol(X))
    alpha.tilde <- as.numeric(t(u.hat) %*% P(S.op, u.hat)) / as.numeric(t(u.hat) %*% u.hat)
    X.tilde <- X - u.hat %*% (t(u.hat) %*% X) / as.numeric(t(u.hat) %*% u.hat)
    
    H.hat <- Xt %*% P(S.op, X) - alpha.tilde * Xt %*% X
    Sigma.hat <- sigma.hat2 * (1 - alpha.tilde)^2 * t(X.tilde) %*% P(S.op, X.tilde) +
      alpha.tilde^2 * t(X.tilde) %*% (X.tilde - P(S.op, X.tilde))
    
    Var <- solve(H.hat, Sigma.hat %*% solve(H.hat))
    SE_delta.hat <- Var %>% diag %>% sqrt
    SE_theta.hat <- SE_delta.hat[2]
    sigma.hat <- sqrt(sigma.hat2)
    
    
    
    
    # I follow the development in [Bekker 1994, p.666] of the
    # estimates of asymptotic covariance
    # dof <- n - 1
    # S.bar <- (Ot %*% P(S.op, O)) / dof
    # S.perp <- Ot %*% (O - P(S.op, O)) / dof
    # S <- S.bar + S.perp
    # S.bar.22 <- S.bar[2:3, 2:3]
    # # S.bar.22 <- S.bar[2, 2]
    # 
    # sigma.hat2 <- t(c(1, -delta.hat)) %*% S %*% c(1, -delta.hat) %>%
    #   as.numeric
    # sigma.hat <- sqrt(sigma.hat2)
    # # print(sigma.hat2)
    # 
    # B <- S.bar - a.hat * S.perp
    # B22 <- B[2:3, 2:3]
    # # B22 <- B[2, 2]
    # C <- S.bar - a.hat *
    #   (S.perp %*% (c(1, -delta.hat) %*% t(c(1, -delta.hat))) %*% S.perp ) /
    #   as.numeric(t(c(1, -delta.hat)) %*% S.perp %*% c(1, -delta.hat))
    # # C22 <- C[2:3, 2:3]
    # C22 <- C[2, 2]
    # D <- a.hat*(C-B)
    # D22 <- D[2:3, 2:3]
    # # D22 <- D[2, 2]
    # 
    # Ahat <- sigma.hat2 * solve(B22, (C22 + D22) %*% solve(B22))
    # # Ahat <- sigma.hat2 * (C22 + D22) / (B22^2)
    # SE_delta.hat <- diag(Ahat) %>% sqrt
    # SE_theta.hat <- SE_delta.hat[2]
    # # SE_theta.hat <- sqrt(Ahat)
    # # print(SE_delta.hat)
    
    list(theta.hat = theta.hat, sigma.hat = sigma.hat, 
         SE_theta.hat = SE_theta.hat, ...)
  })
}

fit.IV <- function(obs, S.op, ...) {
  with(obs, {
    num <- Xt %*% P(S.op, y) 
    denom <- Xt %*% P(S.op, X)
    delta.hat <- solve(denom, num) %>%
      as.numeric
    # theta.hat <- (num/denom) %>% as.numeric
    
    sigma.hat2 <- sum((y - X %*% delta.hat)^2) / (n-1)
    sigma.hat <- sqrt(sigma.hat2)
    SE_delta.hat <- (sigma.hat2 * solve(denom)) %>% diag %>% sqrt
    
    list(theta.hat = delta.hat[2], sigma.hat = sigma.hat, 
         SE_theta.hat = SE_delta.hat[2], ...)
  })
}
