# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#########################################################################
# Dependencies

#########################################################################
# Lasso estimator

lasso <- function(y, X, epsilon = NULL, K = 100, thresh = 1e-7, lambdas = NULL, lambda = NULL,
                  standardize = T, L = 200) {
  y <- as.vector(y)
  n <- nrow(X); p <- ncol(X)

  y.bar <- mean(y); beta0.hat <- y.bar
  y <- y - y.bar
  mu.X <- map_dbl(1:p, ~ { mean(X[,.]) })
  X <- map(1:p, ~ { X[,.] - mu.X[.] }) %>% reduce(cbind)
  sigma2.X <- map_dbl(1:p, ~ { mean(X[,.]^2) })

  if ( standardize == TRUE ) {
    # X <- scale(X, scale = sqrt(sigma2.X)) %>% as.matrix
    X <- map(1:p, ~ { X[,.] / sqrt(sigma2.X[.]) }) %>% reduce(cbind)
    sigma2.X <- rep(1, p)
  }

  epsilon <- ifelse(is.null(epsilon), ifelse(n > p, .0001, .01), epsilon)

  if ( is.null(lambdas) ) {
    # lambda.max <- max(abs(colSums(cbind(ones(n), X) * (y + y.bar))))/n
    lambda.max <- max(abs(colSums(X * y))) / n
    lambda.min <- epsilon * lambda.max
    lambdas <- exp(seq(log(lambda.max), log(lambda.min), length.out = K))
  }

  if ( !is.null(lambda) ) {
    lambda.id <- map(lambdas, ~abs(.-lambda)) %>% which.min %>% as.numeric
    lambdas[lambda.id] <- lambda
  } else {
    lambda.id <- 1 # not meant to be used
  }

  betas <- matrix(nrow = p, ncol = K)
  betas[,1] <- 0
  Ls <- map(1:p, ~ setdiff(1:p, .)) %>% reduce(cbind)

  Xy <- (t(X) %*% y) %>% as.vector

  XX <- t(X) %*% X

  # outer loop over lambda path
  for ( k in 2:K ) {
    lambda <- lambdas[k]
    beta <- betas[,k-1] %>% as.vector
    r <- 0
    # print(sprintf("k = %1.0f", k))
    for ( l in 1:L ) {
      # print(length(S))
      # inner loop over coordinates
      for ( j in sample(1:p) ) {
        xjr <- Xy[j] - sum(XX[j,] * beta)
        # print(xjr)
        Gj <- xjr/n + sigma2.X[j] * beta[j]
        # Gj <- xjr/n + sigma2.X[j] * beta[j]
        beta[j] <- beta.j <- st(Gj, lambda) / sigma2.X[j]
      }
    }
    betas[,k] <- beta
  }
  if ( is.null(lambda) ) { beta <- NULL } else { beta <- betas[, lambda.id] }
  list(lambdas = lambdas, betas = betas, beta0.hat = beta0.hat, beta = beta,
       lambda.id = lambda.id)
}

#########################################################################
# CV Lasso

cv.lasso <- function(y, X, nfolds = 10, epsilon = NULL, K = 100, thresh = 1e-7,
                     standardize = TRUE, ...) {
  y <- as.vector(y)
  n <- nrow(X); p <- ncol(X)
  nperfold <- n %/% nfolds
  is <- sample(1:n)
  Is <- map(0:(nfolds-1), ~ is[(.*nperfold+1):min(((.+1)*nperfold), n)])

  y.bar <- mean(y); beta0.hat <- y.bar
  y <- y - y.bar
  mu.X <- map_dbl(1:p, ~ { mean(X[,.]) })
  X <- map(1:p, ~ { X[,.] - mu.X[.] }) %>% reduce(cbind)
  sigma2.X <- map_dbl(1:p, ~ { mean(X[,.]^2) })

  # y <- scale(y, scale = sd.y) %>% as.vector
  if ( standardize == TRUE ) {
    # X <- scale(X, scale = sqrt(sigma2.X)) %>% as.matrix
    X <- map(1:p, ~ { X[,.] / sqrt(sigma2.X[.]) }) %>% reduce(cbind)
    sigma2.X <- rep(1, p)
  }
  # l2.X <- map_dbl(1:p, ~ { sum(X[,.]^2) %>% as.numeric })
  # print(l2.X)

  epsilon <- ifelse(is.null(epsilon), ifelse(n > p, .0001, .01), epsilon)

  lambda.max <- max(abs(colSums(X * y))) / n
  lambda.min <- epsilon * lambda.max
  lambdas <- exp(seq(log(lambda.max), log(lambda.min), length.out = K))

  fits <- map(1:nfolds, ~ { Ic <- setdiff(1:n, Is[.]); fit <- lasso(y[Ic], X[Ic,], lambdas = lambdas) ;
                      print(sprintf("%2.f%% done fitting", 100*./nfolds)); fit })
  # mses <- map(1:nfolds, ~ { I <- Is[[.]]; fit <- fits[[.]];
  #                           map_dbl(1:K, ~ mean((y[I] - fit$beta0.hat - my.scale(X[I,]) %*% fit$betas[,.])^2))
  #                         })
  mses <- map(1:nfolds, ~ { I <- Is[[.]]; fit <- fits[[.]];
                            map_dbl(1:K, ~ mean((y[I] - fit$beta0.hat - X[I,] %*% fit$betas[,.])^2))
                          })
  min.id <- mses %>% reduce(cbind) %>% apply(1, mean) %>% which.min
  lambda.min <- lambdas[min.id]
  fit <- lasso(y, X, lambdas = lambdas, lambda = lambda.min)
  list(lambdas = lambdas, betas = fit$betas, beta0.hat = fit$beta0.hat, lambda.min = lambda.min,
       beta = fit$betas[, min.id])
}

test <- function(n, p){
  beta0 <- c(rep(1, 5), rep(0, p-5))
  X <- rmvnorm(n, rep(0,p), diag(1, p))
  y <- X %*% beta0 + rnorm(n, 0, 1)
  # list(fit = cv.lasso(y, X), fit.net = cv.glmnet(X, y))
  list(y=y, X=X)
}
#
# map(1:100, ~ { which((fits$fit$betas[,.] != 0) != (fits$fit.net$beta[,.] != 0)) } )
# map(1:100, ~ { sum(abs(fits$fit$betas[,.] - fits$fit.net$beta[,.])) } )
