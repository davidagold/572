# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#########################################################################
# Dependencies

#########################################################################
# Lasso estimator

lasso <- function(y, X, epsilon = .001, K = 100, thresh = 1e-7, lambdas = NULL, lambda = NULL) {
  y <- as.vector(y)
  n <- nrow(X); p <- ncol(X)
  
  y.bar <- mean(y); beta0.hat <- y.bar
  y <- y - y.bar
  X <- scale(X, scale = apply(X, 2, mysd)) %>% as.matrix
  
  if ( is.null(lambdas) ) {
    lambda.max <- max(abs(colSums(cbind(ones(n), X) * y)))/n
    lambda.min <- epsilon * lambda.max
    lambdas <- exp(seq(log(lambda.max), log(lambda.min), length.out = K))
  }
  
  if ( !is.null(lambda) ) {
    lambda.id <- map(lambdas, ~abs(.-lambda)) %>% which.min %>% as.numeric
    lambdas[lambda.id] <- lambda
  } else { 
    lambda.id <- 1 # not meant to be used 
  }

  null.dev <- deviance(lm(y ~ 1))
  null.dev
  beta.hats <- matrix(nrow = p, ncol = K)
  beta.hats[,1] <- 0
  Ls <- map(1:p, ~ setdiff(1:p, .)) %>% reduce(cbind)
  
  # outer loop over lambda path
  for ( k in 2:K ) {
    lambda <- lambdas[k]
    beta.hat <- beta.hats[,k-1] %>% as.numeric
    # print(beta.hat)
    obj0 <- Inf
    # obj1 <- sum((y - beta0.hat - X %*% beta.hat)^2)/(2*n) + lambda*sum(abs(beta.hat))
    obj1 <- sum((y - X %*% beta.hat)^2)/(2*n) + lambda*sum(abs(beta.hat))
    r <- 0
    while ( abs(obj0 - obj1) > thresh*null.dev ) {
      # print(r <- r+1)
      # inner loop over coordinates
      for ( j in 1:p ) {
        L.j <- Ls[,j]
        # y.j <- map_dbl(1:n, ~ beta0.hat + sum(X[. , L.j]*beta.hat[L.j]))
        y.j <- map_dbl(1:n, ~ sum(X[. , L.j]*beta.hat[L.j]))
        beta.hat[j] <- (sum(X[,j] * (y - y.j))/n) %>% st(lambda)
      }
      obj0 <- obj1
      # obj1 <- sum((y - beta0.hat - X %*% beta.hat)^2)/(2*n) + lambda*sum(abs(beta.hat))
      obj1 <- sum((y - X %*% beta.hat)^2)/(2*n) + lambda*sum(abs(beta.hat))
      # print(sprintf("objective function: %f", obj1))
    }
    beta.hats[,k] <- beta.hat
  }
  if ( is.null(lambda) ) { beta.hat <- NULL } else { beta.hat <- beta.hats[, lambda.id] }
  list(lambdas = lambdas, beta.hats = beta.hats, beta0.hat = beta0.hat, beta.hat = beta.hat,
       lambda.id = lambda.id)
}

#########################################################################
# CV Lasso 

cv.lasso <- function(y, X, nfolds = 10, epsilon = .001, K = 100, thresh = 1e-7, ...) {
  # ys <-  
  n <- nrow(X)
  nperfold <- n %/% nfolds 
  is <- sample(1:n)
  Is <- map(0:(nfolds-1), ~ is[(.*nperfold+1):min(((.+1)*nperfold), n)])
  
  lambda.max <- max(abs(colSums(cbind(ones(n), X) * y)))/n
  lambda.min <- epsilon * lambda.max
  lambdas <- exp(seq(log(lambda.max), log(lambda.min), length.out = K))
  
  fits <- map(1:nfolds, ~ { I <- setdiff(1:n, Is[.]); fit <- lasso(y[I], X[I,], lambdas = lambdas) ;
                      print(sprintf("%2.f%% done fitting", 100*./nfolds)); fit })
  mses <- map(1:nfolds, ~ { I <- Is[[.]]; fit <- fits[[.]];
  map_dbl(1:K, 
          ~ mean((y[I] - fit$beta0.hat - X[I,] %*% fit$beta.hats[,.])^2)) 
  } )
  min.id <- mses %>% reduce(cbind) %>% apply(1, mean) %>% which.min
  lambda.min <- lambdas[min.id]
  fit <- lasso(y, X, lambdas = lambdas, lambda = NULL)
  list(lambdas = lambdas, beta.hats = fit$beta.hats, beta0.hat = fit$beta0.hat, lambda.min = lambda.min,
       beta.hat = fit$beta.hats[, min.id])
}

