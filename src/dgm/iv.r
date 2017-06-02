# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#########################################################################
# Dependencies
suppressWarnings(library(mvtnorm))

#########################################################################
# Data-generating mechanism

Z. <- function(n, pz, Sigma_z) {
  Z <- rmvnorm(n, zeros(pz), Sigma_z)
  Z
}

hv. <- function(n, sigma0_h, sigma0_v, sigma0_hv) {
  Sigma0_hv <- matrix(c(sigma0_h^2, sigma0_hv*sigma0_h*sigma0_v,
                        sigma0_hv*sigma0_h*sigma0_v, sigma0_v^2), nrow=2)
  hv <- rmvnorm(n, zeros(2), Sigma0_hv)
  list(h = hv[,1], v = hv[,2])
}

yx. <- function(Z, h, v, beta0, theta0) {
  x <- Z %*% beta0 + v
  y <- x * theta0 + h
  list(x = x, y = y)
}

obs. <- function(.config_id) {
  configs <- read.csv("config/configs.csv")
  config <- filter(configs, config_id == .config_id)
  n <- config$n; pz <- config$pz; Fstar <- config$Fstar; theta0 <- config$theta0
  
  Sigma_z <- .5^abs(outer(1:pz, 1:pz, "-"))
  if ( Fstar == 0 ) {
    beta0 = rep(0, pz)
    sigma0_v = 1
  } else {
    beta0 <- map_dbl(1:pz, ~.7^(.-1))
    sigma0_v <- ((n * t(beta0) %*% Sigma_z %*% beta0 ) / (Fstar * sum(beta0^2) )) %>% sqrt
  }
  sigma0_h <- 1; sigma0_hv <- .3
  
  Z <- Z.(n, pz, Sigma_z)
  hv <- hv.(n, sigma0_h, sigma0_v, sigma0_hv)
  yx <- yx.(Z, hv$h, hv$v, beta0, theta0)
  y <- as.vector(yx$y); x <- as.vector(yx$x)
  
  obs <- list(n = n, pz = pz, y = y, x = x, X = X, Z = Z)
  obs
}

