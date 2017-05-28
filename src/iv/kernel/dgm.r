# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#########################################################################
# Dependencies

suppressWarnings(library(dplyr))
suppressWarnings(library(tidyr))
suppressWarnings(library(purrr))
suppressWarnings(library(mvtnorm))

#########################################################################
# Data-generating mechanism

zeros <- function(p) rep(0, p)
ones <- function(p) rep(1, p)

Z. <- function(n, pz) {
  Sigma_z <- .5^abs(outer(1:pz, 1:pz, "-"))
  Z <- rmvnorm(n, zeros(pz), Sigma_z)
  Z
}

hv. <- function(n, pz) {
  sigma0_h <- 1; sigma0_v <- 1; sigma0_hv <- .3
  Sigma0_hv <- matrix(c(sigma0_h^2, sigma0_hv,
                        sigma0_hv, sigma0_v^2), nrow=2)
  hv <- rmvnorm(n, zeros(2), Sigma0_hv)
  list(h = hv[,1], v = hv[,2])
}

yx. <- function(Z, h, v, beta0, theta0) {
  x <- Z %*% beta0 + v
  y <- x * theta0 + h
  list(x = x, y = y)
}

