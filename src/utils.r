# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#########################################################################
# Dependencies
suppressWarnings(library(dplyr))
suppressWarnings(library(tidyr))
suppressWarnings(library(purrr))
suppressWarnings(library(methods))
suppressWarnings(library(mvtnorm))

#########################################################################
# Utility functions

zeros <- function(p) rep(0, p)
ones <- function(p) rep(1, p)
P <- function(Z, b, Zt = t(Z)) { Z %*% solve(Zt %*% Z, Zt %*% b) }
my.sd <- function(x) { mean((x - mean(x))^2) %>% sqrt }
st <- function(x, lambda) { sign(x)*max(abs(x) - lambda, 0) }

my.scale <- function(X) { scale(X, scale = apply(X, 2, my.sd))}