# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#########################################################################
# Dependencies
suppressWarnings(library(dplyr))
suppressWarnings(library(tidyr))
suppressWarnings(library(purrr))

#########################################################################
# Utility functions

zeros <- function(p) rep(0, p)
ones <- function(p) rep(1, p)
