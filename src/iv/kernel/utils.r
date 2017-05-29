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

betas <- function(.fit) {
  if ( is(.fit, "glmnet") ) { 
    .as <- .fit$a0; .betas <- .fit$beta 
  } else if ( is(.fit, "cv.glmnet" )) { 
    .as <- .fit$glmnet.fit$a0; .betas <- .fit$glmnet.fit$beta 
  } else {}
  dim(.as) <- c(1, length(.as))
  rbind(.as, .betas)
}