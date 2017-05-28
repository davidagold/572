#!/usr/bin/env Rscript
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#########################################################################
# Dependencies
suppressMessages(library(purrr))
suppressMessages(library(dplyr))

#########################################################################
# Configurations

configs. <- function() {
  args = commandArgs(trailingOnly=TRUE)
  config_dir <- args[1] %>% as.character
  
  n <- c(100, 500)
  pz <- c(100)
  Fstar <- c(0, 10, 40, 160)
  theta0 <- c(1)
  
  data.frame(
    n = n,
    pz = pz,
    Fstar = Fstar,
    theta0 = theta0
  ) %>%
    tibble::rownames_to_column(var = "config_id") %>%
    write.csv(paste(config_dir, "configs.csv", sep="/"))
}

configs.()