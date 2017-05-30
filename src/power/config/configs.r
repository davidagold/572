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
  
  .n <- data.frame(n = c(100))
  .pz <- data.frame(pz = c(100))
  .Fstar <- data.frame(Fstar = c(0, 10, 40, 160))
  .theta0 <- data.frame(theta0 = seq(.5, 1.5, by = 1/100))
  
  .configs <- list(.n, .pz, .Fstar, .theta0)
  .configs %>%
    map(~ mutate(., id = 1)) %>%
    reduce(inner_join, by = "id") %>%
    tibble::rownames_to_column(var = "config_id") %>%
    write.csv(paste(config_dir, "configs.csv", sep="/"))
}

configs.()