#!/usr/bin/env Rscript

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(dplyr)
l2 <- function(x) sqrt(sum(x^2))

args <- commandArgs(trailingOnly = TRUE)
args[1] %>% 
  read.csv %>%
  group_by(j, estimator, sigma0) %>%
  summarize(avg_bias = mean(betahat_j - beta0_j)) %>%
  group_by(estimator, sigma0) %>%
  summarize(l2_avg_bias = l2(avg_bias)) %>%
  print