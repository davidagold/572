#!/usr/bin/env Rscript
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#########################################################################
# Dependencies
suppressWarnings(library(dplyr))
suppressWarnings(library(tidyr))
suppressWarnings(library(purrr))

#########################################################################
# Analysis

test <- function(estimate, truth, SE, level = .05) {
  z_star <- qnorm(1-level/2, lower.tail = TRUE)
  ifelse(abs(estimate-truth)/SE > z_star, 1, 0)
}

analyze <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  res_dir <- args[1]
  src_dir <- args[2]
  
  res <- paste(res_dir, "res.csv", sep = "/") %>% read.csv
  configs <- paste(src_dir, "config/configs.csv", sep = "/") %>% read.csv
  
  level = .05
  z_star <- qnorm(1-level/2, lower.tail = TRUE)
  
  res %>%
    inner_join(configs, by = "config_id") %>%
    mutate(rp05 = ifelse(abs(estimate-theta0)/SE > z_star, 1, 0)) %>%
    group_by(config_id, estimator) %>% 
    summarize(avg_rp05 = mean(rp05),
              RMSE = mean((estimate-theta0)^2) %>% sqrt,
              med_bias = quantile(estimate-theta0, .5),
              noinsts = sum(s.hat == 0)) %>%
    { inner_join(configs, ., by = "config_id") } %>%
    print
}

analyze()
# 
# df_est <- read.csv("res/est.csv")
# df_stats <- read.csv("res/stats.csv")
