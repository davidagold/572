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
  
  df_est <- paste(res_dir, "est.csv", sep = "/") %>% read.csv
  df_stats <- paste(res_dir, "stats.csv", sep = "/") %>% read.csv
  configs <- paste(src_dir, "config/configs.csv", sep = "/") %>% read.csv
  
  level = .05
  z_star <- qnorm(1-level/2, lower.tail = TRUE)
  
  df_est %>%
    inner_join(configs, by = "config_id") %>%
    group_by(config_id) %>% 
    summarize(rp_05 = ifelse(abs(estimate-truth)/SE > z_star, 1, 0)) %>%
    print
}

analyze()