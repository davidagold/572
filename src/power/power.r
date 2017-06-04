#!/usr/bin/env Rscript
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(dplyr)

run <- function(){
  trial_id <- Sys.getenv('SLURM_ARRAY_TASK_ID') %>% as.numeric
  args = commandArgs(trailingOnly=TRUE)
  config_id <- args[1] %>% as.numeric
  src_dir <- args[2]
  res_dir <- args[3]
  # trial_id <- 1
  # config_id <- 6
  # src_dir <- "/Users/David/Documents/UW/2016-7/spring/572/src"
  # res_dir <- "/Users/David/Documents/UW/2016-7/spring/572/res"
  source(paste(src_dir, "iv/trial.r", sep="/"), chdir=TRUE)
  source(paste(src_dir, "utils.r", sep="/"), chdir=TRUE)
  source(paste(src_dir, "dgm/iv.r", sep="/"), chdir=TRUE)
  source(paste(src_dir, "est/est.r", sep="/"), chdir=TRUE)
  
  trial(config_id, trial_id, res_dir)
}

run()