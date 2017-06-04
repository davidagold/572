#!/usr/bin/env Rscript
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(dplyr)
source("trial.r", chdir = TRUE)

run <- function(){
  trial_id <- Sys.getenv('SLURM_ARRAY_TASK_ID') %>% as.numeric
  args = commandArgs(trailingOnly=TRUE)
  config_id <- args[1] %>% as.numeric
  src_dir <- args[2]
  res_dir <- args[3]
  # trial_id <- 1
  # config_id <- 3
  # src_dir <- "/Users/David/Documents/UW/2016-7/spring/572/src"
  # res_dir <- "/Users/David/Documents/UW/2016-7/spring/572/res"
  source(paste(src_dir, "utils.r", sep="/"), chdir=T)
  source(paste(src_dir, "dgm/iv.r", sep="/"), chdir=T)
  source(paste(src_dir, "est/est.r", sep="/"), chdir=T)
  
  trial(config_id, trial_id, res_dir)
}

run()

