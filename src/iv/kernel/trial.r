# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#########################################################################
# Dependencies

suppressWarnings(library(dplyr))
suppressWarnings(library(tidyr))
suppressWarnings(library(purrr))
source("dgm.r")
source("estimation.r")

#########################################################################
# Simulation trial

trial <- function(res_dir) {
  # obtain config and trial ids and results directory
  args = commandArgs(trailingOnly=TRUE)
  config_id <- args[1] %>% as.numeric
  res_dir <- args[2]
  trial_id <- Sys.getenv('SLURM_ARRAY_TASK_ID') %>% as.numeric
  
  # set up containers
  R <- 1
  .config_id <- rep(config_id, R)
  .trial_id <- rep(trial_id, R)
  .estimator <- numeric(R)
  .estimate <- numeric(R)
  .SE <- numeric(R)
  .rmse <- numeric(R)
  .shat <- numeric(R)
  
  # Generate data
  obs <- obs.(config_id)
  y <- obs$y; x <- obs$x; Z <- obs$Z
  n <- obs$n; pz <- obs$pz
  
  # 2SLS, CV Lasso
  fit_fs <- cv.glmnet(Z, x, intercept = FALSE)
  beta0hat <- predict(fit_fs, type = "coefficients", s = "lambda.min")[2:(pz+1)] %>% as.numeric
  Shat_CV <- which(beta0hat != 0)
  shat_CV <- length(Shat_CV)
  print(Shat_CV)
  Z_ps_CV <- Z[, Shat_CV]
  fit_tsls_CV <- theta0_tsls.(y, x, Z_ps_CV)
  theta0_tsls_CV <- fit_tsls_CV$theta0_hat
  sigma0_htsls_CV <- fit_tsls_CV$sigma0_hhat
  var_theta0_tsls_CV <- fit_tsls_CV$var_theta0_hat
  
  SE_tsls_CV <- diag(var_theta0_tsls_CV) * n %>% sqrt
  rmse_tsls_CV <- (y - x %*% theta0_tsls_CV)^2 %>% mean %>% sqrt
  
  r <- 1
  .estimator[r] <- "IV-Lasso-CV"
  .estimate[r] <- theta0_tsls_CV
  .SE[r] <- SE_tsls_CV
  .rmse[r] <- rmse_tsls_CV
  .shat[r] <- shat_CV
  
  df_est <- data.frame(
    config_id = .config_id,
    trial_id = .trial_id,
    estimator = .estimator,
    estimate = .estimate,
    SE = .SE
  )
  df_stats <- data.frame(
    config_id = .config_id,
    trial_id = .trial_id,
    estimator = .estimator,
    rmse = .rmse,
    shat = .shat
  )
  
  # list(df_est = df_est, df_stats = df_stats)
  write.csv(df_est, paste(res_dir, "/", config_id, "/est/est", trial_id, ".csv", sep=""))
  write.csv(df_stats, paste(res_dir, "/", config_id, "/stats/stats", trial_id, ".csv", sep=""))
}
