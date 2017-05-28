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
  args = commandArgs(trailingOnly=TRUE)
  # config_id. <- args[1] %>% as.numeric
  trial <- Sys.getenv('SLURM_ARRAY_TASK_ID') %>% as.numeric
  n <- 30
  pz <- 10
  theta0 <- 1
  
  beta0 <- map_dbl(1:pz, ~.7^(.-1))
  
  # Generate data
  Z <- Z.(n, pz)
  hv <- hv.(n, pz)
  h <- hv$h; v <- hv$v
  yx <- yx.(Z, h, v, beta0, theta0)
  y <- yx$y; x <- yx$x
  
  R <- 1
  .trial <- numeric(R)
  .estimator <- numeric(R)
  .estimate <- numeric(R)
  .SE <- numeric(R)
  .rmse <- numeric(R)
  .shat <- numeric(R)
  
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
  print(var_theta0_tsls_CV)
  
  SE_tsls_CV <- diag(var_theta0_tsls_CV) * n %>% sqrt
  print(SE_tsls_CV)
  
  r <- 1
  .trial[r] <- .trial
  .estimator[r] <- "IV-Lasso-CV"
  .estimate[r] <- theta0_tsls_CV
  .SE[r] <- SE_tsls_CV
  
  rmse_tsls_CV <- (y - x %*% theta0_tsls_CV)^2 %>% mean %>% sqrt
  .rmse[r] <- rmse_tsls_CV
  .shat[r] <- shat_CV
  
  df_est <- data.frame(
    trial = .trial,
    estimator = .estimator,
    estimate = .estimate,
    SE = .SE
  )
  df_stats <- data.frame(
    trial = .trial,
    estimator = .estimator,
    rmse = .rmse,
    shat = .shat
  )
  
  # list(df_est = df_est, df_stats = df_stats)
  write.csv(df_est, paste(res_dir, "/est/est", trial, ".csv", sep=""))
  write.csv(df_stats, paste(res_dir, "/stats/stats", trial, ".csv", sep=""))
}
