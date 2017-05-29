# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#########################################################################
# Dependencies
suppressWarnings(library(dplyr))
suppressWarnings(library(tidyr))
suppressWarnings(library(purrr))
source("utils.r")
source("dgm.r")
source("estimation.r")
# source("kernel/utils.r")
# source("kernel/estimation.r")
# source("kernel/dgm.r")

#########################################################################
# Simulation trial

trial <- function(res_dir) {
  # obtain config and trial ids and results directory
  args = commandArgs(trailingOnly=TRUE)
  config_id <- args[1] %>% as.numeric
  trial_id <- Sys.getenv('SLURM_ARRAY_TASK_ID') %>% as.numeric
  # config_id <- 1
  # trial_id <- 1
  res_dir <- args[2]
  
  # set up containers
  R <- 2
  res <- data.frame(
    config_id = rep(config_id, R),
    trial_id = rep(trial_id, R),
    estimator = numeric(R),
    estimate = numeric(R),
    SE = numeric(R),
    sigma0h.hat = numeric(R),
    s.hat = numeric(R),
    s.hat_mod = numeric(R)    
  )
  # Generate data
  obs <- obs.(config_id)
  y <- obs$y; x <- obs$x; Z <- obs$Z
  n <- obs$n; pz <- obs$pz
  
  # helper function for recording results
  record <- function(res, r, estimator, first_stage) {
    Lasso.fit <- first_stage$Lasso.fit; 
    beta0.hat <- first_stage$beta0.hat; 
    lambda <- first_stage$lambda
    S.hat <- which(beta0.hat[2:(pz+1)] != 0); S.hat_mod <- S.hat
    s.hat <- length(S.hat); s.hat_mod <- s.hat
    if ( s.hat_mod == 0 ) {
      min_id <- map_dbl(1:length(Lasso.fit$lambda), 
                      ~ Lasso.fit$beta[,.] %>% { which(. != 0) } %>% length) %>%
        { which(. != 0) } %>% min 
      
      beta0.hat_mod <- c(Lasso.fit$a0[min_id], Lasso.fit$beta[,min_id] %>% as.numeric)
      S.hat_mod <- which(beta0.hat_mod[2:(pz+1)] != 0)
      s.hat_mod <- length(S.hat_mod)
    }
    Z.S <- Z[, S.hat_mod]
    IV.fit <- fit.IV(y, x, Z.S)
    
    res$estimator[r] <- estimator
    res$estimate[r] <- IV.fit$theta0.hat
    res$SE[r] <- IV.fit$SE_theta0.hat
    res$sigma0h.hat[r] <- IV.fit$sigma0h.hat
    res$s.hat[r] <- s.hat
    res$s.hat_mod[r] <- s.hat_mod
    
    assign('res', res, envir=environment(record))
  }
  
  # IV-Lasso
  sigma0_v.hat <- sigma0_v.hat_iter(x, Z)
  lambda_IL <- .lambda(sigma0.hat = sigma0_v.hat, Z = Z)
  Lasso.fit_IL <- glmnet(Z, x, intercept = TRUE)
  beta0.hat_IL <- predict(Lasso.fit_IL, type="coefficients", 
                          s=lambda_IL, exact=TRUE, x=Z, y=x) %>% as.numeric
  record(res, r=1, estimator = "IV-Lasso", 
         first_stage = list(Lasso.fit = Lasso.fit_IL, 
                            beta0.hat = beta0.hat_IL,
                            lambda = lambda_IL))
  
  # IV-Lasso-CV
  Lasso.fit_CV <- cv.glmnet(Z, x, intercept = TRUE, lambda.min.ratio = 0.0001)
  beta0.hat_CV <- predict(Lasso.fit_CV, type = "coefficients",
                          s = "lambda.min") %>% as.numeric
  lambda_CV <- Lasso.fit_CV$lambda.min
  record(res, r=2, estimator="IV-Lasso-CV",
         first_stage = list(Lasso.fit = Lasso.fit_CV$glmnet.fit,
                            beta0.hat = beta0.hat_CV,
                            lambda = lambda_CV))
  res %>%
    write.csv(paste(res_dir, config_id, sprintf("res%d.csv", trial_id), sep = "/"))
  # print
}
