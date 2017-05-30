# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#########################################################################
# Dependencies
suppressWarnings(library(dplyr))
suppressWarnings(library(tidyr))
suppressWarnings(library(purrr))
suppressWarnings(library(methods))
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
  # config_id <- 4
  # trial_id <- 1
  res_dir <- args[2]
  
  # set up containers
  R <- 6
  res <- data.frame(
    config_id = rep(config_id, R),
    trial_id = rep(trial_id, R),
    estimator = numeric(R),
    theta.hat = numeric(R),
    SE_theta.hat = numeric(R),
    sigma.hat = numeric(R),
    s.hat = numeric(R),
    s.op = numeric(R)    
  )
  # Generate data
  obs <- obs.(config_id); n <- obs$n; pz <- obs$pz
  
  # helper function for recording results
  record <- function(res, r, obs, .fit.beta=NULL, S.op=NULL, .fit.delta, ...) {
    y <- obs$y; x <- obs$x; Z <- obs$Z; pz <- obs$pz
    if ( !is.null(.fit.beta) ) {
      S.op <- 2:(pz+1)
      fit. <- .fit.beta$fit.; lambda. <- .fit.beta$lambda.
      fit.beta <- .fit.beta$fit.(x = Z(S.op), y = x, intercept = TRUE, lambda.min.ratio = 0.0001)
      beta.hat <- fit.beta %>% 
        predict(type = "coefficients", lambda = lambda.(obs), exact = TRUE, 
                x = Z(S.op), y = x) %>%
        as.numeric
      
      S.hat <- which(beta.hat[2:(pz+1)] != 0); S.op <- S.hat
      s.hat <- length(S.hat); s.op <- s.hat
      
      if ( s.op == 0 ) {
        id.mod <- map_dbl(1:length(fit.beta$lambda), 
                          ~ betas(fit.beta)[,.] %>% { which(.[2:(pz+1)] != 0) } %>% length) %>%
                          { which(. > 0) } %>% min 
        beta.hat_mod <- betas(fit.beta)[,id.mod] %>% as.numeric
        S.op <- which(beta.hat_mod != 0)
        s.op <- length(S.op)
      }
    } else if ( !is.null(S.op) ) {
      s.hat <- NA; s.op <- length(S.op)
    }
    fit.delta <- .fit.delta(obs, S.op, ...)
    
    res$estimator[r] <- fit.delta$estimator
    res$theta.hat[r] <- fit.delta$theta.hat
    res$SE_theta.hat[r] <- fit.delta$SE_theta.hat
    res$sigma.hat[r] <- fit.delta$sigma.hat
    res$s.hat[r] <- s.hat
    res$s.op[r] <- s.op
    assign('res', res, envir=environment(record))
  }
  
  r <- 0
  # 2SLS(All)
  if ( n == 100 ) { S.op_2SLS <- sample(1:pz, 98, replace = FALSE) } else { S.op_2SLS <- 1:pz }
  record(res, r<-r+1, obs, S.op = S.op_2SLS, .fit.delta = fit.IV, estimator = "2SLS(All)")
  
  # Fuller(All)
  if ( n == 100 ) { S.op_FL <- sample(1:pz, 98, replace = FALSE) } else { S.op_FL <- 1:pz }
  record(res, r<-r+1, obs, S.op = S.op_FL, 
         .fit.delta = fit.Fuller, estimator = "Fuller(All)")
  
  # IV(Lasso-IL)
  record(res, r<-r+1, obs, .fit.beta = list(fit. = glmnet, lambda. = lambda.IL),
         .fit.delta = fit.IV, estimator = "IV(Lasso-IL)")
  
  # Fuller(Lasso-IL)
  record(res, r<-r+1, obs, .fit.beta = list(fit. = glmnet, lambda. = lambda.IL),
         .fit.delta = fit.Fuller, estimator = "Fuller(Lasso-IL)")
  
  # IV(Lasso-CV)
  record(res, r<-r+1, obs, .fit.beta = list(fit. = cv.glmnet, lambda. = lambda.CV),
         .fit.delta = fit.IV, estimator = "IV(Lasso-CV)")
  
  # Fuller(Lasso-CV)
  record(res, r<-r+1, obs, .fit.beta = list(fit. = cv.glmnet, lambda. = lambda.CV),
         .fit.delta = fit.Fuller, estimator = "Fuller(Lasso-CV)")
 
  res %>%
    write.csv(paste(res_dir, config_id, sprintf("res%d.csv", trial_id), sep = "/"))
    # print
}


