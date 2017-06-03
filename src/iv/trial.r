# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#!/usr/bin/env Rscript

#########################################################################
# Simulation trial

trial <- function(config_id, trial_id, res_dir) {
  # set up containers
  R <- 7
  res <- data.frame(
    config_id = rep(config_id, R),
    trial_id = rep(trial_id, R),
    estimator = numeric(R),
    statistic = numeric(R),
    SE_theta.hat = numeric(R),
    sigma.hat = numeric(R),
    s.hat = numeric(R),
    s.op = numeric(R)
  )
  # Generate data
  obs <- obs.(config_id)
  n <- obs$n; pz <- obs$pz
  
  # helper function for recording results
  record <- function(res, r, obs, fit.beta=NULL, S.op=NULL, est.delta, ...) {
    print(r)
    y <- obs$y; x <- obs$x; Z <- obs$Z; pz <- obs$pz
    if ( !is.null(fit.beta) ) {
      beta.hat <- fit.beta$beta; beta.hats <- fit.beta$betas; lambdas <- fit.beta$lambdas
      S.hat <- which(beta.hat != 0); S.op <- S.hat
      s.hat <- length(S.hat); s.op <- s.hat
      
      if ( s.op == 0 ) {
        id.op <- map_dbl(1:length(lambdas), 
                         ~ beta.hats[,.] %>% { which(. != 0) } %>% length) %>%
                         { which(. > 0) } %>% min 
        beta.hat.op <- beta.hats[,id.op] %>% as.numeric
        S.op <- which(beta.hat.op != 0)
        s.op <- length(S.op)
        # fit.delta <- est.delta(obs, S.op, ...)
      }
    } else if ( !is.null(S.op) ) {
      s.hat <- NA; s.op <- length(S.op)
    }
    fit.delta <- est.delta(y = y, x = x, Z = as.matrix(Z[, S.op]), ...)
    
    res$estimator[r] <- fit.delta$estimator
    res$statistic[r] <- fit.delta$theta.hat
    res$SE_theta.hat[r] <- fit.delta$SE_theta.hat
    res$sigma.hat[r] <- fit.delta$sigma.hat
    res$s.hat[r] <- s.hat
    res$s.op[r] <- s.op
    assign('res', res, envir=environment(record))
  }
  record_sup.score <- function(res, r, obs, level = .05) {
    res$estimator[r] <- "Sup-Score"
    res$statistic[r] <- test.sup.score(obs, a = 1, level = level)
    res$SE_theta.hat[r] <- NA
    res$sigma.hat[r] <- NA
    res$s.hat[r] <- NA
    res$s.op[r] <- NA
    assign('res', res, envir=environment(record))
  }
  
  r <- 0
  # 2SLS(All)
  if ( n == 100 ) { S.op_2SLS <- sample(1:pz, 98, replace = FALSE) } else { S.op_2SLS <- 1:pz }
  record(res, r<-r+1, obs, S.op = S.op_2SLS, est.delta = tsls, estimator = "2SLS(All)")
  # Fuller(All)
  if ( n == 100 ) { S.op_FL <- sample(1:pz, 98, replace = FALSE) } else { S.op_FL <- 1:pz }
  record(res, r<-r+1, obs, S.op = S.op_FL, est.delta = fuller, estimator = "Fuller(All)")

  lambda <- lambda_(x = obs$x, Z = obs$Z)
  fit.beta.IL <- lasso(y = obs$x, X = obs$Z, lambda = lambda)
  fit.beta.CV <- cv.lasso(y = obs$x, X = obs$Z, epsilon = .1)
  fit.beta.CV
  # list(x = obs$x, y = obs$y, Z = obs$Z)

  # IV(Lasso-IL)
  record(res, r<-r+1, obs, fit.beta = fit.beta.IL, est.delta = tsls, estimator = "IV(Lasso-IL)")
  # Fuller(Lasso-IL)
  record(res, r<-r+1, obs, fit.beta = fit.beta.IL, est.delta = fuller, estimator = "Fuller(Lasso-IL)")

  # IV(Lasso-CV)
  record(res, r<-r+1, obs, fit.beta = fit.beta.CV, est.delta = tsls, estimator = "IV(Lasso-CV)")
  # Fuller(Lasso-CV)
  record(res, r<-r+1, obs, fit.beta = fit.beta.CV, est.delta = fuller, estimator = "Fuller(Lasso-CV)")

  # Sup-Score
  record_sup.score(res, r<-r+1, obs)

  res %>%
    write.csv(paste(res_dir, config_id, sprintf("res%d.csv", trial_id), sep = "/"))
  # print
}


