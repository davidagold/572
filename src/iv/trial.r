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
    s.op = numeric(R),
    lambda.min = numeric(R),
    lambda.max = numeric(R),
    lambda = numeric(R)
  )
  # Generate data
  obs <- obs.(config_id)
  n <- obs$n; pz <- obs$pz
  
  # helper function for recording results
  record <- function(res, r, obs, fit.beta=NULL, S.op=NULL, est.delta, lambda.min=NA,lambda.max=NA, lambda=NA, ...) {
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
    fit.delta <- est.delta(y = y, X = x, Z = as.matrix(Z[, S.op]), ...)
    
    res$estimator[r] <- fit.delta$estimator
    res$statistic[r] <- fit.delta$theta.hat
    res$SE_theta.hat[r] <- fit.delta$SE_theta.hat
    res$sigma.hat[r] <- fit.delta$sigma.hat
    res$s.hat[r] <- s.hat
    res$s.op[r] <- s.op
    res$lambda.min[r] <- lambda.min
    res$lambda.max[r] <- lambda.max
    res$lambda[r] <- lambda
    assign('res', res, envir=environment(record))
  }
  record_sup.score <- function(res, r, obs, level = .05) {
    res$estimator[r] <- "Sup-Score"
    res$statistic[r] <- test.sup.score(obs, a = 1, level = level)
    res$SE_theta.hat[r] <- NA
    res$sigma.hat[r] <- NA
    res$s.hat[r] <- NA
    res$s.op[r] <- NA
    res$lambda.min[r] <- NA
    res$lambda.max[r] <- NA
    res$lambda[r] <- NA
    assign('res', res, envir=environment(record))
  }
  
  r <- 0
  # 2SLS(All)
  if ( n == 100 ) { S.op_2SLS <- sample(1:pz, 98, replace = FALSE) } else { S.op_2SLS <- 1:pz }
  record(res, r<-r+1, obs, S.op = S.op_2SLS, est.delta = tsls, estimator = "2SLS(All)")
  # Fuller(All)
  if ( n == 100 ) { S.op_FL <- sample(1:pz, 98, replace = FALSE) } else { S.op_FL <- 1:pz }
  record(res, r<-r+1, obs, S.op = S.op_FL, est.delta = fuller, estimator = "Fuller(All)")

  lambdas <- lambda_(x = obs$x, Z = obs$Z)
  sprintf("lambda.dat = %2.6f", lambdas$lambda.dat) %>% print
  sprintf("lambda.thr = %2.6f", lambdas$lambda.thr) %>% print
  
  fit.beta.IL <- lasso(y = obs$x, X = obs$Z, lambda = lambdas$lambda.thr, standardize=FALSE)
  lambdas.IL <- fit.beta.IL$lambdas
  lambda.max.IL <- lambdas.IL[1]
  lambda.min.IL <- lambdas.IL[length(lambdas.IL)]
  lambda.IL <- lambdas$lambda.thr
  sprintf("lambda.max (IR) = %2.6f", fit.beta.IL$lambdas[1]) %>% print
  sprintf("lambda.min (IR) = %2.6f", fit.beta.IL$lambdas[length(fit.beta.IL$lambdas)]) %>% print
  
  fit.beta.CV <- cv.lasso(y = obs$x, X = obs$Z, epsilon = .4, standardize=FALSE)
  lambda.CV <- fit.beta.CV$lambda.min
  sprintf("lambda.max (CV) = %2.6f", fit.beta.CV$lambdas[1]) %>% print
  sprintf("lambda.min (CV) = %2.6f", fit.beta.CV$lambdas[length(fit.beta.IL$lambdas)]) %>% print
  print("S.hat (CV) = "); which(fit.beta.CV$beta != 0) %>% print
  
  lambdas.CV <- fit.beta.CV$lambdas
  lambda.max.CV <- lambdas.CV[1]
  lambda.min.CV <- lambdas.CV[length(lambdas.CV)]

  # IV(Lasso-IL)
  record(res, r<-r+1, obs, fit.beta = fit.beta.IL, est.delta = tsls, lambda.min = lambda.min.IL, lambda.max = lambda.max.IL,
         lambda = lambda.IL, estimator = "IV(Lasso-IL)")
  # Fuller(Lasso-IL)
  record(res, r<-r+1, obs, fit.beta = fit.beta.IL, est.delta = fuller, lambda.min = lambda.min.IL, lambda.max = lambda.max.IL, 
         lambda = lambda.IL, estimator = "Fuller(Lasso-IL)")

  # IV(Lasso-CV)
  record(res, r<-r+1, obs, fit.beta = fit.beta.CV, est.delta = tsls, lambda.min = lambda.min.CV, lambda.max = lambda.max.CV,
         lambda = lambda.CV, estimator = "IV(Lasso-CV)")
  # Fuller(Lasso-CV)
  record(res, r<-r+1, obs, fit.beta = fit.beta.CV, est.delta = fuller, lambda.min = lambda.min.CV, lambda.max = lambda.max.CV,
         lambda = lambda.CV, estimator = "Fuller(Lasso-CV)")

  # Sup-Score
  record_sup.score(res, r<-r+1, obs)

  res %>%
    # write.csv(paste(res_dir, config_id, sprintf("res%d.csv", trial_id), sep = "/"))
  print
}


