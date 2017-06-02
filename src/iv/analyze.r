#!/usr/bin/env Rscript
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#########################################################################
# Dependencies
suppressWarnings(library(dplyr))
suppressWarnings(library(tidyr))
suppressWarnings(library(purrr))

#########################################################################
# Analysis

test <- function(theta.hat, truth, SE, level = .05) {
  z_star <- qnorm(1-level/2, lower.tail = TRUE)
  ifelse(abs(theta.hat-truth)/SE > z_star, 1, 0)
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
    filter(estimator != "Sup-Score") %>%
    inner_join(configs, by = "config_id") %>%
    mutate(rp05 = ifelse(abs(statistic-theta0)/SE_theta.hat > z_star, 1, 0)) %>%
    group_by(config_id, estimator) %>% 
    summarize(avg_rp05 = mean(na.omit(rp05)),
              RMSE = mean((statistic-theta0)^2) %>% sqrt,
              med_bias = quantile(statistic-theta0, .5),
              noinsts = sum(s.hat == 0)) %>%
    xtable(digits = c(rep(0,3), rep(3,4))) %>%
    print
    
    { inner_join(configs, ., by = "config_id") } %>%
    xtable(digits = c(rep(0,9), rep(3,4))) %>%
    print(include.rownames == FALSE)
  
  my_order <- c("2SLS(All)", "Fuller(All)", "IV(Lasso-IL)", "Fuller(Lasso-IL)",
                "IV(Lasso-CV)", "Fuller(Lasso-CV)", "Sup-Score")
  

}

res.BCH <- list(
Fstar0 = matrix(c(
# "2SLS(All)"
0.318, 0.305, 0.862, NA,
0.312, 0.297, 0.852, NA,
# "Fuller(All)" 
2.398, 0.248, 0.704, NA,
1.236, 0.318, 0.066, NA,
# "IV(Lasso-IL)" 
0.511, 0.338, 0.014, 455,
0.477, 0.296, 0.012, 486,
# "Fuller(Lasso-IL)"
0.509, 0.338, 0.010, 455,
0.477, 0.296, 0.012, 486,
# "IV(Lasso-CV)"
0.329, 0.301, 0.652, 0,
0.478, 0.299, 0.064, 348,
# "Fuller(Lasso-CV)"
0.359, 0.305, 0.384, 0,
0.474, 0.299, 0.054, 348,
# "Sup-Score" 
NA, NA, 0.004, NA,
NA, NA, 0.010, NA
), ncol=4, byrow = T),

Fstar10 = matrix(c(
0.058,0.058,0.806, NA,
0.026,0.025,0.808, NA,
# FULL(All),
0.545,0.050,0.690,NA,
0.816,0.006,0.052,NA,
# IV-Lasso,
0.055,0.020,0.042,147,
0.027,0.009,0.056,160,
# FULL-Lasso,
0.054,0.020,0.032,147,
0.027,0.009,0.044,160,
# IV-Lasso-CV,
0.052,0.024,0.072,10,
0.027,0.009,0.054,202,
# FULL-Lasso-CV,
0.051,0.022,0.068,10,
0.027,0.009,0.044,202,
# Sup-Score,
NA, NA, 0.006, NA,
NA, NA, 0.004, NA
), ncol=4, byrow=T),

Fstar40 = matrix(c(
0.081,0.072,0.626,NA,
0.036,0.032,0.636,NA,
# FULL(All),
0.951,0.050,0.690, NA,
0.038,0.000,0.036,NA,
# IV-Lasso,
0.051,0.012,0.048,1,
0.022,0.003,0.048,0,
# FULL-Lasso,
0.051,0.011,0.046,1,
0.022,0.002,0.038,0,
# IV-Lasso-CV,
0.048,0.016,0.058,0,
0.022,0.004,0.052,0,
# FULL-Lasso-CV,
0.049,0.014,0.050,0,
0.022,0.003,0.042,0,
# Sup-Score,
NA,NA,0.004,NA,
NA,NA,0.006,NA),
ncol=4, byrow=T),

Fstar160=matrix(c(
0.075,0.062,0.306,NA,
0.034,0.029,0.334,NA,
# FULL(All),
1.106,0.023,0.622,NA,
0.026,0.002,0.044,NA,
# IV-Lasso,
0.049,0.005,0.064,0,
0.022,0.002,0.044,0,
# FULL-Lasso,
0.049,0.002,0.056,0,
0.022,0.001,0.040,0,
# IV-Lasso-CV,
0.048,0.006,0.054,0,
0.022,0.002,0.040,0,
# FULL-Lasso-CV,
0.049,0.003,0.048,0,
0.022,0.000,0.038,0,
# Sup-Score,
NA,NA,0.004,NA,
NA,NA,0.010,NA
), ncol=4, byrow=T)) %>%
  map2(c(0, 10, 40, 160), ~ { as.data.frame(.x) %>% { colnames(.)=c("rmse_bch", "med.bias_bch", "rp05_bch", "noinsts_bch"); . } %>%
         mutate(Fstar = .y, estimator = rep(my_order, each=2)) }) %>%
  map(~ mutate(., n = rep(c(100, 500), 7))) %>%
  reduce(rbind)

tbl <- res %>%
  inner_join(configs, by = "config_id") %>%
  # filter(n == 100) %>%
  mutate(estimator = factor(estimator, levels = my_order)) %>%
  select(n, estimator, statistic, SE_theta.hat, s.hat, theta0, Fstar) %>%
  mutate(reject = ifelse(estimator == "Sup-Score",
                         statistic,
                         ifelse(abs(statistic-theta0)/SE_theta.hat > z_star, 1, 0))) %>%
  group_by(n, Fstar, estimator) %>%
  summarize(rmse = mean((statistic-theta0)^2) %>% sqrt,
            med.bias = quantile(statistic-theta0, .5),
            rp05 = mean(na.omit(reject)),
            noinsts = sum(s.hat == 0)) %>%
  mutate_at(vars(rmse, med.bias, noinsts), funs(ifelse(estimator == "Sup-Score", NA, .))) %>%
  inner_join(res.BCH, by=c("n", "estimator", "Fstar")) %>%
  select(n, Fstar, estimator, rmse, rmse_bch, med.bias, med.bias_bch, rp05, rp05_bch,
         noinsts, noinsts_bch)

tbl %>% ungroup %>%
  filter(n==100) %>%
  arrange(Fstar, estimator) %>%
  select(-n, -Fstar) %>%
  xtable(digits = c(0,0,rep(3, 7),0)) %>%
  print(include.rownames = F)

arrange(Fstar, estimator) %>% ungroup %>%
  select(-Fstar) %>%
  xtable(digits = c(0, 0, 3, 3, 3, 3)) %>% 
  print(include.rownames = F)
