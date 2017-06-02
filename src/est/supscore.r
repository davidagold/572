# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#########################################################################
# Dependencies

#########################################################################
# Sup-Score

test.sup.score <- function(obs, a = 1, level = 0.5, m = 500) {
  y <- obs$y; x <- obs$x; Z <- obs$Z
  n <- nrow(Z); pz <- ncol(Z)
  Z <- scale(Z, scale = apply(Z, 2, mysd))
  sup.score <- map_dbl(1:pz, ~ abs(sum((y - x*a) * Z[,.])) / 
                         sqrt(mean((y - x*a)^2 * Z[,.]^2))) %>% max
  
  q.star <- map_dbl(1:m, ~ (Z * rnorm(n, 0, 1)) %>%
            apply(2, function(x) { abs(sum(x))/sqrt(mean(x^2)) }) %>%
            max) %>%
    as.numeric %>%
    quantile(1 - level)

  ifelse(sup.score > q.star, 1, 0)
}