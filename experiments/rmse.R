# Experiment the approximations MSE, MAE and RMSE errors
library(rstanarm)
library(bayesboot)
library(tidyverse)
source('experiments/util.R')

n = c(20, 50, 100, 200)
p = c(2, 3, 5)
rho = 0
sigma = 1
n_reps = 4

sim = expand_grid(n, p, rho, sigma, n_reps) %>% split(seq(nrow(.)))

for (i in seq_along(sim)) {
  print(paste0('Iteration :', i))
  conf = sim[[i]]
  var_diff_bb = numeric()
  var_diff_taylor = numeric()
  var_bb = numeric()
  var_taylor = numeric()
  for (j in 1:n_reps) {
    X = get_X(conf$n, conf$p, conf$rho)
    eps = rnorm(conf$n)*conf$sigma
    beta = rnorm(conf$p)
    y = X %*% beta + eps
    df = data.frame(y, X)
    
    fit1 = stan_glm(y ~ X1, data = df, refresh=0)
    fit2 = stan_glm(y ~ X1 + X2, data = df, refresh=0)
    epred1 = apply(posterior_epred(fit1), 2, mean)
    epred2 = apply(posterior_epred(fit2), 2, mean)
    
    var_diff_bb[j] = var_rmse_diff_bb(y, epred1, epred2)
    var_diff_taylor[j] = var_rmse_diff_taylor(y, epred1, epred2)
    var_bb[j] = var_rmse_bb(y, epred2)
    var_taylor[j] = var_rmse_taylor(y, epred2)
  }
  
  sim[[i]]$var_rmse_diff_bb = list(var_diff_bb)
  sim[[i]]$var_rmse_diff_taylor = list(var_diff_taylor)
  sim[[i]]$var_rmse_bb = list(var_bb)
  sim[[i]]$var_rmse_taylor = list(var_taylor)
}

# qplot(var_rmse_diff, var_rmse_diff_taylor) + geom_abline()

