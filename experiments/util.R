############ Data Generation ##############

get_X = function(n, p, rho) {
  Sigma = matrix(rho, nrow = p, ncol = p)
  for (i in 1:p) {
    Sigma[i, i] = 1
  }
  MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
}

############ Loo predictions ################
loo_epred = function(fit) {
  ypred = posterior_epred(fit)
  ll = log_lik(fit)
  M = length(fit$stanfit@sim$n_save)
  S = fit$stanfit@sim$n_save[[1]]
  r_eff = relative_eff(exp(ll), chain_id = rep(1:M, each = S))
  psis_object = psis(log_ratios = -ll, r_eff = r_eff)
  ypredloo = E_loo(ypred, psis_object, log_ratios = -ll)$value
  ypredloo
}

epred = function(fit) {
  apply(posterior_epred(fit), 2, mean)
}

############ R2 ################
var_r2_taylor = function(y, yhat) {
  n = length(y)
  eloo = y - yhat
  mse_pred = mean(eloo^2)
  
  mean_y = mean(y)
  mse_y = mean((y - mean_y)^2)
  
  ratio = mse_pred / mse_y
  var_mse_pred = var(eloo^2) / n
  var_mse_y = var((y -mean_y)^2) / n
  cov_mse = cov(eloo^2, (y - mean_y)^2) / n
  
  
  (var_mse_pred - 2*ratio*cov_mse + ratio^2*var_mse_y) / mse_y^2
}

############ RMSE ##############

var_rmse_bb = function(y, yhat) {
  e = (y - yhat)^2
  n = length(y)
  rd=rudirichlet(4000, n)
  RMSE = sqrt(rowSums(sweep(rd, 2, e, FUN = "*")))
  var(RMSE)
}


var_rmse_diff_bb = function(y, yhat1, yhat2) {
  e1 = (y - yhat1)^2
  e2 = (y - yhat2)^2
  
  n = length(y)
  rd=rudirichlet(4000, n)
  
  RMSE1 = sqrt(rowSums(sweep(rd, 2, e1, FUN = "*")))
  RMSE2 = sqrt(rowSums(sweep(rd, 2, e2, FUN = "*")))
  var(RMSE1 - RMSE2)
}

var_rmse_taylor = function(y, yhat) {
  n = length(y)
  MSE = mean((y-yhat)^2)
  var_MSE = var((y-yhat)^2) / n 
  var_MSE / MSE / 4
}

var_rmse_diff_taylor = function(y, yhat1, yhat2) {
  n = length(y)
  MSE1 = mean((y-yhat1)^2)
  varMSE1 = var((y-yhat1)^2) / n 
  MSE2 = mean((y-yhat2)^2)
  varMSE2 = var((y-yhat2)^2) / n
  covMSE = cov((y-yhat1)^2, (y-yhat2)^2) / n
  
  varMSE1 / MSE1 / 4 + varMSE2 / MSE2 / 4 - covMSE / (2 * sqrt(MSE1)*sqrt(MSE2))
}
