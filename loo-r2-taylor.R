library(loo)
library(rstanarm)
library(foreach)


n_cores <- as.integer(Sys.getenv('SLURM_CPUS_PER_TASK', parallel::detectCores() - 2))
cl <- parallel::makeCluster(n_cores)
doParallel::registerDoParallel(cl)


n <- c(100, 500, 1000)
p <- c(1, 5, 10)
sigma <- c(2, 1, 0.5)
K <- c(100)

sim_params <- expand.grid(n=n, p=p, sigma=sigma, K=K)
sim_params <- split(sim_params, seq_len(nrow(sim_params)))


system.time({results <- foreach(sim_param = sim_params,
                               .packages = c('rstanarm', 'loo')) %dopar% 
  {
    var_taylor <- numeric()
    var_bb <- numeric()
    r2_means <- numeric()
    for (k in 1:sim_param$K) {
      
      X <- matrix(rnorm(sim_param$n * sim_param$p), nrow = sim_param$n)
      beta <- rnorm(sim_param$p)
      eps <- rnorm(sim_param$n) * sim_param$sigma
      y <- X %*% beta + eps
      data <- data.frame(y = y, X)
      
      fit <- stan_glm(y ~ ., data = data, refresh = 0)
      loo_r2 <- loo_R2(fit)
      
      
      ypred <- posterior_epred(fit)
      ll <- log_lik(fit)
      M <- length(fit$stanfit@sim$n_save)
      S <- fit$stanfit@sim$n_save[[1]]
      r_eff <- relative_eff(exp(ll), chain_id = rep(1:M, each = S))
      psis_object <- psis(log_ratios = -ll, r_eff = r_eff)
      ypredloo <- E_loo(ypred, psis_object, log_ratios = -ll)$value
      
      eloo <- y - ypredloo
      mse_pred <- mean(eloo^2)
      
      mean_y <- mean(y)
      mse_y <- mean((y - mean_y)^2)
      
      ratio <- mse_pred / mse_y
      var_mse_pred <- var(eloo^2) / sim_param$n
      var_mse_y <- var((y -mean_y)^2) / sim_param$n
      cov_mse <- cov(eloo^2, (y - mean_y)^2) / sim_param$n
      
      
      var_taylor[k] <- (var_mse_pred - 2*ratio*cov_mse + ratio^2*var_mse_y) / mse_y^2
      var_bb[k] <- var(loo_r2)
      r2_means[k] <- mean(loo_r2)
    }
    result <- list()
    result$sim_params <- sim_param
    result$df <- data.frame(var_taylor, var_bb, R2 = r2_means, sd_taylor = sqrt(var_taylor), sd_bb = sqrt(var_bb))
    result
  }
})


saveRDS(results, 'simres.Rds')
