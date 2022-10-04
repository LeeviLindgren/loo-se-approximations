library(rstanarm)
source('experiments/config.R')


for (i in seq_along(sim)) {
  print(paste0('Iteration: ', i))
  conf = sim[[i]]
  
  # Initialize result vectors
  var_diff_bb = numeric()
  var_diff_taylor = numeric()
  var_bb = numeric()
  var_taylor = numeric()
  
  var_diff_bb_loo = numeric()
  var_diff_taylor_loo = numeric()
  var_bb_loo = numeric()
  var_taylor_loo = numeric()
  
  for (j in 1:n_reps) {
    X = get_X(conf$n, conf$p, conf$rho)
    eps = rnorm(conf$n)*conf$sigma
    beta = rnorm(conf$p)
    y = X %*% beta + eps
    df = data.frame(y, X)
    
    # Fit the models
    fit1 = stan_glm(y ~ X1, data = df, refresh=0)
    fit2 = stan_glm(y ~ X1 + X2, data = df, refresh=0)
    
    # Make mean predictions E[y_i|y] and E[y_i|y_-i]
    epred1 = epred(fit1)
    epred2 = epred(fit2)
    loo_epred1 = loo_epred(fit1)
    loo_epred2 = loo_epred(fit2)
    
    # R2
    r21 = bayes_R2(fit1)
    r22 = bayes_R2(fit2)
    loo_r21 = loo_R2(fit1)
    loo_r22 = loo_R2(fit2)
    
    # var_diff_bb[j] = var(r22 - r21)
    # var_diff_taylor[j] = IMPLEMENT
    var_bb[j] = var(r22)
    var_taylor[j] = var_r2_taylor(y, epred2)
    
    # LOO R2
    # var_diff_bb_loo[j] = var(loo_r22 - loo_r21)
    # var_diff_taylor_loo[j] = IMPLEMENT
    var_bb_loo[j] = var(loo_r22)
    var_taylor_loo[j] = var_r2_taylor(y, loo_epred2)
  }
  
  # sim[[i]]$var_r2_diff_bb = list(var_diff_bb)
  # sim[[i]]$var_r2_diff_taylor = list(var_diff_taylor)
  sim[[i]]$var_r2_bb = list(var_bb)
  sim[[i]]$var_r2_taylor = list(var_taylor)
  
  # sim[[i]]$var_loo_r2_diff_bb = list(var_diff_bb_loo)
  # sim[[i]]$var_loo_r2_diff_taylor = list(var_diff_taylor_loo)
  sim[[i]]$var_loo_r2_bb = list(var_bb_loo)
  sim[[i]]$var_loo_r2_taylor = list(var_taylor_loo)
}

sim = sim %>% 
  map(function(df) df %>% unnest_longer(starts_with('var_'))) %>%
  do.call(rbind, .) %>% 
  pivot_longer(starts_with('var_')) %>% 
  mutate(approximation = sapply(str_split(name, '_'), tail, 1),
         data_reuse = grepl('loo', name),
         data_reuse = if_else(data_reuse, 'LOO', 'In-sample'),
         metric = grepl('diff', name),
         metric = if_else(metric, 'R2 diff', 'R2')) %>%
  pivot_wider(values_from = value, 
              names_from = 'approximation', 
              id_cols = -name) %>%
  unnest_longer(c(bb, taylor))


theme_set(theme_bw(base_size = 16))
p_loo = sim %>%
  filter(data_reuse == 'LOO') %>%
  mutate(p = as.character(p),
         taylor = sqrt(taylor),
         bb = sqrt(bb)) %>%
  ggplot(aes(x = taylor, y = bb)) +
  geom_point(aes(color = p), size=2, alpha=0.5) +
  geom_abline() +
  facet_grid(rows = vars(n), cols = vars(metric)) +
  xlab('SE Taylor approximation') +
  ylab('SE Bayesian bootstrap approximation')

p_in_sample = sim %>%
  filter(data_reuse == 'In-sample') %>%
  mutate(p = as.character(p),
         taylor = sqrt(taylor),
         bb = sqrt(bb)) %>%
  ggplot(aes(x = taylor, y = bb)) +
  geom_point(aes(color = p), size=2, alpha=0.5) +
  geom_abline() +
  facet_grid(rows = vars(n), cols = vars(metric)) +
  xlab('SE Taylor approximation') +
  ylab('SE Bayesian bootstrap approximation')

ggsave('figures/loo_r2.pdf', plot = p_loo, device = 'pdf')

ggsave('figures/r2.pdf', plot = p_in_sample, device = 'pdf')
