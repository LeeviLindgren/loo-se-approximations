library(tidyverse)


results <- readRDS('simres_high_p_n50.Rds')

tibble(results) %>% 
  unnest_wider(results) %>%
  unnest(c(sim_params, df)) %>%
  mutate(facet_label = paste0('n = ', n, ' sigma = ', sigma, ' p = ',  p)) %>%
  arrange(n, sigma, p) %>%
  ggplot(aes(x = sd_taylor, y = sd_bb, color = R2)) +
  geom_abline() +
  geom_point() +
  xlab('standard error taylor approximation') +
  ylab('standard error bayesian bootstrap') +
  facet_wrap(~ facet_label) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

w = 280
#ggsave('figures/simres_high_p_n50.png', width = w, height = w * 1.414, units = "mm")
ggsave('figures/simres_high_p_n50.png', width = w, height = w * 0.7, units = "mm", dpi = 150)

