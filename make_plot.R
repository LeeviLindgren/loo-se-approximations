library(tidyverse)


results <- readRDS('simres.Rds')

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

ggsave('figures/simres.png', width = 15, height = 15, dpi = 300)


for (res in results) {
  
}