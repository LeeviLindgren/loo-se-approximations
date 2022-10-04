library(bayesboot)
library(tidyverse)
library(loo)
source('experiments/util.R')
n = c(20, 50, 100, 200)
p = c(2, 3, 5)
rho = 0
sigma = 1
n_reps = 100

sim = expand_grid(n, p, rho, sigma, n_reps) %>% split(seq(nrow(.)))
