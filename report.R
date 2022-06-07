#' ---
#' title: "Analytical approximation for LOO-R2 standard error"
#' author: "Leevi Lindgren"
#' date: "`r format(Sys.Date())`"
#' output:
#'   html_document:
#'     toc: true
#'     toc_depth: 3
#'     toc_float: true
#' ---

#+ libraries, warnings = FALSE, message = FALSE
library(rstanarm)
library(loo)
library(ggplot2)
theme_set(bayesplot::theme_default(base_family = "sans"))
set.seed(123456)

#' # Taylor approximation of function of two random variables
#' Consider the general problem of estimating variance of some function 
#' $f: \mathbb{R}^2 \rightarrow \mathbb{R}$ of two random variables. 
#' 
#' Let's first make an approximation for the mean: 
#' \begin{align}
#' \text{var}(f(X, Y)) &\approx E\left[ \left( \frac{\partial f}{\partial x} (X - \mu_X) + \frac{\partial f}{\partial y} (Y - \mu_Y) \right)^2 \right] \nonumber \\
#' &= E\left[ \left( \frac{\partial f}{\partial x}  \right)^2 (X - \mu_X)^2 + 2 \frac{\partial f}{\partial x} \frac{\partial f}{\partial y} (X - \mu_X)(Y - \mu_Y) + \left( \frac{\partial f}{\partial y}  \right)^2 (Y - \mu_Y)^2  \right] \nonumber \\
#' &= \left( \frac{\partial f}{\partial x}  \right)^2 \text{var}(X) + 2 \frac{\partial f}{\partial x} \frac{\partial f}{\partial y} \text{cov}(X, Y) + \left( \frac{\partial f}{\partial y}  \right)^2 \text{var}(Y)
#' \end{align}