### plot making ###
# Binary MRT sample size paper
#
# see _plot making log.md for simulation set up
#
# Tianchen Qian
# 2022.06.20


rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(tidyverse)
library(latex2exp)
library(latex2exp)
source("../../functions/myggplot.R")

result <- readRDS("simulation9.1/result_df_collected_nsim2000.RDS")

result$gt_shape_and_gt_t_theta <- paste0(result$gt_t_shape, result$gt_t_theta)
result <- result %>% filter(! gt_shape_and_gt_t_theta %in% c("linear_theta0", "quadratic_theta0"))

result$ft_shape_and_ft_t_theta <- paste0(result$ft_t_shape, result$ft_t_theta)
result <- result %>% filter(! ft_shape_and_ft_t_theta %in% c("linear_theta0", "quadratic_theta0"))

result$tau_t_shape_and_tau_t_theta <- paste0(result$tau_t_shape, result$tau_t_theta)

result <- result %>% filter(! tau_t_shape_and_tau_t_theta %in% c("linear0", "sine0"))

result$tau_t_theta[result$tau_t_shape == "constant"] <- 99

result$tau_t_shape_name <- result$tau_t_shape
result$tau_t_shape_name[result$tau_t_shape_name == "sine"] <- "periodic"
result$tau_t_shape_name <- factor(result$tau_t_shape_name)
result$tau_t_theta_name <- factor(result$tau_t_theta)

levels(result$tau_t_shape_name) <- c("constant"=TeX(r'($\tau(t)$: constant)'),
                                     "periodic"=TeX(r'($\tau(t)$: periodic)'),
                                     "linear"=TeX(r'($\tau(t)$: linear)'))
levels(result$tau_t_theta_name) <- c("-0.2" = TeX(r'($\theta^*_\tau = \theta^\w_\tau = -0.2$)'),
                                     "0.2" = TeX(r'($\theta^*_\tau = \theta^\w_\tau = 0.2$)'),
                                     "99" = " ")

# tau_t_shape_and_tau_t_theta_names <- list(
#     'constant0'=TeX(r'(constant)'),
#     'linear-0.2'=TeX(r'(linear, $\theta^*_\tau = \theta^\w_\tau = -0.2$)'),
#     'linear0.2'=TeX(r'(linear, $\theta^*_\tau = \theta^\w_\tau = 0.2$)'),
#     'sine-0.2'=TeX(r'(periodic, $\theta^*_\tau = \theta^\w_\tau = -0.2$)'),
#     'sine0.2'=TeX(r'(periodic, $\theta^*_\tau = \theta^\w_\tau = 0.2$)')
# )
# tau_t_shape_and_tau_t_theta_labeller <- function(variable,value){
#     return(tau_t_shape_and_tau_t_theta_names[value])
# }
# 
# tau_t_shape_names <- list(
#     'constant'="constant",
#     'linear'="linear",
#     'sine'="periodic")
# tau_t_shape_labeller <- function(variable,value){
#     return(tau_t_shape_names[value])
# }
# 
# tau_t_theta_names <- list(
#     '-0.2'=TeX(r'($\theta^*_\tau = \theta^\w_\tau = -0.2$)'),
#     '0.2'=TeX(r'($\theta^*_\tau = \theta^\w_\tau = 0.2$)'),
#     '99'="")
# tau_t_theta_labeller <- function(variable,value){
#     return(tau_t_theta_names[value])
# }
# 
# tau_t_shape_and_tau_t_theta_labeller <- function(variable,value){
#     return(tau_t_shape_and_tau_t_theta_names[value])
# }

result %>% 
    ggplot(aes(x = AvgTau_tw_ratio, y = power_adj, linetype = factor(gt_shape_and_gt_t_theta),
               color = factor(ft_shape_and_ft_t_theta))) +
    geom_line() +
    # facet_grid(~ tau_t_shape_and_tau_t_theta, labeller = tau_t_shape_and_tau_t_theta_labeller) +
    facet_grid(~ tau_t_shape_name + tau_t_theta_name, labeller=label_parsed) +
    geom_hline(yintercept = 0.8, linetype = 2) +
    scale_linetype_discrete(name = TeX(r'(Pattern of $SPNC^*(t) = SPNC^w(t)$)'),
                            labels = c(TeX(r'(constant)'),
                                       TeX(r'(linear, $\theta_g^* = \theta_g^w = -0.3$)'),
                                       TeX(r'(linear, $\theta_g^* = \theta_g^w = 0.3$)'),
                                       TeX(r'(quadratic, $\theta_g^* = \theta_g^w = -0.3$)'),
                                       TeX(r'(quadratic, $\theta_g^* = \theta_g^w = 0.3$)'))) +
    scale_color_discrete(name = TeX(r'(Pattern of $MEE^*(t) = MEE^w(t)$)'),
                            labels = c(TeX(r'(constant)'),
                                       TeX(r'(linear, $\theta_f^* = \theta_f^w = -0.3$)'),
                                       TeX(r'(linear, $\theta_f^* = \theta_f^w = 0.3$)'),
                                       TeX(r'(quadratic, $\theta_f^* = \theta_f^w = -0.3$)'),
                                       TeX(r'(quadratic, $\theta_f^* = \theta_f^w = 0.3$)'))) +
    xlab(TeX(r'($AA^*/AA^w$)')) +
    ylab("power") +
    coord_cartesian(xlim = c(0.78, 1.22)) +
    # scale_y_continuous(breaks = seq(from = 0.1, to = 1, by = 0.1)) +
    scale_x_continuous(breaks = c(0.8, 0.9, 1.0, 1.1, 1.2)) +
    # facet_grid(gt_t_theta ~ ASPNC_t, labeller = label_both) +
    theme_bw() + 
    myggfont() +
    theme(legend.text.align = 0, plot.title = element_text(hjust = 0.5))
ggsave(paste0("figure 8 - power.pdf"), width = 15, height = 6)

# result %>% 
#     ggplot(aes(x = AvgTau_tw_ratio, y = power_adj, color = tau_t_shape,
#                linetype = factor(tau_t_theta))) +
#     geom_line() +
#     geom_hline(yintercept = 0.8) +
#     facet_grid(ft_t_shape + ft_t_theta ~ gt_t_shape + gt_t_theta,
#                labeller = label_both) +
#     theme_bw() + myggfont()
# ggsave(paste0("figure 8 - power.pdf"), width = 10, height = 7.5)
# 
# result %>% 
#     ggplot(aes(x = AvgTau_tw_ratio, y = typeierror_adj, color = tau_t_shape,
#                linetype = factor(tau_t_theta))) +
#     geom_line() +
#     geom_hline(yintercept = 0.05) +
#     facet_grid(ft_t_shape + ft_t_theta ~ gt_t_shape + gt_t_theta,
#                labeller = label_both) +
#     theme_bw() + myggfont()
# ggsave(paste0("figure 8 - typeierror.pdf"), width = 10, height = 7.5)