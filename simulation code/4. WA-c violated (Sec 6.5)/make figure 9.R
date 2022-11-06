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
source("../../functions/myggplot.R")

result <- readRDS("simulation8.1/result_df_collected_nsim2000.RDS")


result$gt_shape_and_gt_t_theta <- paste0(result$gt_t_shape, result$gt_t_theta)
result <- result %>% filter(! gt_shape_and_gt_t_theta %in% c("linear_theta0", "quadratic_theta0"))

result$ft_shape_and_ft_t_theta <- paste0(result$ft_t_shape, result$ft_t_theta)
result <- result %>% filter(! ft_shape_and_ft_t_theta %in% c("linear_theta0", "quadratic_theta0"))

result$tau_t_shape_and_tau_t_theta <- paste0(result$tau_t_shape, result$tau_t_theta)

result <- result %>% filter(! tau_t_shape_and_tau_t_theta %in% c("linear0", "sine0"))

result <- result %>% filter(tau_t_shape != "constant")

result$tau_t_shape_name <- result$tau_t_shape
result$tau_t_shape_name[result$tau_t_shape_name == "sine"] <- "periodic"
result$tau_t_shape_name <- factor(result$tau_t_shape_name)

levels(result$tau_t_shape_name) <- c("periodic"=TeX(r'($\tau^*(t)$: periodic)'),
                                     "linear"=TeX(r'($\tau^*(t)$: linear)'))


result %>% 
    ggplot(aes(x = tau_t_theta, y = power_adj, linetype = factor(gt_shape_and_gt_t_theta),
               color = factor(ft_shape_and_ft_t_theta))) +
    geom_line() +
    # facet_grid(~ tau_t_shape_and_tau_t_theta, labeller = tau_t_shape_and_tau_t_theta_labeller) +
    facet_grid(~ tau_t_shape_name, labeller = label_parsed) +
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
    xlab(TeX(r'($\theta_\tau^*$)')) +
    ylab("power") +
    coord_cartesian(ylim = c(0.7, 0.85), xlim = c(-0.43, 0.43)) +
    scale_x_continuous(breaks = c(-0.4, -0.2, 0, 0.2, 0.4)) +
    ggtitle(TeX(r'($\tau^w(t)$: constant)')) +
    theme_bw() + 
    myggfont() +
    theme(legend.text.align = 0, plot.title = element_text(hjust = 0.5))
ggsave(paste0("figure 9 - power.pdf"), width = 10, height = 6)

# result %>% 
#     ggplot(aes(x = tau_t_theta, y = power_adj, color = tau_t_shape)) +
#     geom_line(data = result %>% filter(tau_t_shape != "constant")) +
#     geom_point(data = result %>% filter(tau_t_shape == "constant")) +
#     geom_hline(yintercept = 0.8) +
#     facet_grid(ft_t_shape + ft_t_theta ~ gt_t_shape + gt_t_theta,
#                labeller = label_both) +
#     theme_bw() + myggfont()
# ggsave(paste0("figure 9 - power.pdf"), width = 10, height = 7.5)
# 
# result %>% 
#     ggplot(aes(x = tau_t_theta, y = typeierror_adj, color = tau_t_shape)) +
#     geom_line(data = result %>% filter(tau_t_shape != "constant")) +
#     geom_point(data = result %>% filter(tau_t_shape == "constant")) +
#     geom_hline(yintercept = 0.05) +
#     facet_grid(ft_t_shape + ft_t_theta ~ gt_t_shape + gt_t_theta,
#                labeller = label_both) +
#     theme_bw() + myggfont()
# ggsave(paste0("figure 9 - typeierror.pdf"), width = 10, height = 7.5)