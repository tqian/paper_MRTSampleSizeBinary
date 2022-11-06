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

result <- readRDS("simulation4.1/result_df_collected_nsim2000.RDS")

# result %>% 
#     ggplot(aes(x = ATE_tw_ratio, y = power_adj, linetype = factor(ASPNC_t),
#                color = factor(gt_t_theta))) +
#     geom_line() +
#     geom_hline(yintercept = 0.8) +
#     facet_grid(ft_t_shape ~ ft_t_theta, labeller = label_both) +
#     theme_bw() + myggfont()
# ggsave(paste0("figure 1 - power.pdf"), width = 10, height = 7.5)
# 
# result %>% 
#     ggplot(aes(x = ATE_tw_ratio, y = typeierror_adj, linetype = factor(ASPNC_t),
#                color = factor(gt_t_theta))) +
#     geom_line() +
#     geom_hline(yintercept = 0.05) +
#     facet_grid(ft_t_shape ~ ft_t_theta, labeller = label_both) +
#     theme_bw() + myggfont()
# ggsave(paste0("figure 1 - typeierror.pdf"), width = 10, height = 7.5)


result$ft_t_shape_and_theta <- paste0(result$ft_t_shape, result$ft_t_theta)
result$ASPNC_t_and_gt_t_theta <- paste0(result$ASPNC_t, result$gt_t_theta)

result %>% 
    filter(! ft_t_shape_and_theta %in% c("linear_theta0", "quadratic_theta0")) %>%
    ggplot(aes(x = ATE_tw_ratio, y = power_adj, linetype = factor(ft_t_shape_and_theta),
               color = factor(ASPNC_t_and_gt_t_theta))) +
    geom_line() +
    geom_hline(yintercept = 0.8, linetype = 2) +
    scale_color_discrete(name = TeX(r'(Quadratic $SPNC^*(t) = SPNC^w(t)$)'),
                         labels = c(TeX(r'($ASPN^* = ASPN^w = 0.2$, $\theta_g^* = \theta_g^w = -0.3$)'),
                                    TeX(r'($ASPN^* = ASPN^w = 0.2$, $\theta_g^* = \theta_g^w = 0$)'),
                                    TeX(r'($ASPN^* = ASPN^w = 0.2$, $\theta_g^* = \theta_g^w = 0.3$)'),
                                    TeX(r'($ASPN^* = ASPN^w = 0.3$, $\theta_g^* = \theta_g^w = -0.3$)'),
                                    TeX(r'($ASPN^* = ASPN^w = 0.3$, $\theta_g^* = \theta_g^w = 0$)'),
                                    TeX(r'($ASPN^* = ASPN^w = 0.3$, $\theta_g^* = \theta_g^w = 0.3$)'))) +
    scale_linetype_discrete(name = TeX(r'(Pattern of $MEE^*(t)$ and $MEE^w(t)$ )'),
        labels = c(TeX(r'(constant)'),
                   TeX(r'(linear, $\theta_f^* = \theta_f^w = -0.3$)'),
                   # TeX(r'(linear, $\theta_f^* = \theta_f^w = 0$)'),
                   TeX(r'(linear, $\theta_f^* = \theta_f^w = 0.3$)'),
                   TeX(r'(quadratic, $\theta_f^* = \theta_f^w = -0.3$)'),
                   # TeX(r'(quadratic, $\theta_f^* = \theta_f^w = 0$)'),
                   TeX(r'(quadratic, $\theta_f^* = \theta_f^w = 0.3$)'))) +
    xlab(TeX(r'($ATE^*/ATE^w$)')) +
    ylab("power") +
    scale_y_continuous(breaks = c(0.2, 0.4, 0.6, 0.8, 1)) +
    # facet_grid(gt_t_theta ~ ASPNC_t, labeller = label_both) +
    theme_bw() + 
    myggfont() +
    theme(legend.text.align = 0)
ggsave(paste0("figure 1 - power.pdf"), width = 13, height = 8)



result %>% 
    filter(! ft_t_shape_and_theta %in% c("linear_theta0", "quadratic_theta0")) %>%
    ggplot(aes(x = ATE_tw_ratio, y = typeierror_adj, linetype = factor(ft_t_shape_and_theta),
               color = factor(ASPNC_t_and_gt_t_theta))) +
    geom_line() +
    geom_hline(yintercept = 0.05, linetype = 2) +
    scale_color_discrete(name = TeX(r'(Quadratic $SPNC^*(t) = SPNC^w(t)$)'),
                            labels = c(TeX(r'($ASPN^* = ASPN^w = 0.2$, $\theta_g^* = \theta_g^w = -0.3$)'),
                                       TeX(r'($ASPN^* = ASPN^w = 0.2$, $\theta_g^* = \theta_g^w = 0$)'),
                                       TeX(r'($ASPN^* = ASPN^w = 0.2$, $\theta_g^* = \theta_g^w = 0.3$)'),
                                       TeX(r'($ASPN^* = ASPN^w = 0.3$, $\theta_g^* = \theta_g^w = -0.3$)'),
                                       TeX(r'($ASPN^* = ASPN^w = 0.3$, $\theta_g^* = \theta_g^w = 0$)'),
                                       TeX(r'($ASPN^* = ASPN^w = 0.3$, $\theta_g^* = \theta_g^w = 0.3$)'))) +
    scale_linetype_discrete(name = TeX(r'(Pattern of $MEE^*(t)$ and $MEE^w(t)$ )'),
                         labels = c(TeX(r'(constant)'),
                                    TeX(r'(linear, $\theta_f^* = \theta_f^w = -0.3$)'),
                                    # TeX(r'(linear, $\theta_f^* = \theta_f^w = 0$)'),
                                    TeX(r'(linear, $\theta_f^* = \theta_f^w = 0.3$)'),
                                    TeX(r'(quadratic, $\theta_f^* = \theta_f^w = -0.3$)'),
                                    # TeX(r'(quadratic, $\theta_f^* = \theta_f^w = 0$)'),
                                    TeX(r'(quadratic, $\theta_f^* = \theta_f^w = 0.3$)'))) +
    xlab(TeX(r'($ATE^*/ATE^w$)')) +
    ylab("type I error") +
    scale_y_continuous(breaks = c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06)) +
    # facet_grid(gt_t_theta ~ ASPNC_t, labeller = label_both) +
    theme_bw() + 
    myggfont() +
    theme(legend.text.align = 0)
ggsave(paste0("figure 1 - typeierror.pdf"), width = 13, height = 8)
