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

result6.6 <- readRDS("simulation6.6/result_df_collected_nsim2000.RDS")
result6.7 <- readRDS("simulation6.7/result_df_collected_nsim2000.RDS")
result6.8 <- readRDS("simulation6.8/result_df_collected_nsim2000.RDS")
result6.9 <- readRDS("simulation6.9/result_df_collected_nsim2000.RDS")
result <- rbind(result6.6, result6.7, result6.8, result6.9)
# result <- result6.6


result$gt_t_shape_name <- result$gt_t_shape
result$gt_t_shape_name <- factor(result$gt_t_shape_name)
levels(result$gt_t_shape_name) <- c("linear_theta"=TeX(r'(linear $SPNC^*(t)$)'),
                                    "quadratic_theta"=TeX(r'(quadratic $SPNC^*(t)$)'))
result$gt_w_shape_name <- result$gt_w_shape
result$gt_w_shape_name <- factor(result$gt_w_shape_name)
levels(result$gt_w_shape_name) <- c("linear_theta"=TeX(r'(linear $SPNC^w(t)$)'),
                                    "quadratic_theta"=TeX(r'(quadratic $SPNC^w(t)$)'))

result$ft_t_shape_and_ft_t_theta <- paste0(result$ft_t_shape, result$ft_t_theta)

result %>%
    filter(!(ft_t_shape_and_ft_t_theta %in% c("linear_theta0", "quadratic_theta0"))) %>%
    ggplot(aes(x = gt_t_theta, y = power_adj, linetype = factor(ft_t_shape_and_ft_t_theta),
                      color = gt_w_theta, group = interaction(gt_w_theta, ft_t_shape_and_ft_t_theta))) +
    geom_line() +
    geom_hline(yintercept = 0.8, linetype = 2) +
    facet_grid(gt_t_shape_name ~ gt_w_shape_name, labeller = label_parsed) +
    scale_linetype_discrete(name = TeX(r'($MEE^*(t) = MEE^w(t)$)'),
                            labels = c("constant",
                                TeX(r'(linear, $\theta_f^* = \theta_f^w = -0.5$)'),
                                TeX(r'(linear, $\theta_f^* = \theta_f^w = 0.5$)'),
                                TeX(r'(quadratic, $\theta_f^* = \theta_f^w = -0.5$)'),
                                TeX(r'(quadratic, $\theta_f^* = \theta_f^w = 0.5$)'))) +
    scale_color_continuous(name = TeX(r'($\theta_g^w$)')) +
    xlab(TeX(r'($\theta_g^*$)')) +
    ylab("power") +
    # coord_cartesian(ylim = c(0.5, 0.83)) +
    # scale_y_continuous(breaks = seq(from = 0.1, to = 1, by = 0.1)) +
    # scale_x_continuous(breaks = c(-1,0,1)) +
    # facet_grid(gt_t_theta ~ ASPNC_t, labeller = label_both) +
    theme_bw() + 
    myggfont() +
    theme(legend.text.align = 0, plot.title = element_text(hjust = 0.5))
ggsave(paste0("figure A2 - power.pdf"), width = 12, height = 8)
