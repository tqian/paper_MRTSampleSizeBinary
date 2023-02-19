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

result2.6 <- readRDS("simulation2.6/result_df_collected_nsim2000.RDS")
result2.7 <- readRDS("simulation2.7/result_df_collected_nsim2000.RDS")
result2.8 <- readRDS("simulation2.8/result_df_collected_nsim2000.RDS")
result2.9 <- readRDS("simulation2.9/result_df_collected_nsim2000.RDS")
result <- rbind(result2.6, result2.7, result2.8, result2.9)
# result <- result2.6


result$ft_t_shape_name <- result$ft_t_shape
result$ft_t_shape_name <- factor(result$ft_t_shape_name)
levels(result$ft_t_shape_name) <- c("linear_theta"=TeX(r'(linear $MEE^*(t)$)'),
                                    "quadratic_theta"=TeX(r'(quadratic $MEE^*(t)$)'))
result$ft_w_shape_name <- result$ft_w_shape
result$ft_w_shape_name <- factor(result$ft_w_shape_name)
levels(result$ft_w_shape_name) <- c("linear_theta"=TeX(r'(linear $MEE^w(t)$)'),
                                    "quadratic_theta"=TeX(r'(quadratic $MEE^w(t)$)'))


result %>%
    ggplot(aes(x = ft_t_theta, y = power_adj, linetype = factor(gt_t_theta),
                      color = ft_w_theta, group = interaction(gt_t_theta, ft_w_theta))) +
    geom_line() +
    geom_hline(yintercept = 0.8, linetype = 2) +
    facet_grid(ft_t_shape_name ~ ft_w_shape_name, labeller = label_parsed) +
    scale_linetype_discrete(name = TeX(r'($SPNC^*(t) = SPNC^w(t)$)'),
                            labels = c(TeX(r'(quadratic, $\theta_g^* = \theta_g^w = -0.5$)'),
                                       TeX(r'(quadratic, $\theta_g^* = \theta_g^w = 0$)'),
                                       TeX(r'(quadratic, $\theta_g^* = \theta_g^w = 0.5$)'))) +
    scale_color_continuous(name = TeX(r'($\theta_f^w$)')) +
    xlab(TeX(r'($\theta_f^*$)')) +
    ylab("power") +
    # coord_cartesian(ylim = c(0.5, 0.83)) +
    # scale_y_continuous(breaks = seq(from = 0.1, to = 1, by = 0.1)) +
    # scale_x_continuous(breaks = c(-1,0,1)) +
    # facet_grid(gt_t_theta ~ ASPNC_t, labeller = label_both) +
    theme_bw() + 
    myggfont() +
    theme(legend.text.align = 0, plot.title = element_text(hjust = 0.5))
    
ggsave(paste0("figure A1 - power.pdf"), width = 12, height = 8)
