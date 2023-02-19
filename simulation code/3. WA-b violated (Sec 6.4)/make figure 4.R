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


result <- readRDS("simulation5.2/result_df_collected_nsim2000.RDS")



result$gt_t_shape_and_theta <- paste0(result$gt_t_shape, result$gt_t_theta)
# result$ASPNC_t_and_gt_t_theta <- paste0(result$ASPNC_t, result$gt_t_theta)

result$ASPNC_tw_diff <- result$ASPNC_t - result$ASPNC_w
# seems to solely depend on the difference, not the ratio

result %>% 
    filter(! gt_t_shape_and_theta %in% c("linear_theta0", "quadratic_theta0")) %>%
    ggplot(aes(x = ASPNC_tw_ratio, y = power_adj, linetype = factor(ASPNC_t),
               color = factor(gt_t_shape_and_theta))) +
    geom_line() +
    geom_hline(yintercept = 0.8, linetype = 2) +
    scale_linetype_discrete(name = TeX(r'(Value of $ASPN^*$)')) +
    scale_color_discrete(name = TeX(r'(Pattern of $SPNC^*(t)$ and $SPNC^w(t)$ )'),
                            labels = c(TeX(r'(constant)'),
                                       TeX(r'(linear, $\theta_g^* = \theta_g^w = -0.3$)'),
                                       # TeX(r'(linear, $\theta_g^* = \theta_g^w = 0$)'),
                                       TeX(r'(linear, $\theta_g^* = \theta_g^w = 0.3$)'),
                                       TeX(r'(quadratic, $\theta_g^* = \theta_g^w = -0.3$)'),
                                       # TeX(r'(quadratic, $\theta_g^* = \theta_g^w = 0$)'),
                                       TeX(r'(quadratic, $\theta_g^* = \theta_g^w = 0.3$)'))) +
    xlab(TeX(r'($ASPN^*/ASPN^w$)')) +
    ylab("power") +
    scale_y_continuous(breaks = (1:10)/10) +
    # facet_grid(gt_t_theta ~ ASPNC_t, labeller = label_both) +
    theme_bw() + 
    myggfont() +
    theme(legend.text.align = 0)
ggsave(paste0("figure 4 - power.pdf"), width = 12, height = 8)

# result %>% 
#     ggplot(aes(x = ASPNC_tw_ratio, y = power_adj, color = factor(ASPNC_t))) +
#     geom_line() +
#     geom_hline(yintercept = 0.8) +
#     facet_grid(gt_t_shape ~ gt_t_theta, labeller = label_both) +
#     theme_bw() + myggfont()
# ggsave("figure 4 - power.pdf", width = 10, height = 7.5)
# 
# result %>% 
#     ggplot(aes(x = ASPNC_tw_ratio, y = typeierror_adj, color = factor(ASPNC_t))) +
#     geom_line() +
#     geom_hline(yintercept = 0.05) +
#     ggtitle(paste0("nsim = ", nsim)) +
#     facet_grid(gt_t_shape ~ gt_t_theta, labeller = label_both) +
#     theme_bw() + myggfont()
# ggsave("figure 4 - typeierror.pdf", width = 10, height = 7.5)
