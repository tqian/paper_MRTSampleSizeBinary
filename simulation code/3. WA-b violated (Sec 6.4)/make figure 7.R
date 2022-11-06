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

result7.1 <- readRDS("simulation7.1/result_df_collected_nsim2000.RDS")
result7.1$correct_ASPN <- "correct ASPN"
result7.2 <- readRDS("simulation7.2/result_df_collected_nsim2000.RDS")
result7.2$correct_ASPN <- "incorrect ASPN"
result <- rbind(result7.1, result7.2)


result$ft_t_shape_and_ft_t_theta <- paste0(result$ft_t_shape, result$ft_t_theta)
result$gt_t_shape_and_gt_t_theta <- paste0(result$gt_t_shape, result$gt_t_theta)


result$correct_ASPN_name <- result$correct_ASPN
result$correct_ASPN_name <- factor(result$correct_ASPN_name)

levels(result$correct_ASPN_name) <- c("correct ASPN"=TeX(r'($ASPN^w = APSN^*$)'),
                                      "incorrect ASPN"=TeX(r'($ASPN^w \neq APSN^*$)'))


result %>% 
    filter(! ft_t_shape_and_ft_t_theta %in% c("linear_theta0", "quadratic_theta0")) %>%
    filter(! gt_t_shape_and_gt_t_theta %in% c("linear_theta0", "quadratic_theta0")) %>%
    ggplot(aes(x = coef_DE, y = power_adj, linetype = factor(ft_t_shape_and_ft_t_theta),
               color = factor(gt_t_shape_and_gt_t_theta))) +
    geom_line() +
    facet_grid( ~ correct_ASPN_name, labeller = label_parsed) +
    geom_hline(yintercept = 0.8, linetype = 2) +
    scale_color_discrete(name = TeX(r'(Pattern of $SPNC^*(t)$ and $SPNC^w(t)$)'),
                         labels = c(TeX(r'(constant)'),
                                    TeX(r'(linear, $\theta_g^* = \theta_g^w = -0.3$)'),
                                    TeX(r'(linear, $\theta_g^* = \theta_g^w = 0.3$)'),
                                    TeX(r'(quadratic, $\theta_g^* = \theta_g^w = -0.3$)'),
                                    TeX(r'(quadratic, $\theta_g^* = \theta_g^w = 0.3$)'))) +
    scale_linetype_discrete(name = TeX(r'(Pattern of $MEE^*(t) = MEE^w(t)$ )'),
                            labels = c(TeX(r'(constant)'),
                                       TeX(r'(linear, $\theta_f^* = \theta_f^w = -0.3$)'),
                                       TeX(r'(linear, $\theta_f^* = \theta_f^w = 0.3$)'),
                                       TeX(r'(quadratic, $\theta_f^* = \theta_f^w = -0.3$)'),
                                       TeX(r'(quadratic, $\theta_f^* = \theta_f^w = 0.3$)'))) +
    xlab(TeX(r'($\gamma_1^*$ (delayed effect))')) +
    ylab("power") +
    # scale_y_continuous(breaks = (1:10)/10) +
    # coord_cartesian(ylim = c(0.7, 0.85)) +
    # facet_grid(gt_t_theta ~ ASPNC_t, labeller = label_both) +
    theme_bw() + 
    myggfont() +
    theme(legend.text.align = 0)
ggsave(paste0("figure 7 - power.pdf"), width = 13, height = 6)


# result %>% 
#     # pivot_longer(c("power_adj", "power_unadj"),
#     # names_to = "ss_correction", values_to = "power") %>%
#     ggplot(aes(x = coef_DE, y = power_adj,
#                color = ft_t_shape, linetype = gt_t_shape)) +
#     # geom_point() +
#     geom_line() +
#     geom_hline(yintercept = 0.8) +
#     facet_grid(ft_t_theta ~ gt_t_theta,
#                labeller = label_both) +
#     theme_bw() + myggfont()
# ggsave(paste0("figure 7 - power.pdf"), width = 10, height = 7.5)
# 
# result %>% 
#     # pivot_longer(c("power_adj", "power_unadj"),
#     # names_to = "ss_correction", values_to = "power") %>%
#     ggplot(aes(x = coef_DE, y = typeierror_adj,
#                color = ft_t_shape, linetype = gt_t_shape)) +
#     # geom_point() +
#     geom_line() +
#     geom_hline(yintercept = 0.05) +
#     facet_grid(ft_t_theta ~ gt_t_theta,
#                labeller = label_both) +
#     theme_bw() + myggfont()
# ggsave(paste0("figure 7 - typeierror.pdf"), width = 10, height = 7.5)