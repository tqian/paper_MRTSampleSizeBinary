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

result11.6 <- readRDS("simulation11.6/result_df_collected_nsim2000.RDS")
result11.6$correct_tau <- "2. correct form"
result11.7 <- readRDS("simulation11.7/result_df_collected_nsim2000.RDS")
result11.7$correct_tau <- "1. incorrect form and incorrect AA"
result11.8 <- readRDS("simulation11.8/result_df_collected_nsim2000.RDS")
result11.8$correct_tau <- "3. incorrect form but correct AA"
result <- rbind(result11.6, result11.7, result11.8)

# result <- result %>% filter(!(gt_t_shape != "constant" & gt_t_theta == 0))


# result$ft_t_shape_and_theta <- paste0(result$ft_t_shape, result$ft_t_theta)
# result$ASPNC_t_and_gt_shape_and_gt_t_theta <- paste0(result$ASPNC_t, result$gt_t_shape, result$gt_t_theta)
# result$gt_shape_and_gt_t_theta <- paste0(result$gt_t_shape, result$gt_t_theta)
# result$ATE_t_and_ASPNC_t <- paste0(result$ATE_t, result$ASPNC_t)

# result$gt_shape_and_gt_t_theta[result$gt_shape_and_gt_t_theta == "linear_theta0"] <- "constant"

# gt_t_shape_names <- list(
#     'linear_theta'=TeX(r'($SPNC^*(t)$: linear)'),
#     'quadratic_theta'=TeX(r'($SPNC^*(t)$: quadratic)')
# )
# gt_t_shape_labeller <- function(variable,value){
#     return(gt_t_shape_names[value])
# }

# result$gt_t_shape_theta <- paste0(result$gt_t_shape, "_", result$gt_t_theta)

result$gamma3_plus_gamma4 <- result$gamma3 + result$gamma4

result %>%
    ggplot(aes(x = gamma3_plus_gamma4, y = power_adj, color = gamma3,
               group = factor(gamma3))) +
    geom_line() +
    facet_grid( ~ correct_tau)
ggsave(paste0("figure 11 - explore1.pdf"), width = 9, height = 5)


result %>%
    ggplot(aes(x = gamma3_plus_gamma4, y = power_adj, color = gamma4,
               group = factor(gamma4))) +
    geom_line() +
    facet_grid( ~ correct_tau)
ggsave(paste0("figure 11 - explore2.pdf"), width = 9, height = 5)

library(metR)
breaks <- seq(from = 0.7, to = 0.9, by = 0.02)
result %>%
    ggplot(aes(x = gamma3, y = gamma4, z = power_adj)) +
    geom_contour_fill(breaks = breaks) +
    geom_text_contour(breaks = breaks, stroke = 0.1, size = 3, skip = 0) +
    facet_grid( ~ correct_tau)
ggsave(paste0("figure 11 - explore3 (contour).pdf"), width = 9, height = 5)


result %>%
    ggplot(aes(x = gamma4, y = power_adj,
               color = gamma3, group = gamma3)) +
    geom_line() +
    facet_grid( ~ correct_tau) +
    # facet_grid(~ gt_t_shape, labeller = gt_t_shape_labeller) +
    geom_hline(yintercept = 0.8, linetype = 2) +
    # scale_color_discrete(guide = "none") +
    scale_color_continuous(name = TeX(r'($\gamma_3^*$)')) +
    xlab(TeX(r'($\gamma_4^*$)')) +
    ylab("power") +
    # ggtitle(TeX(r'($SPNC^w(t)$: constant)')) +
    # coord_cartesian(ylim = c(0.5, 0.83), xlim = c(-1, 1)) +
    # scale_y_continuous(breaks = seq(from = 0.1, to = 1, by = 0.1)) +
    # scale_x_continuous(breaks = c(-1,0,1)) +
    # facet_grid(gt_t_theta ~ ASPNC_t, labeller = label_both) +
    theme_bw()
ggsave(paste0("figure 11 - explore4.pdf"), width = 9, height = 5)


result %>%
    ggplot(aes(x = gamma3, y = power_adj,
               color = gamma4, group = gamma4)) +
    geom_line() +
    facet_grid( ~ correct_tau) +
    # facet_grid(~ gt_t_shape, labeller = gt_t_shape_labeller) +
    geom_hline(yintercept = 0.8, linetype = 2) +
    # scale_color_discrete(guide = "none") +
    scale_color_continuous(name = TeX(r'($\gamma_4^*$)')) +
    xlab(TeX(r'($\gamma_3^*$)')) +
    ylab("power") +
    # ggtitle(TeX(r'($SPNC^w(t)$: constant)')) +
    # coord_cartesian(ylim = c(0.5, 0.83), xlim = c(-1, 1)) +
    # scale_y_continuous(breaks = seq(from = 0.1, to = 1, by = 0.1)) +
    # scale_x_continuous(breaks = c(-1,0,1)) +
    # facet_grid(gt_t_theta ~ ASPNC_t, labeller = label_both) +
    theme_bw()
ggsave(paste0("figure 11 - explore5.pdf"), width = 9, height = 5)


result$correct_tau_name <- result$correct_tau
result$correct_tau_name <- factor(result$correct_tau_name)

levels(result$correct_tau_name) <- c("1. incorrect form and incorrect AA"=TeX(r'($\tau^w(t) \neq \tau^*(t)$, $AA^w \neq AA^*$)'),
                                     "2. correct form"=TeX(r'($\tau^w(t) = \tau^*(t)$, $AA^w = AA^*$)'),
                                     "3. incorrect form but correct AA"=TeX(r'($\tau^w(t) \neq \tau^*(t)$, $AA^w = AA^*$)'))

result %>%
    filter(correct_tau != "3. incorrect form but correct AA") %>%
    ggplot(aes(x = gamma4, y = power_adj,
               color = gamma3, group = gamma3)) +
    geom_line() +
    facet_grid( ~ correct_tau_name, labeller = label_parsed) +
    # facet_grid(~ gt_t_shape, labeller = gt_t_shape_labeller) +
    geom_hline(yintercept = 0.8, linetype = 2) +
    # scale_color_discrete(guide = "none") +
    scale_color_continuous(name = TeX(r'($\gamma_3^*$)')) +
    xlab(TeX(r'($\gamma_4^*$)')) +
    ylab("power") +
    # ggtitle(TeX(r'($SPNC^w(t)$: constant)')) +
    # coord_cartesian(ylim = c(0.5, 0.83), xlim = c(-1, 1)) +
    # scale_y_continuous(breaks = seq(from = 0.1, to = 1, by = 0.1)) +
    # scale_x_continuous(breaks = c(-1,0,1)) +
    # facet_grid(gt_t_theta ~ ASPNC_t, labeller = label_both) +
    theme_bw() + 
    myggfont(legend_text_size = 14, legend_title_size = 16) +
    theme(legend.text.align = 0, plot.title = element_text(hjust = 0.5))
ggsave(paste0("figure 11 - power.pdf"), width = 9, height = 5)

# result %>%
#     ggplot(aes(x = gamma3, y = power_adj,
#                color = gamma4, group = gamma4)) +
#     geom_line() +
#     facet_grid( ~ correct_tau) +
#     # facet_grid(~ gt_t_shape, labeller = gt_t_shape_labeller) +
#     geom_hline(yintercept = 0.8, linetype = 2) +
#     # scale_color_discrete(guide = "none") +
#     scale_color_continuous(name = TeX(r'($\gamma_4^*$)')) +
#     xlab(TeX(r'($\gamma_3^*$)')) +
#     ylab("power") +
#     # ggtitle(TeX(r'($SPNC^w(t)$: constant)')) +
#     # coord_cartesian(ylim = c(0.5, 0.83), xlim = c(-1, 1)) +
#     # scale_y_continuous(breaks = seq(from = 0.1, to = 1, by = 0.1)) +
#     # scale_x_continuous(breaks = c(-1,0,1)) +
#     # facet_grid(gt_t_theta ~ ASPNC_t, labeller = label_both) +
#     theme_bw() + 
#     myggfont(legend_text_size = 14, legend_title_size = 16) +
#     theme(legend.text.align = 0, plot.title = element_text(hjust = 0.5))
# ggsave(paste0("figure 11 - power 2.pdf"), width = 9, height = 5)
