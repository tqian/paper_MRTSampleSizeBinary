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

result6.4 <- readRDS("simulation6.4/result_df_collected_nsim2000.RDS")
result6.5 <- readRDS("simulation6.5/result_df_collected_nsim2000.RDS")
result <- rbind(result6.4, result6.5)

# result$ft_t_shape_and_theta <- paste0(result$ft_t_shape, result$ft_t_theta)
# result$ASPNC_t_and_gt_shape_and_gt_t_theta <- paste0(result$ASPNC_t, result$gt_t_shape, result$gt_t_theta)
# result$gt_shape_and_gt_t_theta <- paste0(result$gt_t_shape, result$gt_t_theta)
# result$ATE_t_and_ASPNC_t <- paste0(result$ATE_t, result$ASPNC_t)
# 
# result$gt_shape_and_gt_t_theta[result$gt_shape_and_gt_t_theta == "linear_theta0"] <- "constant"

gt_w_shape_names <- list(
    'linear_theta'=TeX(r'($SPNC^w(t)$: linear)'),
    'quadratic_theta'=TeX(r'($SPNC^w(t)$: quadratic)')
)
gt_w_shape_labeller <- function(variable,value){
    return(gt_w_shape_names[value])
}

result %>% 
    ggplot(aes(x = gt_w_theta, y = power_adj, linetype = factor(ATE_t),
               color = factor(ASPNC_t))) +
    geom_line() +
    facet_grid(~ gt_w_shape, labeller = gt_w_shape_labeller) +
    geom_hline(yintercept = 0.8, linetype = 2) +
    scale_linetype_discrete(name = TeX(r'(Value of $ATE^* = ATE^w$)')) +
    scale_color_discrete(name = TeX(r'(Value of $ASPN^* = ASPN^w$)')) +
    xlab(TeX(r'($\theta_g^w$)')) +
    ylab("power") +
    ggtitle(TeX(r'($SPNC^*(t)$: constant)')) +
    coord_cartesian(ylim = c(0.5, 0.83), xlim = c(-1, 1)) +
    scale_y_continuous(breaks = seq(from = 0.1, to = 1, by = 0.1)) +
    scale_x_continuous(breaks = c(-1,0,1)) +
    # facet_grid(gt_t_theta ~ ASPNC_t, labeller = label_both) +
    theme_bw() + 
    myggfont() +
    theme(legend.text.align = 0, plot.title = element_text(hjust = 0.5))
ggsave(paste0("figure 6 - power.pdf"), width = 9, height = 8)

# result %>% 
#     filter(ASPNC_t %in% c(0.2, 0.4)) %>%
#     mutate(ATE_ASPNC = paste0(ATE_t, ", ", ASPNC_t)) %>%
#     ggplot(aes(x = gt_w_theta, y = power_adj,
#                color = ATE_ASPNC)) +
#     # geom_point() +
#     geom_line() +
#     geom_hline(yintercept = 0.8) +
#     facet_grid( ~ gt_w_shape,
#                labeller = label_both) +
#     theme_bw() + myggfont()
# ggsave(paste0("figure 6 - power.pdf"), width = 10, height = 7.5)
# 
# result %>% 
#     filter(ASPNC_t %in% c(0.2, 0.4)) %>%
#     mutate(ATE_ASPNC = paste0(ATE_t, ", ", ASPNC_t)) %>%
#     ggplot(aes(x = gt_w_theta, y = typeierror_adj,
#                color = ATE_ASPNC)) +
#     # geom_point() +
#     geom_line() +
#     geom_hline(yintercept = 0.05) +
#     facet_grid( ~ gt_w_shape,
#                 labeller = label_both) +
#     theme_bw() + myggfont()
# ggsave(paste0("figure 6 - typeierror.pdf"), width = 10, height = 7.5)
