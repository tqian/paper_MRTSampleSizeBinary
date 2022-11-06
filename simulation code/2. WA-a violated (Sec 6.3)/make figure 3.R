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

result2.4 <- readRDS("simulation2.4/result_df_collected_nsim2000.RDS")
result2.5 <- readRDS("simulation2.5/result_df_collected_nsim2000.RDS")
result <- rbind(result2.4, result2.5)

# result$ft_t_shape_and_theta <- paste0(result$ft_t_shape, result$ft_t_theta)
result$ASPNC_t_and_gt_shape_and_gt_t_theta <- paste0(result$ASPNC_t, result$gt_t_shape, result$gt_t_theta)
result$gt_shape_and_gt_t_theta <- paste0(result$gt_t_shape, result$gt_t_theta)
result$ATE_t_and_ASPNC_t <- paste0(result$ATE_t, result$ASPNC_t)

result$gt_shape_and_gt_t_theta[result$gt_shape_and_gt_t_theta == "linear_theta0"] <- "constant"

ft_w_shape_names <- list(
    'linear_theta'=TeX(r'($MEE^w(t)$: linear)'),
    'quadratic_theta'=TeX(r'($MEE^w(t)$: quadratic)')
)
ft_w_shape_labeller <- function(variable,value){
    return(ft_w_shape_names[value])
}

p3 <- result %>% 
    filter(gt_shape_and_gt_t_theta != "quadratic_theta0") %>%
    ggplot(aes(x = ft_w_theta, y = power_adj, linetype = factor(gt_shape_and_gt_t_theta),
               color = factor(ATE_t_and_ASPNC_t))) +
    geom_line() +
    facet_grid(~ ft_w_shape, labeller = ft_w_shape_labeller) +
    geom_hline(yintercept = 0.8, linetype = 2) +
    scale_linetype_discrete(name = TeX(r'(Pattern of $SPNC^*(t) = SPNC^w(t)$)'),
                            labels = c(TeX(r'(constant)'),
                                       TeX(r'(linear, $\theta_g^* = \theta_g^w = -0.3$)'),
                                       TeX(r'(linear, $\theta_g^* = \theta_g^w = 0.3$)'),
                                       TeX(r'(quadratic, $\theta_g^* = \theta_g^w = -0.3$)'),
                                       TeX(r'(quadratic, $\theta_g^* = \theta_g^w = 0.3$)'))) +
    # scale_color_discrete(name = TeX(r'(ATE and ASPN)'),
    #                      labels = c(TeX(r'($ATE^* = ATE^w = 1.2$, $ASPN^* = ASPN^w = 0.2$)'),
    #                                 TeX(r'($ATE^* = ATE^w = 1.2$, $ASPN^* = ASPN^w = 0.3$)'),
    #                                 TeX(r'($ATE^* = ATE^w = 1.2$, $ASPN^* = ASPN^w = 0.4$)'),
    #                                 TeX(r'($ATE^* = ATE^w = 1.3$, $ASPN^* = ASPN^w = 0.2$)'),
    #                                 TeX(r'($ATE^* = ATE^w = 1.3$, $ASPN^* = ASPN^w = 0.3$)'),
    #                                 TeX(r'($ATE^* = ATE^w = 1.3$, $ASPN^* = ASPN^w = 0.4$)'))) +
    scale_color_discrete(name = TeX(r'($ATE^* = ATE^w$, $ASPN^* = ASPN^w$)'),
                         labels = c(TeX(r'($ATE = 1.2$, $ASPN = 0.2$)'),
                                    TeX(r'($ATE = 1.2$, $ASPN = 0.3$)'),
                                    TeX(r'($ATE = 1.2$, $ASPN = 0.4$)'),
                                    TeX(r'($ATE = 1.3$, $ASPN = 0.2$)'),
                                    TeX(r'($ATE = 1.3$, $ASPN = 0.3$)'),
                                    TeX(r'($ATE = 1.3$, $ASPN = 0.4$)'))) +
    xlab(TeX(r'($\theta_f^w$)')) +
    ylab("power") +
    ggtitle(TeX(r'($MEE^*(t)$: constant)')) +
    coord_cartesian(ylim = c(0.5, 0.83)) +
    scale_y_continuous(breaks = seq(from = 0.1, to = 1, by = 0.1)) +
    scale_x_continuous(breaks = c(-1,0,1)) +
    # facet_grid(gt_t_theta ~ ASPNC_t, labeller = label_both) +
    theme_bw() + 
    myggfont() +
    theme(legend.text.align = 0, plot.title = element_text(hjust = 0.5))
p3
ggsave(paste0("figure 3 - power.pdf"), width = 10, height = 8)

# result %>% 
#     filter(ASPNC_t %in% c(0.2, 0.4)) %>%
#     mutate(ATE_ASPNC = paste0(ATE_t, ", ", ASPNC_t)) %>%
#     ggplot(aes(x = ft_w_theta, y = power_adj, linetype = gt_t_shape,
#                color = ATE_ASPNC)) +
#     # geom_point() +
#     geom_line() +
#     geom_hline(yintercept = 0.8) +
#     facet_grid(ft_w_shape ~ gt_t_theta,
#                labeller = label_both) +
#     theme_bw() + myggfont()
# ggsave(paste0("figure 3 - power.pdf"), width = 10, height = 7.5)
# 
# result %>% 
#     filter(ASPNC_t %in% c(0.2, 0.4)) %>%
#     mutate(ATE_ASPNC = paste0(ATE_t, ", ", ASPNC_t)) %>%
#     ggplot(aes(x = ft_w_theta, y = typeierror_adj, linetype = gt_t_shape,
#                color = ATE_ASPNC)) +
#     # geom_point() +
#     geom_line() +
#     geom_hline(yintercept = 0.05) +
#     facet_grid(ft_w_shape ~ gt_t_theta,
#                labeller = label_both) +
#     theme_bw() + myggfont()
# ggsave(paste0("figure 3 - typeierror.pdf"), width = 10, height = 7.5)
