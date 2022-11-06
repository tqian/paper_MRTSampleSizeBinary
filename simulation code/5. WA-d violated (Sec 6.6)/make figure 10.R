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

result10.3 <- readRDS("simulation10.3/result_df_collected_nsim2000.RDS")
result10.3$correct_SPNC <- "3. correct form"
result10.4 <- readRDS("simulation10.4/result_df_collected_nsim2000.RDS")
result10.4$correct_SPNC <- "1. incorrect form and incorrect ASPN"
result10.5 <- readRDS("simulation10.5/result_df_collected_nsim2000.RDS")
result10.5$correct_SPNC <- "2. incorrect form but correct ASPN"
result <- rbind(result10.3, result10.4, result10.5)

result <- result %>% filter(!(gt_t_shape != "constant" & gt_t_theta == 0))


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

result$gt_t_shape_theta <- paste0(result$gt_t_shape, "_", result$gt_t_theta)

result %>% 
    ggplot(aes(x = coef_SC, y = power_adj, linetype = correct_SPNC,
               color = gt_t_shape_theta)) +
    geom_line() +
    # facet_grid(~ gt_t_shape, labeller = gt_t_shape_labeller) +
    geom_hline(yintercept = 0.8, linetype = 2) +
    scale_linetype_discrete(name = TeX(r'(Correctness of $SPNC^w(t)$ and $ASPN^w$)'),
                            labels = c(TeX(r'($SPNC^w(t) \neq SPNC^*(t)$, $ASPN^w \neq ASPN^*$)'),
                                       TeX(r'($SPNC^w(t) \neq SPNC^*(t)$, $ASPN^w = ASPN^*$)'),
                                       TeX(r'($SPNC^w(t) = SPNC^*(t)$, $ASPN^w = ASPN^*$)'))) +
    # scale_color_discrete(guide = "none") +
    scale_color_discrete(name = TeX(r'(Pattern of $g^*(t)^T \alpha^*$)'),
                         labels = c(TeX(r'(constant)'),
                                    TeX(r'(linear, $\theta_g^* = -0.3$)'),
                                    TeX(r'(linear, $\theta_g^* = 0.3$)'),
                                    TeX(r'(quadratic, $\theta_g^* = -0.3$)'),
                                    TeX(r'(quadratic, $\theta_g^* = 0.3$)'))) +
    xlab(TeX(r'($\gamma_2^*$ (serial correlation))')) +
    ylab("power") +
    # ggtitle(TeX(r'($SPNC^w(t)$: constant)')) +
    # coord_cartesian(ylim = c(0.5, 0.83), xlim = c(-1, 1)) +
    # scale_y_continuous(breaks = seq(from = 0.1, to = 1, by = 0.1)) +
    # scale_x_continuous(breaks = c(-1,0,1)) +
    # facet_grid(gt_t_theta ~ ASPNC_t, labeller = label_both) +
    theme_bw() + 
    myggfont(legend_text_size = 14, legend_title_size = 16) +
    theme(legend.text.align = 0, plot.title = element_text(hjust = 0.5))
ggsave(paste0("figure 10 - power.pdf"), width = 9, height = 5)

# result %>% 
#     filter(ASPNC_t %in% c(0.2, 0.4)) %>%
#     mutate(ATE_ASPNC = paste0(ATE_t, ", ", ASPNC_t)) %>%
#     ggplot(aes(x = gt_t_theta, y = power_adj,
#                color = ATE_ASPNC)) +
#     # geom_point() +
#     geom_line() +
#     geom_hline(yintercept = 0.8) +
#     facet_grid( ~ gt_t_shape,
#                labeller = label_both) +
#     theme_bw() + myggfont()
# ggsave(paste0("figure 5 - power.pdf"), width = 10, height = 7.5)
# 
# result %>% 
#     filter(ASPNC_t %in% c(0.2, 0.4)) %>%
#     mutate(ATE_ASPNC = paste0(ATE_t, ", ", ASPNC_t)) %>%
#     ggplot(aes(x = gt_t_theta, y = typeierror_adj,
#                color = ATE_ASPNC)) +
#     # geom_point() +
#     geom_line() +
#     geom_hline(yintercept = 0.05) +
#     facet_grid( ~ gt_t_shape,
#                 labeller = label_both) +
#     theme_bw() + myggfont()
# ggsave(paste0("figure 5 - typeierror.pdf"), width = 10, height = 7.5)
