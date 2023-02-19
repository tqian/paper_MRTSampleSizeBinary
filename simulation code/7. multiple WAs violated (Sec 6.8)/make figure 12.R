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

result12 <- readRDS("simulation12/result_df_collected_nsim2000.RDS")

result12 %>%
    ggplot(aes(x = power_adj)) +
    geom_histogram(aes(y = after_stat(count)), binwidth = 0.002, color = "white") +
    # geom_vline(xintercept = 0.8, linetype = 2, color = "blue", linewidth = 1.2) +
    # stat_function(fun = dnorm,
    #               args = list(mean = mean(result$power_adj),
    #                           sd = sd(result$power_adj)),
    #               color = "orange",
    #               size = 1.5) +
    xlab("power") + 
    theme_bw() + myggfont()
ggsave(paste0("figure 12 - power.pdf"), width = 7, height = 5)
