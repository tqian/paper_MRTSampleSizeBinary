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
library(gridExtra) 
source("../../functions/myggplot.R")

nsim <- 2000
result <- readRDS("simulation1.2/result_df_collected_nsim2000.RDS")

power_1 <-
    result %>%
    ggplot(aes(x = power_adj)) +
    geom_histogram(aes(y = after_stat(count)), binwidth = 0.002, color = "white") +
    geom_vline(xintercept = 0.8, linetype = 2, color = "blue", size = 1.2) +
    # stat_function(fun = dnorm,
    #               args = list(mean = mean(result$power_adj),
    #                           sd = sd(result$power_adj)),
    #               color = "orange",
    #               size = 1.5) +
    xlab("power") + 
    theme_bw() + myggfont()

power_2 <- 
    result %>%
    ggplot(aes(x = n, y = power_adj)) +
    geom_point(alpha = 0.6, size = 0.8) +
    geom_hline(yintercept = 0.8, linetype = 2, color = "blue", size = 1.2) +
    # geom_smooth(method = "loess", size = 1.5, se = FALSE) +
    xlab("output sample size") +
    ylab("power") +
    scale_x_continuous(breaks = c(20, 60, 100, 200, 300)) +
    theme_bw() + myggfont()

pdf("figure 0 - power.pdf", width = 10, height = 5)
grid.arrange(power_1, power_2, ncol = 2)
dev.off()


typeierror_1 <-
    result %>%
    ggplot(aes(x = typeierror_adj)) +
    geom_histogram(aes(y = ..density..), binwidth = 0.002, color = "white") +
    # stat_function(fun = dnorm,
    #               args = list(mean = mean(result$typeierror_adj),
    #                           sd = sd(result$typeierror_adj)),
    #               color = "orange",
    #               size = 1.5) +
    xlab("type I error") + 
    theme_bw() + myggfont()

typeierror_2 <- 
    result %>%
    ggplot(aes(x = n, y = typeierror_adj)) +
    geom_point(alpha = 0.6, size = 0.8) +
    # geom_smooth(method = "loess", size = 1.5) +
    xlab("output sample size") +
    ylab("type I error") +
    scale_x_continuous(breaks = c(20, 60, 100, 200, 300)) +
    theme_bw() + myggfont()

pdf("figure 0 - typeirror.pdf", width = 10, height = 5)
grid.arrange(typeierror_1, typeierror_2, ncol = 2)
dev.off()
