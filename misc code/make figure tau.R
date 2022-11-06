### plot making ###
# Binary MRT sample size paper
#
# see _plot making log.md for simulation set up
#
# Tianchen Qian
# 2022.09.20


rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(tidyverse)
library(gridExtra) 
library(rootSolve)
library(latex2exp)
source("../../functions/myggplot.R")
# source("../../functions/dgm_noDE.R")

m <- 30

df_all <- data.frame()

df <- data.frame(t = 1:m, type = "constant")
df$tau <- 0.7
df_all <- rbind(df_all, df)

df <- data.frame(t = 1:m, type = "linear, theta = -0.2")
df$tau <- 0.7 + seq(from = -0.2, to = 0.2, length.out = m)
df_all <- rbind(df_all, df)

# df <- data.frame(t = 1:m, type = "linear, theta = 0.2")
# df$tau <- 0.7 + seq(from = 0.2, to = -0.2, length.out = m)
# df_all <- rbind(df_all, df)

df <- data.frame(t = 1:m, type = "linear, theta = 0.3")
df$tau <- 0.7 + seq(from = 0.3, to = -0.3, length.out = m)
df_all <- rbind(df_all, df)

df <- data.frame(t = 1:m, type = "linear, theta = -0.3")
df$tau <- 0.7 + seq(from = -0.3, to = 0.3, length.out = m)
df_all <- rbind(df_all, df)

df <- data.frame(t = (100:(100*m))/100, type = "sine, theta = 0.2")
df$tau <- 0.7 + 0.2 * sin(df$t)
df_all <- rbind(df_all, df)

df <- data.frame(t = (100:(100*m))/100, type = "sine, theta = 0.1")
df$tau <- 0.7 + 0.1 * sin(df$t)
df_all <- rbind(df_all, df)



df_all %>%
    ggplot(aes(x = t, y = tau, linetype = type)) +
    geom_line() +
    theme_bw() +
    xlab("t (decision point)") +
    ylab(TeX(r'($\tau(t)$)')) +
    annotate("text", x = 30.5, y = 0.4, label = TeX(r'(linear, $\theta_\tau = 0.3$)'), hjust = 0, size = 5) +
    annotate("text", x = 30.5, y = 0.5, label = TeX(r'(periodic, $\theta_\tau = 0.2$)'), hjust = 0, size = 5) +
    annotate("text", x = 30.5, y = 0.6, label = TeX(r'(periodic, $\theta_\tau = 0.1$)'), hjust = 0, size = 5) +
    annotate("text", x = 30.5, y = 0.7, label = TeX(r'(constant)'), hjust = 0, size = 5) +
    annotate("text", x = 30.5, y = 0.9, label = TeX(r'(linear, $\theta_\tau = -0.2$)'), hjust = 0, size = 5) +
    annotate("text", x = 30.5, y = 1, label = TeX(r'(linear, $\theta_\tau = -0.3$)'), hjust = 0, size = 5) +
    scale_linetype_discrete(guide = "none") +
    coord_cartesian(xlim = c(1, 34)) +
    scale_x_continuous(breaks = c(1, 10, 20, 30)) +
    myggfont()
ggsave("theta-tau.pdf", width = 12, height = 4)
