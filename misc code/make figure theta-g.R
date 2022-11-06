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



# Making plot for constant f ----------------------------------------------

ASPN <- 0.4
m <- 30

df_all <- data.frame()

df <- data.frame(t = 1:m, type = "constant")
df$SPNC <- ASPN
df_all <- rbind(df_all, df)

df_all$logSPNC <- log(df$SPNC)

df_all %>%
    ggplot(aes(x = t, y = logSPNC)) +
    geom_line() +
    theme_bw() +
    xlab("t (decision point)") +
    ylab(TeX(r'($\log\{\textrm{SPNC}(t)\} = g(t)^T\alpha$)')) +
    scale_linetype_discrete(guide = "none") +
    coord_cartesian(ylim = c(-1.7, -0.4), xlim = c(1, 30)) +
    scale_x_continuous(breaks = c(1, 10, 20, 30)) +
    myggfont()
ggsave("theta-g-constant.pdf", width = 6, height = 4)

# Making plot for linear f ------------------------------------------------


compute_alpha_linear <- function(theta, ASPN, m){
    if (theta == 1) {
        compute_ASPN_from_alpha1 <- function(alpha1) {
            alpha0 <- - m * alpha1
            logSPNC <- alpha0 + (1:m) * alpha1
            return(mean(exp(logSPNC)) - ASPN)
        }
        alpha1 <- uniroot(compute_ASPN_from_alpha1, interval = c(-1, 1))$root
        alpha0 <- - m * alpha1
    } else if (theta == 0) {
        alpha1 <- 0
        alpha0 <- log(ASPN)
    } else {
        ratio <- (1+theta) / (1-theta)
        compute_ASPN_from_alpha1 <- function(alpha1) {
            alpha0 <- (ratio * m - 1) / (1 - ratio) * alpha1
            logSPNC <- alpha0 + (1:m) * alpha1
            return(mean(exp(logSPNC)) - ASPN)
        }
        alpha1 <- uniroot(compute_ASPN_from_alpha1, interval = c(-1, 1))$root
        alpha0 <- (ratio * m - 1) / (1 - ratio) * alpha1
    }
    
    return(list(alpha0 = alpha0, alpha1 = alpha1))
}

ASPN <- 0.4
m <- 30

df_all <- data.frame()

df <- data.frame(t = 1:m, type = "1. linear, theta = 0")
theta <- 0
alphas <- compute_alpha_linear(theta, ASPN, m)
alpha0 <- alphas$alpha0
alpha1 <- alphas$alpha1
df$logSPNC <- alpha0 + df$t * alpha1
df_all <- rbind(df_all, df)

df <- data.frame(t = 1:m, type = "2. linear, theta = 0.5")
theta <- 0.5
alphas <- compute_alpha_linear(theta, ASPN, m)
alpha0 <- alphas$alpha0
alpha1 <- alphas$alpha1
df$logSPNC <- alpha0 + df$t * alpha1
df_all <- rbind(df_all, df)

df <- data.frame(t = 1:m, type = "3. linear, theta = -0.5")
theta <- -0.5
alphas <- compute_alpha_linear(theta, ASPN, m)
alpha0 <- alphas$alpha0
alpha1 <- alphas$alpha1
df$logSPNC <- alpha0 + df$t * alpha1
df_all <- rbind(df_all, df)

df <- data.frame(t = 1:m, type = "4. linear, theta = 0.2")
theta <- 0.2
alphas <- compute_alpha_linear(theta, ASPN, m)
alpha0 <- alphas$alpha0
alpha1 <- alphas$alpha1
df$logSPNC <- alpha0 + df$t * alpha1
df_all <- rbind(df_all, df)

df_all %>%
    ggplot(aes(x = t, y = logSPNC, linetype = type)) +
    geom_line() +
    theme_bw() +
    xlab("t (decision point)") +
    ylab(TeX(r'($\log\{\textrm{SPNC}(t)\} = g(t)^T\alpha$)')) +
    annotate("text", x = 0, y = -0.4, label = TeX(r'($\theta_g = -0.5$)'), hjust = 0, size = 5) +
    annotate("text", x = 0, y = -0.85, label = TeX(r'($\theta_g = 0$)'), hjust = 0, size = 5) +
    annotate("text", x = 0, y = -1.02, label = TeX(r'($\theta_g = 0.2$)'), hjust = 0, size = 5) +
    annotate("text", x = 0, y = -1.28, label = TeX(r'($\theta_g = 0.5$)'), hjust = 0, size = 5) +
    scale_linetype_discrete(guide = "none") +
    coord_cartesian(ylim = c(-1.7, -0.4), xlim = c(1, 30)) +
    scale_x_continuous(breaks = c(1, 10, 20, 30)) +
    myggfont()
ggsave("theta-g-linear.pdf", width = 6, height = 4)


# Making plot for quadratic f ------------------------------------------------


compute_alpha_quadratic <- function(theta, ASPN, m){
    if (theta == 1) {
        compute_ASPN_from_alpha2 <- function(alpha2) {
            alpha0 <- m * alpha2
            alpha1 <- - (m+1) * alpha2
            logSPNC <- alpha0 + (1:m) * alpha1 + (1:m)^2 * alpha2
            return(mean(exp(logSPNC)) - ASPN)
        }
        alpha2 <- uniroot(compute_ASPN_from_alpha2, interval = c(-1, 1))$root
        alpha0 <- m * alpha2
        alpha1 <- - (m+1) * alpha2
    } else if (theta == 0) {
        alpha0 <- log(ASPN)
        alpha1 <- 0
        alpha2 <- 0
    } else {
        ratio <- (1+theta) / (1-theta)
        compute_ASPN_from_alpha2 <- function(alpha2) {
            alpha0 <- ((m+1)^2/4 - ratio * (m+1) + ratio) / (1-ratio) * alpha2
            alpha1 <- - (m+1) * alpha2
            logSPNC <- alpha0 + (1:m) * alpha1 + (1:m)^2 * alpha2
            return(mean(exp(logSPNC)) - ASPN)
        }
        alpha2 <- uniroot(compute_ASPN_from_alpha2, interval = c(-1, 1))$root
        alpha0 <- ((m+1)^2/4 - ratio * (m+1) + ratio) / (1-ratio) * alpha2
        alpha1 <- - (m+1) * alpha2
    }
    
    return(list(alpha0 = alpha0, alpha1 = alpha1, alpha2 = alpha2))
}


ASPN <- 0.4
m <- 30

df_all <- data.frame()

# df <- data.frame(t = 1:m, type = "constant")
# df$explogSPNC <- ASPN
# df_all <- rbind(df_all, df)

df <- data.frame(t = 1:m, type = "1. quadratic, theta = 0")
theta <- 0
alphas <- compute_alpha_quadratic(theta, ASPN, m)
alpha0 <- alphas$alpha0
alpha1 <- alphas$alpha1
alpha2 <- alphas$alpha2
df$logSPNC <- alpha0 + (1:m) * alpha1 + (1:m)^2 * alpha2
df_all <- rbind(df_all, df)

df <- data.frame(t = 1:m, type = "2. quadratic, theta = 0.5")
theta <- 0.5
alphas <- compute_alpha_quadratic(theta, ASPN, m)
alpha0 <- alphas$alpha0
alpha1 <- alphas$alpha1
alpha2 <- alphas$alpha2
df$logSPNC <- alpha0 + (1:m) * alpha1 + (1:m)^2 * alpha2
df_all <- rbind(df_all, df)

df <- data.frame(t = 1:m, type = "3. quadratic, theta = -0.5")
theta <- -0.5
alphas <- compute_alpha_quadratic(theta, ASPN, m)
alpha0 <- alphas$alpha0
alpha1 <- alphas$alpha1
alpha2 <- alphas$alpha2
df$logSPNC <- alpha0 + (1:m) * alpha1 + (1:m)^2 * alpha2
df_all <- rbind(df_all, df)

df <- data.frame(t = 1:m, type = "4. quadratic, theta = 0.2")
theta <- 0.2
alphas <- compute_alpha_quadratic(theta, ASPN, m)
alpha0 <- alphas$alpha0
alpha1 <- alphas$alpha1
alpha2 <- alphas$alpha2
df$logSPNC <- alpha0 + (1:m) * alpha1 + (1:m)^2 * alpha2
df_all <- rbind(df_all, df)

df_all %>%
    ggplot(aes(x = t, y = logSPNC, linetype = type)) +
    geom_line() +
    theme_bw() +
    xlab("t (decision point)") +
    ylab(TeX(r'($\log\{\textrm{SPNC}(t)\} = g(t)^T\alpha$)')) +
    annotate("text", x = 13, y = -0.5, label = TeX(r'($\theta_g = -0.5$)'), hjust = 0, size = 5) +
    annotate("text", x = 13, y = -0.83, label = TeX(r'($\theta_g = 0$)'), hjust = 0, size = 5) +
    annotate("text", x = 13, y = -1.13, label = TeX(r'($\theta_g = 0.2$)'), hjust = 0, size = 5) +
    annotate("text", x = 13, y = -1.32, label = TeX(r'($\theta_g = 0.5$)'), hjust = 0, size = 5) +
    scale_linetype_discrete(guide = "none") +
    coord_cartesian(ylim = c(-1.7, -0.4), xlim = c(1, 30)) +
    scale_x_continuous(breaks = c(1, 10, 20, 30)) +
    myggfont()
ggsave("theta-g-quadratic.pdf", width = 6, height = 4)
