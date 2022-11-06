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

ATE <- 1.2
m <- 30

df_all <- data.frame()

df <- data.frame(t = 1:m, type = "constant")
df$MEE <- log(ATE)
df_all <- rbind(df_all, df)

df_all %>%
    ggplot(aes(x = t, y = MEE)) +
    geom_line() +
    theme_bw() +
    xlab("t (decision point)") +
    ylab(TeX(r'($\textrm{MEE}(t) = f(t)^T\beta$)')) +
    scale_linetype_discrete(guide = "none") +
    coord_cartesian(ylim = c(0, 0.5), xlim = c(1, 30)) +
    scale_x_continuous(breaks = c(1, 10, 20, 30)) +
    myggfont()
ggsave("theta-f-constant.pdf", width = 6, height = 4)

# Making plot for linear f ------------------------------------------------


compute_beta_linear <- function(theta, ATE, m){
    if (theta == 1) {
        compute_ATE_from_beta1 <- function(beta1) {
            beta0 <- - m * beta1
            MEE <- beta0 + (1:m) * beta1
            return(mean(exp(MEE)) - ATE)
        }
        beta1 <- uniroot(compute_ATE_from_beta1, interval = c(-1, 1))$root
        beta0 <- - m * beta1
    } else if (theta == 0) {
        beta1 <- 0
        beta0 <- log(ATE)
    } else {
        ratio <- (1+theta) / (1-theta)
        compute_ATE_from_beta1 <- function(beta1) {
            beta0 <- (ratio * m - 1) / (1 - ratio) * beta1
            MEE <- beta0 + (1:m) * beta1
            return(mean(exp(MEE)) - ATE)
        }
        beta1 <- uniroot(compute_ATE_from_beta1, interval = c(-1, 1))$root
        beta0 <- (ratio * m - 1) / (1 - ratio) * beta1
    }
    
    return(list(beta0 = beta0, beta1 = beta1))
}

ATE <- 1.2
m <- 30

df_all <- data.frame()

df <- data.frame(t = 1:m, type = "1. linear, theta = 0")
theta <- 0
betas <- compute_beta_linear(theta, ATE, m)
beta0 <- betas$beta0
beta1 <- betas$beta1
df$MEE <- beta0 + df$t * beta1
df_all <- rbind(df_all, df)

df <- data.frame(t = 1:m, type = "2. linear, theta = 1")
theta <- 1
betas <- compute_beta_linear(theta, ATE, m)
beta0 <- betas$beta0
beta1 <- betas$beta1
df$MEE <- beta0 + df$t * beta1
df_all <- rbind(df_all, df)

df <- data.frame(t = 1:m, type = "3. linear, theta = -1")
theta <- -1
betas <- compute_beta_linear(theta, ATE, m)
beta0 <- betas$beta0
beta1 <- betas$beta1
df$MEE <- beta0 + df$t * beta1
df_all <- rbind(df_all, df)

df <- data.frame(t = 1:m, type = "4. linear, theta = 0.4")
theta <- 0.4
betas <- compute_beta_linear(theta, ATE, m)
beta0 <- betas$beta0
beta1 <- betas$beta1
df$MEE <- beta0 + df$t * beta1
df_all <- rbind(df_all, df)

df_all %>%
    ggplot(aes(x = t, y = MEE, linetype = type)) +
    geom_line() +
    theme_bw() +
    xlab("t (decision point)") +
    ylab(TeX(r'($\textrm{MEE}(t) = f(t)^T\beta$)')) +
    annotate("text", x = 0, y = 0.06, label = TeX(r'($\theta_f = -1$)'), hjust = 0, size = 5) +
    annotate("text", x = 0, y = 0.16, label = TeX(r'($\theta_f = 0$)'), hjust = 0, size = 5) +
    annotate("text", x = 0, y = 0.29, label = TeX(r'($\theta_f = 0.4$)'), hjust = 0, size = 5) +
    annotate("text", x = 0, y = 0.39, label = TeX(r'($\theta_f = 1$)'), hjust = 0, size = 5) +
    scale_linetype_discrete(guide = "none") +
    coord_cartesian(ylim = c(0, 0.5), xlim = c(1, 30)) +
    scale_x_continuous(breaks = c(1, 10, 20, 30)) +
    myggfont()
ggsave("theta-f-linear.pdf", width = 6, height = 4)


# Making plot for quadratic f ------------------------------------------------


compute_beta_quadratic <- function(theta, ATE, m){
    if (theta == 1) {
        compute_ATE_from_beta2 <- function(beta2) {
            beta0 <- m * beta2
            beta1 <- - (m+1) * beta2
            MEE <- beta0 + (1:m) * beta1 + (1:m)^2 * beta2
            return(mean(exp(MEE)) - ATE)
        }
        beta2 <- uniroot(compute_ATE_from_beta2, interval = c(-1, 1))$root
        beta0 <- m * beta2
        beta1 <- - (m+1) * beta2
    } else if (theta == 0) {
        beta0 <- log(ATE)
        beta1 <- 0
        beta2 <- 0
    } else {
        ratio <- (1+theta) / (1-theta)
        compute_ATE_from_beta2 <- function(beta2) {
            beta0 <- ((m+1)^2/4 - ratio * (m+1) + ratio) / (1-ratio) * beta2
            beta1 <- - (m+1) * beta2
            MEE <- beta0 + (1:m) * beta1 + (1:m)^2 * beta2
            return(mean(exp(MEE)) - ATE)
        }
        beta2 <- uniroot(compute_ATE_from_beta2, interval = c(-1, 1))$root
        beta0 <- ((m+1)^2/4 - ratio * (m+1) + ratio) / (1-ratio) * beta2
        beta1 <- - (m+1) * beta2
    }
    
    return(list(beta0 = beta0, beta1 = beta1, beta2 = beta2))
}


ATE <- 1.2
m <- 30

df_all <- data.frame()

# df <- data.frame(t = 1:m, type = "constant")
# df$expMEE <- ATE
# df_all <- rbind(df_all, df)

df <- data.frame(t = 1:m, type = "1. quadratic, theta = 0")
theta <- 0
betas <- compute_beta_quadratic(theta, ATE, m)
beta0 <- betas$beta0
beta1 <- betas$beta1
beta2 <- betas$beta2
df$MEE <- beta0 + (1:m) * beta1 + (1:m)^2 * beta2
df_all <- rbind(df_all, df)

df <- data.frame(t = 1:m, type = "2. quadratic, theta = 1")
theta <- 1
betas <- compute_beta_quadratic(theta, ATE, m)
beta0 <- betas$beta0
beta1 <- betas$beta1
beta2 <- betas$beta2
df$MEE <- beta0 + (1:m) * beta1 + (1:m)^2 * beta2
df_all <- rbind(df_all, df)

df <- data.frame(t = 1:m, type = "3. quadratic, theta = -1")
theta <- -1
betas <- compute_beta_quadratic(theta, ATE, m)
beta0 <- betas$beta0
beta1 <- betas$beta1
beta2 <- betas$beta2
df$MEE <- beta0 + (1:m) * beta1 + (1:m)^2 * beta2
df_all <- rbind(df_all, df)

df <- data.frame(t = 1:m, type = "4. quadratic, theta = 0.4")
theta <- 0.4
betas <- compute_beta_quadratic(theta, ATE, m)
beta0 <- betas$beta0
beta1 <- betas$beta1
beta2 <- betas$beta2
df$MEE <- beta0 + (1:m) * beta1 + (1:m)^2 * beta2
df_all <- rbind(df_all, df)

df_all %>%
    ggplot(aes(x = t, y = MEE, linetype = type)) +
    geom_line() +
    theme_bw() +
    xlab("t (decision point)") +
    ylab(TeX(r'($\textrm{MEE}(t) = f(t)^T\beta$)')) +
    annotate("text", x = 13, y = 0.04, label = TeX(r'($\theta_f = -1$)'), hjust = 0, size = 5) +
    annotate("text", x = 13, y = 0.16, label = TeX(r'($\theta_f = 0$)'), hjust = 0, size = 5) +
    annotate("text", x = 13, y = 0.25, label = TeX(r'($\theta_f = 0.4$)'), hjust = 0, size = 5) +
    annotate("text", x = 13, y = 0.31, label = TeX(r'($\theta_f = 1$)'), hjust = 0, size = 5) +
    scale_linetype_discrete(guide = "none") +
    coord_cartesian(ylim = c(0, 0.5), xlim = c(1, 30)) +
    scale_x_continuous(breaks = c(1, 10, 20, 30)) +
    myggfont()
# \theta_{\textrm{quad}}
ggsave("theta-f-quadratic.pdf", width = 6, height = 4)
