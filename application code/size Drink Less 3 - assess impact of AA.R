# Tianchen Qian
# 2022.09.28

rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)
library(rootSolve)
library(latex2exp)
library(MRTSampleSizeBinary)

source("../functions/dgm_noDE.R")
source("../functions/utility.R")
source("../functions/myggplot.R")

# functions copied from the "MRTAnalysisBinary" R package (not published)
# make sure to update R package folder in "Carrie Cheng" if revised these functions
# source("../functions/estimator_EMEE_alwaysCenterA.R")
# source("../functions/find_change_location.R")
# source("../functions/get_alpha_beta_from_multiroot_result.R")

# 1. n vs AvgTau_t ----------------------------------------------

# SD: simulation design
SD <- expand.grid(gamma = 0.05,  # type I error
                  b = 0.2, # type II error
                  m = c(30),
                  pt_shape = c("constant 0.6"),
                  ft_t_shape = c("constant"),
                  ft_t_theta = 0,
                  gt_t_shape = c("constant"),
                  gt_t_theta = c(0), 
                  ATE_t = 1.15,
                  ATE_t_empirical = NA,
                  ASPNC_t = 0.3,
                  ASPNC_t_empirical = NA,
                  AvgTau_t = seq(from = 0.5, to = 1, by = 0.01),
                  tau_t_shape = c("constant"),
                  tau_t_theta = 0,
                  beta_t0 = NA,
                  beta_t1 = NA,
                  beta_t2 = NA,
                  alpha_t0 = NA,
                  alpha_t1 = NA,
                  alpha_t2 = NA,
                  p = NA,
                  q = NA,
                  
                  n = NA,
                  
                  stringsAsFactors = FALSE)

SD <- SD %>% filter(! (ft_t_shape == "constant" & ft_t_theta != 0))

##### fill the SD #####

for (i in 1:nrow(SD)) {
    
    ## gamma, b, m
    gamma <- SD$gamma[i]
    b <- SD$b[i]
    m <- SD$m[i]
    
    ## pt
    pt <- construct_pt(m, SD$pt_shape[i])
    
    ## taut_t
    taut_t <- construct_taut_theta(m, SD$AvgTau_t[i], SD$tau_t_shape[i], SD$tau_t_theta[i])
    
    ## gt_t, alpha_t
    gt_t <- construct_ftgt(m, SD$gt_t_shape[i])
    SD$q[i] <- ncol(gt_t)
    alpha_t <- solve_alpha(SD$gt_t_shape[i], m, SD$ASPNC_t[i], gt_t, taut_t, 
                           theta = SD$gt_t_theta[i])
    SD[i, c("alpha_t0", "alpha_t1", "alpha_t2")[1:SD$q[i]]] <- alpha_t
    SD$ASPNC_t_empirical[i] <- compute_ASPNC(gt_t, alpha_t, taut_t)
    stopifnot(all.equal(SD$ASPNC_t[i], SD$ASPNC_t_empirical[i]))
    
    ## ft_t, beta_t
    ft_t <- construct_ftgt(m, SD$ft_t_shape[i])
    SD$p[i] <- ncol(ft_t)
    beta_t <- solve_beta(SD$ft_t_shape[i], m, SD$ATE_t[i], ft_t, gt_t, alpha_t, 
                         taut_t, theta = SD$ft_t_theta[i])
    SD[i, c("beta_t0", "beta_t1", "beta_t2")[1:SD$p[i]]] <- beta_t
    SD$ATE_t_empirical[i] <- compute_ATE(gt_t, alpha_t, ft_t, beta_t, taut_t)
    stopifnot(all.equal(SD$ATE_t[i], SD$ATE_t_empirical[i]))
    
    ## sample size n
    tryCatch({
        SD$n[i] <- mrt_binary_ss(taut_t, ft_t, gt_t, beta_t, alpha_t, pt, gamma, b,
                                 less_than_10_possible = TRUE)
        ## generate a data set to make sure there are no probabilities outside [0,1]
        dta <- dgm_noDE(SD$n[i], m, ft_t, beta_t, gt_t, alpha_t, taut_t, pt)
    },
    error = function(cond){
        message("iteration ", i, ": ", cond)
        return(NA)
    })
}

summary(SD)

n_constant_MEE <- SD$n[SD$ft_t_shape == "constant"]

SD %>% 
    ggplot(aes(x = AvgTau_t, y = n)) +
    geom_line() +
    xlab("AA") +
    ylab("n (sample size)") +
    # ggtitle(TeX(r'(constant $\tau(t)$)')) +
    # coord_cartesian(ylim = c(60, 160)) +
    # scale_x_continuous(breaks = c(1.05, 1.1, 1.2, 1.3, 1.4)) +
    theme_bw() + 
    myggfont() +
    theme(legend.text.align = 0, plot.title = element_text(hjust = 0.5))
ggsave("size Drink Less 3 - n vs AA.pdf", height = 4, width = 5)



# 2. n vs theta_g ----------------------------------------------

# SD: simulation design
SD <- expand.grid(gamma = 0.05,  # type I error
                  b = 0.2, # type II error
                  m = c(30),
                  pt_shape = c("constant 0.6"),
                  ft_t_shape = c("constant", "linear_theta", "quadratic_theta"),
                  ft_t_theta = c(-0.5, 0, 0.5),
                  gt_t_shape = c("quadratic_theta"),
                  gt_t_theta = c(0), 
                  ATE_t = 1.15,
                  ATE_t_empirical = NA,
                  ASPNC_t = 0.3,
                  ASPNC_t_empirical = NA,
                  AvgTau_t = 0.7,
                  tau_t_shape = c("constant", "linear", "sine"),
                  tau_t_theta = seq(from = -0.3, by = 0.01, to = 0.3),
                  beta_t0 = NA,
                  beta_t1 = NA,
                  beta_t2 = NA,
                  alpha_t0 = NA,
                  alpha_t1 = NA,
                  alpha_t2 = NA,
                  p = NA,
                  q = NA,
                  
                  n = NA,
                  
                  stringsAsFactors = FALSE)

SD <- SD %>% filter(! (ft_t_shape == "constant" & ft_t_theta != 0))
SD <- SD %>% filter(! (tau_t_shape == "constant" & tau_t_theta != 0))

##### fill the SD #####

for (i in 1:nrow(SD)) {
    
    ## gamma, b, m
    gamma <- SD$gamma[i]
    b <- SD$b[i]
    m <- SD$m[i]
    
    ## pt
    pt <- construct_pt(m, SD$pt_shape[i])
    
    ## taut_t
    taut_t <- construct_taut_theta(m, SD$AvgTau_t[i], SD$tau_t_shape[i], SD$tau_t_theta[i])
    
    ## gt_t, alpha_t
    gt_t <- construct_ftgt(m, SD$gt_t_shape[i])
    SD$q[i] <- ncol(gt_t)
    alpha_t <- solve_alpha(SD$gt_t_shape[i], m, SD$ASPNC_t[i], gt_t, taut_t, 
                           theta = SD$gt_t_theta[i])
    SD[i, c("alpha_t0", "alpha_t1", "alpha_t2")[1:SD$q[i]]] <- alpha_t
    SD$ASPNC_t_empirical[i] <- compute_ASPNC(gt_t, alpha_t, taut_t)
    stopifnot(all.equal(SD$ASPNC_t[i], SD$ASPNC_t_empirical[i]))
    
    ## ft_t, beta_t
    ft_t <- construct_ftgt(m, SD$ft_t_shape[i])
    SD$p[i] <- ncol(ft_t)
    beta_t <- solve_beta(SD$ft_t_shape[i], m, SD$ATE_t[i], ft_t, gt_t, alpha_t, 
                         taut_t, theta = SD$ft_t_theta[i])
    SD[i, c("beta_t0", "beta_t1", "beta_t2")[1:SD$p[i]]] <- beta_t
    SD$ATE_t_empirical[i] <- compute_ATE(gt_t, alpha_t, ft_t, beta_t, taut_t)
    stopifnot(all.equal(SD$ATE_t[i], SD$ATE_t_empirical[i]))
    
    ## sample size n
    tryCatch({
        SD$n[i] <- mrt_binary_ss(taut_t, ft_t, gt_t, beta_t, alpha_t, pt, gamma, b,
                                 less_than_10_possible = TRUE, exact = TRUE)
        ## generate a data set to make sure there are no probabilities outside [0,1]
        dta <- dgm_noDE(SD$n[i], m, ft_t, beta_t, gt_t, alpha_t, taut_t, pt)
    },
    error = function(cond){
        message("iteration ", i, ": ", cond)
        return(NA)
    })
}

summary(SD)

n_constant_tau <- SD$n[SD$tau_t_shape == "constant"]

SD$ft_t_shape_and_ft_t_theta <- paste0(SD$ft_t_shape, SD$ft_t_theta)

SD %>% 
    filter(tau_t_shape != "constant") %>%
    filter(! (ft_t_shape == "linear_theta" & ft_t_theta == 0)) %>%
    filter(! (ft_t_shape == "quadratic_theta" & ft_t_theta == 0)) %>%
    ggplot(aes(x = tau_t_theta, y = n, linetype = tau_t_shape, color = ft_t_shape_and_ft_t_theta)) +
    geom_line() +
    scale_color_discrete(name = "Pattern of MEE(t)",
                             labels = c(TeX(r'(constant)'),
                                        TeX(r'(linear, $\theta_f = -0.5$)'),
                                        TeX(r'(linear, $\theta_f = 0.5$)'),
                                        TeX(r'(quadratic, $\theta_f = -0.5$)'),
                                        TeX(r'(quadratic, $\theta_f = 0.5$)'))) +
    scale_linetype_discrete(name = TeX(r'(Pattern of $\tau(t)$)'),
                         labels = c("linear", "periodic")) +
    xlab(TeX(r'($\theta_\tau$)')) +
    ylab("n (sample size)") +
    theme_bw() + 
    myggfont() +
    theme(legend.text.align = 0, plot.title = element_text(hjust = 0.5))
ggsave("size Drink Less 3 - n vs theta_tau.pdf", height = 4, width = 7)
