# Tianchen Qian
# 2022.09.28

rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)
library(rootSolve)
library(MRTSampleSizeBinary)

source("../functions/dgm_noDE.R")
source("../functions/utility.R")
source("../functions/myggplot.R")

# functions copied from the "MRTAnalysisBinary" R package (not published)
# make sure to update R package folder in "Carrie Cheng" if revised these functions
# source("../functions/estimator_EMEE_alwaysCenterA.R")
# source("../functions/find_change_location.R")
# source("../functions/get_alpha_beta_from_multiroot_result.R")

# 1. n vs ATE ----------------------------------------------

# SD: simulation design
SD <- expand.grid(gamma = 0.05,  # type I error
                  b = 0.2, # type II error
                  m = c(30),
                  pt_shape = c("constant 0.6"),
                  ft_t_shape = c("constant"),
                  ft_t_theta = 0,
                  gt_t_shape = c("constant"),
                  gt_t_theta = c(0), 
                  ATE_t = seq(from = 1.05, to = 1.4, by = 0.01),
                  ATE_t_empirical = NA,
                  ASPNC_t = c(0.1, 0.3, 0.5, 0.7),
                  ASPNC_t_empirical = NA,
                  AvgTau_t = c(1),
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

SD %>%
    ggplot(aes(x = ATE_t, y = n, color = factor(ASPNC_t))) +
    geom_line() +
    scale_color_discrete(name = "ASPN") +
    xlab("ATE") +
    ylab("n (sample size)") +
    # ggtitle("constant MEE(t) and SPNC(t)") +
    coord_cartesian(ylim = c(0, 1000)) +
    scale_x_continuous(breaks = c(1.05, 1.1, 1.2, 1.3, 1.4)) +
    theme_bw() + 
    myggfont() +
    theme(legend.text.align = 0, plot.title = element_text(hjust = 0.5))
ggsave("size Drink Less 1 - n vs ATE.pdf", height = 6, width = 6)



# 2. n vs ASPN ----------------------------------------------

# SD: simulation design
SD <- expand.grid(gamma = 0.05,  # type I error
                  b = 0.2, # type II error
                  m = c(30),
                  pt_shape = c("constant 0.6"),
                  ft_t_shape = c("constant"),
                  ft_t_theta = 0,
                  gt_t_shape = c("constant"),
                  gt_t_theta = c(0), 
                  ATE_t = c(1.05, 1.1, 1.2, 1.3, 1.4),
                  ATE_t_empirical = NA,
                  ASPNC_t = seq(from = 0.1, to = 0.7, by = 0.01),
                  ASPNC_t_empirical = NA,
                  AvgTau_t = c(1),
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

SD %>%
    ggplot(aes(x = ASPNC_t, y = n, color = factor(ATE_t))) +
    geom_line() +
    scale_color_discrete(name = "ATE") +
    xlab("ASPN") +
    ylab("n (sample size)") +
    # ggtitle("constant MEE(t) and SPNC(t)") +
    coord_cartesian(ylim = c(0, 1000)) +
    scale_x_continuous(breaks = c(0.1, 0.3, 0.5, 0.7)) +
    theme_bw() + 
    myggfont() +
    theme(legend.text.align = 0, plot.title = element_text(hjust = 0.5))
ggsave("size Drink Less 1 - n vs ASPN.pdf", height = 6, width = 6)


# 3. ATE and ASPN (contours) ----------------------------------------------

# SD: simulation design
SD <- expand.grid(gamma = 0.05,  # type I error
                  b = 0.2, # type II error
                  m = c(30),
                  pt_shape = c("constant 0.6"),
                  ft_t_shape = c("constant"),
                  ft_t_theta = 0,
                  gt_t_shape = c("constant"),
                  gt_t_theta = c(0), 
                  ATE_t = seq(from = 1.1, to = 1.4, by = 0.01),
                  ATE_t_empirical = NA,
                  ASPNC_t = seq(from = 0.2, to = 0.6, by = 0.01),
                  ASPNC_t_empirical = NA,
                  AvgTau_t = c(1),
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

##### fill the SD #####


for (i in 1:nrow(SD)) {
    if (i == 1) {
        print(paste0("for loop total length: ", nrow(SD)))
    }
    if (i %% 100 == 0) {
        cat(i, "")
    }
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

library(metR)
SD %>%
    filter(ATE_t <= 1.3 & ATE_t >= 1.1, ASPNC_t <= 0.6 & ASPNC_t >= 0.2) %>%
    ggplot(aes(x = ATE_t, y = ASPNC_t, z = n)) +
    geom_contour(binwidth = 50) +
    geom_text_contour(skip = 0, stroke = 0.2, size = 5) +
    xlab("ATE") +
    ylab("ASPN") +
    # ggtitle("contour: sample size n") +
    theme_bw() + 
    myggfont() +
    theme(legend.text.align = 0, plot.title = element_text(hjust = 0.5))
ggsave("size Drink Less 1 - n vs ASPN, ATE (contour).pdf", height = 6, width = 6)
