### simulation ###
# Binary MRT sample size paper
#
# see _simulation observations.md for simulation set up
#
# Tianchen Qian
# 2022.05.07

rm(list = ls())

# library(tidyverse)
library(rootSolve)
library(MRTSampleSizeBinary)

source("../../../functions/dgm_noDE.R")
source("../../../functions/utility.R")

# functions copied from the "MRTAnalysisBinary" R package (not published)
# make sure to update R package folder in "Carrie Cheng" if revised these functions
source("../../../functions/estimator_EMEE_alwaysCenterA.R")
source("../../../functions/find_change_location.R")
source("../../../functions/get_alpha_beta_from_multiroot_result.R")

# 1. Simulation design setup ----------------------------------------------

# SD: simulation design
SD <- expand.grid(gamma = 0.05,  # type I error
                  b = 0.2, # type II error
                  # m = c(30, 60),
                  m = c(30),
                  pt_shape = c("constant 0.5", "linear 0.5 pm 0.1"),
                  # pt_shape = c("constant 0.5"),
                  ft_t_shape = c("quadratic_theta", "linear_theta", "constant"),
                  ft_t_theta = c(-0.5, 0, 0.5),
                  gt_t_shape = c("quadratic_theta", "linear_theta", "constant"),
                  gt_t_theta = c(-0.3, 0, 0.3), 
                  # for specifying gt when gt_t_shape is linear or quadratic
                  # ATE_t = c(1.2, 1.3, 1.4, 1.5),
                  ATE_t = c(1.2, 1.3),
                  ATE_t_empirical = NA,
                  ASPNC_t = c(0.2, 0.4),
                  ASPNC_t_empirical = NA,
                  AvgTau_t = c(0.5, 0.8, 1),
                  # AvgTau_t = c(1),
                  tau_t_shape = c("constant", "linear pm 0.2", "sin pm 0.2"),
                  beta_t0 = NA,
                  beta_t1 = NA,
                  beta_t2 = NA,
                  alpha_t0 = NA,
                  alpha_t1 = NA,
                  alpha_t2 = NA,
                  p = NA,
                  q = NA,
                  
                  ft_w_shape = NA,
                  ft_w_theta = NA,
                  gt_w_shape = NA,
                  gt_w_theta = NA,
                  beta_w0 = NA,
                  beta_w1 = NA,
                  beta_w2 = NA,
                  alpha_w0 = NA,
                  alpha_w1 = NA,
                  alpha_w2 = NA,
                  ATE_w = NA,
                  ATE_w_empirical = NA,
                  ASPNC_w = NA,
                  ASPNC_w_empirical = NA,
                  AvgTau_w = NA,
                  tau_w_shape = NA,
                  p_w = NA,
                  q_w = NA,
                  n = NA,
                  
                  stringsAsFactors = FALSE)

# remove simulation settings where pt*ft is not in the linear span of gt
SD <- SD[!(SD$pt_shape == "linear 0.5 pm 0.1" &
               SD$ft_t_shape == "quadratic_theta"), ]
SD <- SD[!(SD$pt_shape == "linear 0.5 pm 0.1" &
               SD$ft_t_shape == "linear_theta" &
               SD$gt_t_shape != "quadratic_theta"), ]
SD <- SD[!(SD$pt_shape == "linear 0.5 pm 0.1" &
               SD$ft_t_shape == "constant" &
               SD$gt_t_shape == "constant"), ]

SD <- SD[!(SD$ft_t_shape == "linear_theta" &
               SD$gt_t_shape == "constant"), ]
SD <- SD[!(SD$ft_t_shape == "quadratic_theta" &
               SD$gt_t_shape != "quadratic_theta"), ]

SD <- SD[!(SD$ft_t_shape == "constant" &
               SD$ft_t_theta != 0), ]
SD <- SD[!(SD$gt_t_shape == "constant" &
               SD$gt_t_theta != 0), ]

SD <- SD[!(SD$AvgTau_t == 1 &
               SD$tau_t_shape != "constant"), ]


## add stuff that are known by simulation design
# taut
SD$AvgTau_w <- SD$AvgTau_t
SD$tau_w_shape <- SD$tau_t_shape
# gt
SD$gt_w_shape <- SD$gt_t_shape
SD$gt_w_theta <- SD$gt_t_theta
SD$ASPNC_w <- SD$ASPNC_t
# ft
SD$ft_w_shape <- SD$ft_t_shape
SD$ft_w_theta <- SD$ft_t_theta
SD$ATE_w <- SD$ATE_t



##### fill the SD #####

options(warn = 2)
for (i in 1:nrow(SD)) {
    
    ## gamma, b, m
    gamma <- SD$gamma[i]
    b <- SD$b[i]
    m <- SD$m[i]
    
    ## pt
    pt <- construct_pt(m, SD$pt_shape[i])
    
    ## taut_t
    taut_t <- construct_taut(m, SD$AvgTau_t[i], SD$tau_t_shape[i])
    
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
    
    ## taut_w
    taut_w <- construct_taut(m, SD$AvgTau_w[i], SD$tau_w_shape[i])
    
    ## gt_w, alpha_w
    gt_w <- construct_ftgt(m, SD$gt_w_shape[i])
    SD$q_w[i] <- ncol(gt_w)
    alpha_w <- solve_alpha(SD$gt_w_shape[i], m, SD$ASPNC_w[i], gt_w, taut_w,
                           theta = SD$gt_w_theta[i])
    SD[i, c("alpha_w0", "alpha_w1", "alpha_w2")[1:SD$q_w[i]]] <- alpha_w
    SD$ASPNC_w_empirical[i] <- compute_ASPNC(gt_w, alpha_w, taut_w)
    stopifnot(all.equal(SD$ASPNC_w[i], SD$ASPNC_w_empirical[i]))
    
    ## ft_w
    ft_w <- construct_ftgt(m, SD$ft_w_shape[i])
    SD$p_w[i] <- ncol(ft_w)
    beta_w <- solve_beta(SD$ft_w_shape[i], m, SD$ATE_w[i], ft_w, gt_w, alpha_w,
                         taut_w, theta = SD$ft_w_theta[i])
    SD[i, c("beta_w0", "beta_w1", "beta_w2")[1:SD$p_w[i]]] <- beta_w
    SD$ATE_w_empirical[i] <- compute_ATE(gt_w, alpha_w, ft_w, beta_w, taut_w)
    stopifnot(all.equal(SD$ATE_w[i], SD$ATE_w_empirical[i]))
    
    ## sample size n
    tryCatch({
        SD$n[i] <- mrt_binary_ss(taut_w, ft_w, gt_w, beta_w, alpha_w, pt, gamma, b,
                                 less_than_10_possible = TRUE)
    },
    error = function(cond){
        message("iteration ", i, ": ", cond)
        return(NA)
    })
    
    if (!is.na(SD$n[i])) {
        ## generate a data set to make sure there are no probabilities outside [0,1]
        dta <- dgm_noDE(SD$n[i], m, ft_t, beta_t, gt_t, alpha_t, taut_t, pt)
    }
}

summary(SD)

##### to debug #####

write.csv(SD, file = "_simulation design.csv")
saveRDS(SD, file = "_simulation design.RDS")



# 2. Simulation (parallel) ------------------------------------------------

rm(list = ls())

options(warn = 0)

nsim <- 2000

# library(tidyverse)
library(rootSolve)
library(MRTSampleSizeBinary)

source("../../../functions/dgm_noDE.R")
source("../../../functions/utility.R")

# functions copied from the "MRTAnalysisBinary" R package (not published)
# make sure to update R package folder in "Carrie Cheng" if revised these functions
source("../../../functions/estimator_EMEE_alwaysCenterA.R")
source("../../../functions/find_change_location.R")
source("../../../functions/get_alpha_beta_from_multiroot_result.R")

SD <- readRDS("_simulation design.RDS")

library(foreach)
library(doRNG)
# library(parallel)

# library(doMC)
# max_cores <- 16
# registerDoMC(min(detectCores() - 1, max_cores))

library(doSNOW)
max_cores <- 24
cl <- makeCluster(min(parallel::detectCores() - 1, max_cores))
registerDoSNOW(cl)

result_df_collected <- data.frame()
result_list_collected <- list()
# result_df_collected <- readRDS("result_df_collected.RDS")
# i <- 1
# for (i in 1:2) {

start_time <- Sys.time()
for (i in 1:nrow(SD)) {
    
    current_time <- Sys.time()
    print(Sys.time())
    print(paste0(round(difftime(current_time, start_time, units = "hours"), 2),
                 " hours has lapsed."))
    cat("i =", i, "/", nrow(SD), "\n")
    
    print(SD[i, ])
    
    ## gamma, b, m
    gamma <- SD$gamma[i]
    b <- SD$b[i]
    m <- SD$m[i]
    
    ## pt
    pt <- construct_pt(m, SD$pt_shape[i])
    
    ## taut_t
    taut_t <- construct_taut(m, SD$AvgTau_t[i], SD$tau_t_shape[i])
    
    ## gt_t, alpha_t
    gt_t <- construct_ftgt(m, SD$gt_t_shape[i])
    alpha_t <- as.numeric(SD[i, c("alpha_t0", "alpha_t1", "alpha_t2")[1:SD$q[i]]])
    
    ## ft_t, beta_t
    ft_t <- construct_ftgt(m, SD$ft_t_shape[i])
    beta_t <- as.numeric(SD[i, c("beta_t0", "beta_t1", "beta_t2")[1:SD$p[i]]])
    
    ## taut_w
    # (not used in the following)
    
    ## gt_w, alpha_w
    gt_w <- construct_ftgt(m, SD$gt_w_shape[i])
    alpha_w <- as.numeric(SD[i, c("alpha_w0", "alpha_w1", "alpha_w2")[1:SD$q_w[i]]])
    
    ## ft_w
    ft_w <- construct_ftgt(m, SD$ft_w_shape[i])
    beta_w <- as.numeric(SD[i, c("beta_w0", "beta_w1", "beta_w2")[1:SD$p_w[i]]])
    
    ## n, p_w, q_w
    n <- SD$n[i]
    p_w <- SD$p_w[i]
    q_w <- SD$q_w[i]
    
    ##### simulation begins #####
    
    beta_t_H1 <- beta_t
    for (truth_hypothesis in c("H1", "H0")) {
        
        if (truth_hypothesis == "H0") {
            beta_t <- rep(0, length(beta_t_H1))
        } else if (truth_hypothesis == "H1") {
            beta_t <- beta_t_H1
        }
        
        set.seed(123)
        
        pb <- txtProgressBar(max = nsim, style = 3)
        progress <- function(n) setTxtProgressBar(pb, n)
        opts <- list(progress = progress)
        
        result <- foreach(isim = 1:nsim, .combine = "c", .options.snow = opts) %dorng% {
            library(rootSolve)
            dta <- dgm_noDE(n, m, ft_t, beta_t, gt_t, alpha_t, taut_t, pt)
            
            dta <- add_analysis_vars(dta, n, ft_w, gt_w)
            control_varname <- get_control_varname(gt_w)
            moderator_varname <- get_moderator_varname(ft_w)
            
            fit <- estimator_EMEE_alwaysCenterA(
                dta = dta,
                id_varname = "id",
                decision_time_varname = "t",
                treatment_varname = "A",
                outcome_varname = "Y",
                control_varname = control_varname, # intercept will be automatically added
                moderator_varname = moderator_varname, # intercept will be automatically added
                rand_prob_varname = "A_prob1",
                avail_varname = "I",
                estimator_initial_value = NULL)
            
            test_stat_unadj <- fit$beta_hat %*%
                solve(fit$varcov[(q_w+1):(q_w+p_w), (q_w+1):(q_w+p_w)]) %*%
                fit$beta_hat
            
            test_stat_adj <- fit$beta_hat %*%
                solve(fit$varcov_adjusted[(q_w+1):(q_w+p_w), (q_w+1):(q_w+p_w)]) %*%
                fit$beta_hat
            
            output <- list(list(fit = fit, test_stat_unadj = test_stat_unadj,
                                test_stat_adj = test_stat_adj))
        }
        close(pb)
        
        test_stat_unadj <- sapply(result, function(l) l$test_stat_unadj)
        test_stat_adj <- sapply(result, function(l) l$test_stat_adj)
        
        cri_val <- p_w * (n - q_w - 1) / (n - q_w - p_w) *
            qf(1 - gamma, df1 = p_w, df2 = n - q_w - p_w)
        
        if (truth_hypothesis == "H1") {
            power_unadj <- sum(test_stat_unadj > cri_val) / nsim
            power_adj <- sum(test_stat_adj > cri_val) / nsim
            power_vec <- c(power_unadj, power_adj)
            names(power_vec) <- c("power_unadj", "power_adj")
            print(power_vec)
        } else if (truth_hypothesis == "H0") {
            typeierror_unadj <- sum(test_stat_unadj > cri_val) / nsim
            typeierror_adj <- sum(test_stat_adj > cri_val) / nsim
            typeierror_vec <- c(typeierror_unadj, typeierror_adj)
            names(typeierror_vec) <- c("typeierror_unadj", "typeierror_adj")
            print(typeierror_vec)
        }
    }
        
    result_df <- cbind(SD[i, ],
                       data.frame(power_unadj = power_unadj,
                                  power_adj = power_adj,
                                  typeierror_unadj = typeierror_unadj,
                                  typeierror_adj = typeierror_adj))
    
    result_df_collected <- rbind(result_df_collected, result_df)
    result_list_collected <- c(result_list_collected, list(list(
        result_df = result_df,
        test_stat_unadj = test_stat_unadj,
        test_stat_adj = test_stat_adj
    )))
    
    saveRDS(result_df_collected, file = paste0("result_df_collected_nsim", nsim, ".RDS"))
    saveRDS(result_list_collected, file = paste0("result_list_collected_nsim", nsim, ".RDS"))
    write.csv(result_df_collected, file = paste0("result_df_collected_nsim", nsim, ".csv"))
    
}

stopCluster(cl) 




# 3. Make plots and tables ------------------------------------------------

if (0) {
    rm(list = ls())
    setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
    library(tidyverse)
    library("gridExtra") 
    source("../../../functions/myggplot.R")
    
    nsim <- 2000
    result <- readRDS(paste0("result_df_collected_nsim", nsim, ".RDS"))
    # result2 <- readRDS(paste0("../simulation2.2.1/result_df_collected_nsim", nsim, ".RDS"))
    # result <- rbind(result, result2)
    # saveRDS(result, paste0("result_df_collected_nsim", nsim, ".RDS"))
    
    hist(result$power_adj)
    hist(result$n)
    
    min_power <- min(c(result$power_adj, result$power_unadj))
    max_power <- max(c(result$power_adj, result$power_unadj))
    min_typeierror <- min(c(result$typeierror_adj, result$typeierror_unadj))
    max_typeierror <- max(c(result$typeierror_adj, result$typeierror_unadj))
    
    power_p1 <- result %>%
        ggplot(aes(x = power_adj)) +
        geom_histogram(aes(y = ..density..), binwidth = 0.002, color = "white") +
        stat_function(fun = dnorm,
                      args = list(mean = 0.8,
                                  sd = sd(result$power_adj)),
                      color = "orange",
                      size = 1.5) +
        xlab("Power (adj)") + 
        coord_cartesian(xlim = c(min_power, max_power)) +
        theme_bw()
    
    power_p2 <- result %>%
        ggplot(aes(x = power_unadj)) +
        geom_histogram(aes(y = ..density..), binwidth = 0.002, color = "white") +
        stat_function(fun = dnorm,
                      args = list(mean = 0.8,
                                  sd = sd(result$power_unadj)),
                      color = "orange",
                      size = 1.5) +
        xlab("Power (unadj)") + 
        coord_cartesian(xlim = c(min_power, max_power)) +
        theme_bw()
    
    typeierror_p1 <- result %>%
        ggplot(aes(x = typeierror_adj)) +
        geom_histogram(aes(y = ..density..), binwidth = 0.002, color = "white") +
        stat_function(fun = dnorm,
                      args = list(mean = 0.05,
                                  sd = sd(result$typeierror_adj)),
                      color = "orange",
                      size = 1.5) +
        xlab("Type I error (adj)") + 
        coord_cartesian(xlim = c(min_typeierror, max_typeierror)) +
        theme_bw()

    typeierror_p2 <- result %>%
        ggplot(aes(x = typeierror_unadj)) +
        geom_histogram(aes(y = ..density..), binwidth = 0.002, color = "white") +
        stat_function(fun = dnorm,
                      args = list(mean = 0.05,
                                  sd = sd(result$typeierror_unadj)),
                      color = "orange",
                      size = 1.5) +
        xlab("Type I error (unadj)") + 
        coord_cartesian(xlim = c(min_typeierror, max_typeierror)) +
        theme_bw()
    
    pdf("plot_simu1.2_power_typeierror.pdf", width = 10, height = 10)
    grid.arrange(power_p1, typeierror_p1, power_p2, typeierror_p2, ncol = 2)
    dev.off()
    
    result %>%
        ggplot(aes(x = n, y = power_adj)) +
        geom_point(alpha = 0.6) +
        geom_smooth(method = "loess") +
        xlab("calculated sample size") +
        ylab("Power (adj)") +
        ggtitle("Blue curve: LOWESS fit with 95% pointwise conf. int.") +
        scale_x_continuous(breaks = c(20, 60, 100, 200, 300)) +
        theme_bw() + myggfont()
    ggsave(paste0("plot_simu1.2_power_scatter.pdf"), width = 10, height = 7.5)
}
