### dgm with delayed effect on Y###
# Binary MRT sample size paper
#
# Tianchen Qian
# 2022.05.07

### WARNING ###
# ATE (and hence solver for beta) is computed the same way as for dgm_noDE(),
# as long as:
# (1) A_{t-1} doesn't depend on A_t
# (2) distribution of A_t is time-homogenous
# 
# ASPNC (and hence solver for alpha) is computed differently as it depends on coef_DE.
# 
# See "DGM" in GoodNotes for details.

dgm_DEonY <- function(n, m, ft, beta, gt, alpha, taut, pt, coef_DE) {
    # n: sample size
    # m: number of decision points per individual
    # coef_DE: coefficient for delayed effect
    
    if (length(unique(pt)) != 1) {
        stop("pt is not constant over time. This will make ATE and ASPNC calculation incorrect.")
    }
    
    df_names <- c("id", "t", "I", "A", "Y",
                  "I_prob1", "A_prob1", "Y_prob1")
    
    dta <- data.frame(matrix(NA, nrow = n * m, 
                             ncol = length(df_names)))
    names(dta) <- df_names
    
    dta$id <- rep(1:n, each = m)
    dta$t <- rep(1:m, times = n)

    for (t in 1:m) {
        rows <- seq(from = t, length.out = n, by = m)
        dta$I_prob1[rows] <- rep(taut[t], n)
        dta$I[rows] <- rbinom(n, 1, dta$I_prob1[rows])
        dta$A_prob1[rows] <- rep(pt[t], n) * dta$I[rows]
        dta$A[rows] <- rbinom(n, 1, dta$A_prob1[rows])
        if (t == 1) {
            A_lag1 <- rbinom(n, 1, pt[1])
        }
        dta$Y_prob1[rows] <- exp(as.numeric(gt[t, ] %*% alpha) + 
                                     dta$A[rows] * as.numeric(ft[t, ] %*% beta) + 
                                     A_lag1 * coef_DE)
        if (!(all(dta$Y_prob1[rows] >= 0 & dta$Y_prob1[rows] <= 1))) {
            message("min(dta$Y_prob1[rows]) = ", min(dta$Y_prob1[rows]), 
                    "; max(dta$Y_prob1[rows]) = ", max(dta$Y_prob1[rows]))
        }
        stopifnot(all(dta$Y_prob1[rows] >= 0 & dta$Y_prob1[rows] <= 1))
        dta$Y[rows] <- rbinom(n, 1, dta$Y_prob1[rows])
        A_lag1 <- dta$A[rows]
    }
    
    return(dta)    
}


##### compute empirical quantities #####

compute_ATE_DEonY <- function(gt, alpha, ft, beta, taut) {
    num <- sum(exp(gt %*% alpha + ft %*% beta) * taut)
    denom <- sum(exp(gt %*% alpha) * taut)
    return(num / denom)
}

compute_ASPNC_DEonY <- function(gt, alpha, taut, coef_DE, pt) {
    pt_lag1 <- c(pt[1], pt[1:(length(pt) - 1)])
    num <- sum(exp(gt %*% alpha) * (exp(coef_DE) * pt_lag1 + (1 - pt_lag1)) * taut)
    denom <- sum(taut)
    return(num / denom)
}

compute_SPNC_DEonY <- function(gt, alpha, taut, coef_DE, pt) {
    pt_lag1 <- c(pt[1], pt[1:(length(pt) - 1)])
    num <- exp(gt %*% alpha) * (exp(coef_DE) * pt_lag1 + (1 - pt_lag1)) * taut
    denom <- taut
    return(num / denom)
}


##### solve beta #####
# These are the same as for dgm_noDE

solve_beta_DEonY <- function(shape, m, ATE, ft, gt, alpha, taut, ...) {
    solve_beta(shape, m, ATE, ft, gt, alpha, taut, ...)
}


##### solve alpha #####
# These are almost dgm_noDE but with an offset on beta0
# see goodnotes "dgm_DEonY"

solve_alpha_DEonY <- function(shape, m, ASPNC, gt, taut, coef_DE, pt, ...) {
    alpha <- solve_alpha(shape, m, ASPNC, gt, taut, ...)
    pt1 <- pt[1]
    offset <- (-1) * log(exp(coef_DE) * pt1 + 1 - pt1)
    alpha[1] <- alpha[1] + offset
    return(alpha)
}