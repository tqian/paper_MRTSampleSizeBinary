### dgm with delayed effect on Y###
# Binary MRT sample size paper
#
# Tianchen Qian
# 2022.10.06


dgm_SerialCorrelationOnY <- function(n, m, ft, beta, gt, alpha, taut, pt, coef_SC) {
    # n: sample size
    # m: number of decision points per individual
    # coef_SC: coefficient for how Y_{t+1} depends on Y_t
    
    if (length(unique(pt)) != 1) {
        stop("pt is not constant over time. This will make ATE and ASPNC calculation incorrect.")
    }
    if (!all(taut == 1)) {
        stop("taut is not all 1. This will make the calculation of SPNC(t) incorrect.")
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
            Y_lag1 <- rep(0, n)
        }
        dta$Y_prob1[rows] <- exp(as.numeric(gt[t, ] %*% alpha) + 
                                     dta$A[rows] * as.numeric(ft[t, ] %*% beta) + 
                                     Y_lag1 * coef_SC)
        dta$Y[rows] <- rbinom(n, 1, dta$Y_prob1[rows])
        Y_lag1 <- dta$Y[rows]
    }
    if (!(all(dta$Y_prob1 >= 0 & dta$Y_prob1 <= 1))) {
        message("min(dta$Y_prob1) = ", min(dta$Y_prob1), 
                "; max(dta$Y_prob1) = ", max(dta$Y_prob1))
        stop("some Y_prob1 is less than 0 or greater than 1")
    }
    return(dta)    
}


##### compute empirical quantities #####

compute_ATE_SConY <- function(gt, alpha, ft, beta, taut) {
    num <- sum(exp(gt %*% alpha + ft %*% beta) * taut)
    denom <- sum(exp(gt %*% alpha) * taut)
    return(num / denom)
}

compute_ASPNC_SConY <- function(SPNC) {
    return(mean(SPNC))
}

compute_SPNC_SConY <- function(gt, alpha, ft, beta, taut, coef_SC, pt) {
    m <- nrow(gt)
    SPNC <- rep(NA, m)
    SPNC[1] <- exp(as.numeric(gt[1,] %*% alpha))
    for (t in 2:m) {
        SPNC[t] <- exp(as.numeric(gt[t,] %*% alpha)) *
            (1 + (exp(coef_SC) - 1) * SPNC[t-1] * 
                 (1 - pt[t-1] + pt[t-1] * exp(as.numeric(ft[t-1,] %*% beta))))
    }
    return(SPNC)
}


##### solve beta #####
# These are the same as for dgm_noDE

solve_beta_SConY <- function(shape, m, ATE, ...) {
    stopifnot(shape == "constant")
    solve_beta(shape, m, ATE, ...)
}


##### solve alpha #####
# Use the same function as dgm_noDE
# This doesn't give us the correct ASPN, but we will use this anyway,
# and report the true ASPN calculated using compute_ASPNC_SConY()

solve_alpha_SConY <- function(shape, m, ASPNC, gt, taut, ...) {
    solve_alpha(shape, m, ASPNC, gt, taut, ...)
}