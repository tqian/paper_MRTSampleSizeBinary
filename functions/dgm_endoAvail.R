### dgm with endogenous availability process ###
# Binary MRT sample size paper
#
# Tianchen Qian
# 2022.10.06


dgm_endoAvail <- function(n, m, ft, beta, gt, alpha, pt, gamma3, gamma4) {
    # n: sample size
    # m: number of decision points per individual
    # gamma3: coefficient of A_{t-1} on I_t
    # gamma4: coefficient of Y_t on I_t
    
    stopifnot(0.5 + gamma3 + gamma4 < 1)
    stopifnot(0.5 + gamma3 + gamma4 > 0)
    
    df_names <- c("id", "t", "I", "A", "Y",
                  "I_prob1", "A_prob1", "Y_prob1")
    
    dta <- data.frame(matrix(NA, nrow = n * m, 
                             ncol = length(df_names)))
    names(dta) <- df_names
    
    dta$id <- rep(1:n, each = m)
    dta$t <- rep(1:m, times = n)

    for (t in 1:m) {
        if (t == 1) {
            A_lag1 <- rep(0, n)
            Y_lag1 <- rep(0, n)
        }
        rows <- seq(from = t, length.out = n, by = m)
        dta$I_prob1[rows] <- 0.5 + gamma3 * A_lag1 + gamma4 * Y_lag1
        dta$I[rows] <- rbinom(n, 1, dta$I_prob1[rows])
        dta$A_prob1[rows] <- rep(pt[t], n) * dta$I[rows]
        dta$A[rows] <- rbinom(n, 1, dta$A_prob1[rows])
        dta$Y_prob1[rows] <- exp(as.numeric(gt[t, ] %*% alpha) + 
                                     dta$A[rows] * as.numeric(ft[t, ] %*% beta))
        dta$Y[rows] <- rbinom(n, 1, dta$Y_prob1[rows])
        A_lag1 <- dta$A[rows]
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


compute_tau_endoAvail <- function(gt, alpha, ft, beta, pt, gamma3, gamma4) {
    m <- nrow(gt)
    tau <- rep(NA, m)
    tau[1] <- 0.5
    for (t in 2:m) {
        tau[t] <- 0.5 + gamma3 * pt[t-1] + 
            (exp(as.numeric(ft[t,] %*% beta)) * pt[t-1] + 1 - pt[t-1]) *
            gamma4 * exp(as.numeric(gt[t,] %*% alpha)) * tau[t-1]
    }
    return(tau)
}

# 
# ##### solve beta #####
# # These are the same as for dgm_noDE
# 
# solve_beta_SConY <- function(shape, m, ATE, ...) {
#     stopifnot(shape == "constant")
#     solve_beta(shape, m, ATE, ...)
# }
# 
# 
# ##### solve alpha #####
# # Use the same function as dgm_noDE
# # This doesn't give us the correct ASPN, but we will use this anyway,
# # and report the true ASPN calculated using compute_ASPNC_SConY()
# 
# solve_alpha_SConY <- function(shape, m, ASPNC, gt, taut, ...) {
#     solve_alpha(shape, m, ASPNC, gt, taut, ...)
# }