### dgm without delayed effect ###
# Binary MRT sample size paper
#
# Tianchen Qian
# 2022.05.07

# a faster way to generate data
dgm_noDE <- function(n, m, ft, beta, gt, alpha, taut, pt) {
    # n: sample size
    # m: number of decision points per individual
    
    df_names <- c("id", "t", "I", "A", "Y",
                  "I_prob1", "A_prob1", "Y_prob1")
    
    dta <- data.frame(matrix(NA, nrow = n * m, 
                             ncol = length(df_names)))
    names(dta) <- df_names
    
    dta$id <- rep(1:n, each = m)
    dta$t <- rep(1:m, times = n)
    
    dta$I_prob1 <- rep(taut, n)
    dta$I <- rbinom(n * m, 1, dta$I_prob1)
    dta$A_prob1 <- rep(pt, n) * dta$I
    dta$A <- rbinom(n * m, 1, dta$A_prob1)
    dta$Y_prob1 <- exp(dta$A * rep(as.vector(ft %*% beta), n) + 
                           rep(as.vector(gt %*% alpha), n))
    if (!(all(dta$Y_prob1 >= 0 & dta$Y_prob1 <= 1))) {
        message("min(dta$Y_prob1) = ", min(dta$Y_prob1), 
                "; max(dta$Y_prob1) = ", max(dta$Y_prob1))
    }
    stopifnot(all(dta$Y_prob1 >= 0 & dta$Y_prob1 <= 1))
    dta$Y <- rbinom(n * m, 1, dta$Y_prob1)
    
    return(dta)    
}

dgm_noDE2 <- function(n, m, ft, beta, gt, alpha, taut, pt) {
    # n: sample size
    # m: number of decision points per individual
    
    df_names <- c("id", "t", "I", "A", "Y",
                  "I_prob1", "A_prob1", "Y_prob1")
    
    dta <- data.frame(matrix(NA, nrow = n * m, 
                             ncol = length(df_names)))
    names(dta) <- df_names
    
    dta$id <- rep(1:n, each = m)
    dta$t <- rep(1:m, times = n)
    
    for (i in 1:n) {
        rows <- ((i-1) * m + 1) : (i * m)
        dta$I_prob1[rows] <- taut
        dta$I[rows] <- rbinom(m, 1, dta$I_prob1[rows])
        dta$A_prob1[rows] <- pt * dta$I[rows]
        dta$A[rows] <- rbinom(m, 1, dta$A_prob1[rows])
        dta$Y_prob1[rows] <- exp(dta$A[rows] * ft %*% beta) * exp(gt %*% alpha)
        dta$Y[rows] <- rbinom(m, 1, dta$Y_prob1[rows])
    }
    
    return(dta)    
}


##### compute empirical quantities #####

compute_ATE <- function(gt, alpha, ft, beta, taut) {
    num <- sum(exp(gt %*% alpha + ft %*% beta) * taut)
    denom <- sum(exp(gt %*% alpha) * taut)
    return(num / denom)
}

compute_ASPNC <- function(gt, alpha, taut) {
    num <- sum(exp(gt %*% alpha) * taut)
    denom <- sum(taut)
    return(num / denom)
}

compute_SPNC <- function(gt, alpha, taut) {
    num <- exp(gt %*% alpha) * taut
    denom <- taut
    return(num / denom)
}


##### solve beta #####

solve_beta <- function(shape, m, ATE, ft, gt, alpha, taut, ...) {
    match.arg(shape, c("constant", "linear", "quadratic",
                       "linear_theta", "quadratic_theta"))
    if (shape == "constant") {
        return(solve_beta_constant(m, ATE, ft, gt, alpha, taut))
    } else if (shape == "linear") {
        return(solve_beta_linear(m, ATE, ft, gt, alpha, taut))
    } else if (shape == "quadratic") {
        return(solve_beta_quadratic(m, ATE, ft, gt, alpha, taut))
    } else if (shape == "linear_theta") {
        return(solve_beta_linear_theta(m, ATE, ft, gt, alpha, taut, ...))
    } else if (shape == "quadratic_theta") {
        return(solve_beta_quadratic_theta(m, ATE, ft, gt, alpha, taut, ...))
    }
}

solve_beta_constant <- function(m, ATE, ft, gt, alpha, taut) {
    return(log(ATE))
}

solve_beta_linear <- function(m, ATE, ft, gt, alpha, taut) {
    # effect at last decision point is 0
    stopifnot(ncol(ft) == 2)
    from_beta1_to_ATE <- function(beta1) {
        beta0 <- -m * beta1
        beta <- c(beta0, beta1)
        return(compute_ATE(gt, alpha, ft, beta, taut) - ATE)
    }
    beta1 <- uniroot(from_beta1_to_ATE, interval = c(-0.5, 0.5), extendInt = "yes",
                     tol = 1e-12)$root
    beta0 <- -m * beta1
    return(c(beta0, beta1))
}

solve_beta_quadratic <- function(m, ATE, ft, gt, alpha, taut) {
    # effect at first and last decision points is 0
    stopifnot(ncol(ft) == 3)
    from_beta2_to_ATE <- function(beta2) {
        beta <- c(m * beta2, -(m+1) * beta2, beta2)
        return(compute_ATE(gt, alpha, ft, beta, taut) - ATE)
    }
    beta2 <- uniroot(from_beta2_to_ATE, interval = c(-0.5, 0.5), extendInt = "yes",
                     tol = 1e-12)$root
    return(c(m * beta2, -(m+1) * beta2, beta2))
}


solve_beta_linear_theta <- function(m, ATE, ft, gt, alpha, taut, theta) {
    # see goodnotes "dgm_noDE"
    stopifnot(ncol(ft) == 2)
    stopifnot(theta >= -1 & theta <= 1)
    
    if (theta == 1) {
        from_beta1_to_ATE <- function(beta1) {
            beta0 <- -m * beta1
            beta <- c(beta0, beta1)
            return(compute_ATE(gt, alpha, ft, beta, taut) - ATE)
        }
        beta1 <- uniroot(from_beta1_to_ATE, interval = c(-0.5, 0.5), extendInt = "yes",
                         tol = 1e-12)$root
        beta0 <- -m * beta1
        return(c(beta0, beta1))
    } else if (theta == 0) {
        beta1 <- 0
        from_beta0_to_ATE <- function(beta0) {
            beta <- c(beta0, 0)
            return(compute_ATE(gt, alpha, ft, beta, taut) - ATE)
        }
        beta0 <- uniroot(from_beta0_to_ATE, interval = c(-0.5, 0.5), extendInt = "yes",
                         tol = 1e-12)$root
        return(c(beta0, beta1))
    } else {
        ratio <- (1+theta) / (1-theta)
        from_beta1_to_ATE <- function(beta1) {
            beta0 <- (ratio * m - 1) / (1 - ratio) * beta1
            beta <- c(beta0, beta1)
            return(compute_ATE(gt, alpha, ft, beta, taut) - ATE)
        }
        beta1 <- uniroot(from_beta1_to_ATE, interval = c(-0.5, 0.5), extendInt = "yes",
                         tol = 1e-12)$root
        beta0 <- (ratio * m - 1) / (1 - ratio) * beta1
        return(c(beta0, beta1))
    }
}


solve_beta_quadratic_theta <- function(m, ATE, ft, gt, alpha, taut, theta) {
    # see goodnotes "dgm_noDE"
    stopifnot(ncol(ft) == 3)
    stopifnot(theta >= -1 & theta <= 1)
    if (theta == 1) {
        from_beta2_to_ATE <- function(beta2) {
            beta0 <- m * beta2
            beta1 <- -(m+1) * beta2
            beta <- c(beta0, beta1, beta2)
            return(compute_ATE(gt, alpha, ft, beta, taut) - ATE)
        }
        beta2 <- uniroot(from_beta2_to_ATE, interval = c(-0.5, 0.5), extendInt = "yes",
                         tol = 1e-12)$root
        beta0 <- m * beta2
        beta1 <- -(m+1) * beta2
        return(c(beta0, beta1, beta2))
    } else if (theta == 0) {
        beta1 <- 0
        beta2 <- 0
        from_beta0_to_ATE <- function(beta0) {
            beta <- c(beta0, 0, 0)
            return(compute_ATE(gt, alpha, ft, beta, taut) - ATE)
        }
        beta0 <- uniroot(from_beta0_to_ATE, interval = c(-0.5, 0.5), extendInt = "yes",
                         tol = 1e-12)$root
        return(c(beta0, beta1, beta2))
    } else {
        ratio <- (1+theta) / (1-theta)
        from_beta2_to_ATE <- function(beta2) {
            beta0 <- ((m+1)^2/4 - ratio*(m+1)+ratio) / (1-ratio) * beta2
            beta1 <- -(m+1) * beta2
            beta <- c(beta0, beta1, beta2)
            return(compute_ATE(gt, alpha, ft, beta, taut) - ATE)
        }
        beta2 <- uniroot(from_beta2_to_ATE, interval = c(-0.5, 0.5), extendInt = "yes",
                         tol = 1e-12)$root
        beta0 <- ((m+1)^2/4 - ratio*(m+1)+ratio) / (1-ratio) * beta2
        beta1 <- -(m+1) * beta2
        return(c(beta0, beta1, beta2))
    }
}

##### solve alpha #####

solve_alpha <- function(shape, m, ASPNC, gt, taut, ...) {
    match.arg(shape, c("constant", "linear", "quadratic",
                       "linear_theta", "quadratic_theta"))
    if (shape == "constant") {
        return(solve_alpha_constant(m, ASPNC, gt, taut))
    } else if (shape == "linear") {
        return(solve_alpha_linear(m, ASPNC, gt, taut, ...))
    } else if (shape == "quadratic") {
        return(solve_alpha_quadratic(m, ASPNC, gt, taut, ...))
    } else if (shape == "linear_theta") {
        return(solve_alpha_linear_theta(m, ASPNC, gt, taut, ...))
    } else if (shape == "quadratic_theta") {
        return(solve_alpha_quadratic_theta(m, ASPNC, gt, taut, ...))
    }
}

solve_alpha_constant <- function(m, ASPNC, gt, taut) {
    return(log(ASPNC))
}

solve_alpha_linear <- function(m, ASPNC, gt, taut, MIF) {
    # MIF: max inflation factor
    # SPNC at first decision point = ASPNC * MIF
    stopifnot(ncol(gt) == 2)
    from_alpha1_to_ASPNC <- function(alpha1) {
        alpha <- c(log(ASPNC * MIF) - alpha1,
                   alpha1)
        return(compute_ASPNC(gt, alpha, taut) - ASPNC)
    }
    alpha1 <- uniroot(from_alpha1_to_ASPNC, interval = c(-0.5, 0.5), extendInt = "yes",
                      tol = 1e-12)$root
    return(c(log(ASPNC * MIF) - alpha1,
             alpha1))
}

solve_alpha_quadratic <- function(m, ASPNC, gt, taut, MIF) {
    # MIF: max inflation factor
    # SPNC at first and last decision points are equal
    # SPNC at max (t = (m+1)/2) is ASPNC * MIF
    stopifnot(ncol(gt) == 3)
    from_alpha2_to_ASPNC <- function(alpha2) {
        alpha <- c((m+1)^2/4 * alpha2 + log(ASPNC * MIF),
                   -(m+1) * alpha2,
                   alpha2)
        return(compute_ASPNC(gt, alpha, taut) - ASPNC)
    }
    alpha2 <- uniroot(from_alpha2_to_ASPNC, interval = c(-0.5, 0.5), extendInt = "yes",
                      tol = 1e-12)$root
    return(c((m+1)^2/4 * alpha2 + log(ASPNC * MIF),
             -(m+1) * alpha2,
             alpha2))
}

solve_alpha_linear_theta <- function(m, ASPNC, gt, taut, theta) {
    # see goodnotes "dgm_noDE" solve_beta_linear_theta
    # the math for alpha and beta are exactly the same
    stopifnot(ncol(gt) == 2)
    stopifnot(theta >= -1 & theta <= 1)
    
    if (theta == 1) {
        from_alpha1_to_ASPNC <- function(alpha1) {
            alpha0 <- -m * alpha1
            alpha <- c(alpha0, alpha1)
            return(compute_ASPNC(gt, alpha, taut) - ASPNC)
        }
        alpha1 <- uniroot(from_alpha1_to_ASPNC, interval = c(-0.5, 0.5), extendInt = "yes",
                          tol = 1e-12)$root
        alpha0 <- -m * alpha1
        return(c(alpha0, alpha1))
    } else if (theta == 0) {
        alpha1 <- 0
        from_alpha0_to_ASPNC <- function(alpha0) {
            alpha <- c(alpha0, 0)
            return(compute_ASPNC(gt, alpha, taut) - ASPNC)
        }
        alpha0 <- uniroot(from_alpha0_to_ASPNC, interval = c(-0.5, 0.5), extendInt = "yes",
                          tol = 1e-12)$root
        return(c(alpha0, alpha1))
    } else {
        ratio <- (1+theta) / (1-theta)
        from_alpha1_to_ASPNC <- function(alpha1) {
            alpha0 <- (ratio * m - 1) / (1 - ratio) * alpha1
            alpha <- c(alpha0, alpha1)
            return(compute_ASPNC(gt, alpha, taut) - ASPNC)
        }
        alpha1 <- uniroot(from_alpha1_to_ASPNC, interval = c(-0.5, 0.5), extendInt = "yes",
                          tol = 1e-12)$root
        alpha0 <- (ratio * m - 1) / (1 - ratio) * alpha1
        return(c(alpha0, alpha1))
    }
}


solve_alpha_quadratic_theta <- function(m, ASPNC, gt, taut, theta) {
    # see goodnotes "dgm_noDE" solve_beta_quadratic_theta
    # the math for alpha and beta are exactly the same
    stopifnot(ncol(gt) == 3)
    stopifnot(theta >= -1 & theta <= 1)
    
    if (theta == 1) {
        from_alpha2_to_ASPNC <- function(alpha2) {
            alpha0 <- m * alpha2
            alpha1 <- -(m+1) * alpha2
            alpha <- c(alpha0, alpha1, alpha2)
            return(compute_ASPNC(gt, alpha, taut) - ASPNC)
        }
        alpha2 <- uniroot(from_alpha2_to_ASPNC, interval = c(-0.5, 0.5), extendInt = "yes",
                          tol = 1e-12)$root
        alpha0 <- m * alpha2
        alpha1 <- -(m+1) * alpha2
        return(c(alpha0, alpha1, alpha2))
    } else if (theta == 0) {
        alpha1 <- 0
        alpha2 <- 0
        from_alpha0_to_ASPNC <- function(alpha0) {
            alpha <- c(alpha0, 0, 0)
            return(compute_ASPNC(gt, alpha, taut) - ASPNC)
        }
        alpha0 <- uniroot(from_alpha0_to_ASPNC, interval = c(-0.5, 0.5), extendInt = "yes",
                          tol = 1e-12)$root
        return(c(alpha0, alpha1, alpha2))
    } else {
        ratio <- (1+theta) / (1-theta)
        from_alpha2_to_ASPNC <- function(alpha2) {
            alpha0 <- ((m+1)^2/4 - ratio*(m+1)+ratio) / (1-ratio) * alpha2
            alpha1 <- -(m+1) * alpha2
            alpha <- c(alpha0, alpha1, alpha2)
            return(compute_ASPNC(gt, alpha, taut) - ASPNC)
        }
        alpha2 <- uniroot(from_alpha2_to_ASPNC, interval = c(-1, 1), extendInt = "yes",
                          tol = 1e-12)$root
        alpha0 <- ((m+1)^2/4 - ratio*(m+1)+ratio) / (1-ratio) * alpha2
        alpha1 <- -(m+1) * alpha2
        return(c(alpha0, alpha1, alpha2))
    }
}