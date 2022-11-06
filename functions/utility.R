### utility functions ###
# Binary MRT sample size paper
#
# Tianchen Qian
# 2022.05.07


##### construct ft, gt, taut, pt from patterns
construct_ftgt <- function(m, shape) {
    # linear_theta is for ft and beta, see goodnotes "dgm_noDE"
    match.arg(shape, c("constant", "linear", "quadratic", 
                       "linear_theta", "quadratic_theta"))
    if (shape == "constant") {
        return(matrix(1, nrow = m, ncol = 1))
    } else if (shape %in% c("linear", "linear_theta")) {
        return(cbind(1, 1:m))
    } else if (shape %in% c("quadratic", "quadratic_theta")) {
        return(cbind(1, 1:m, (1:m)^2))
    }
}

construct_taut <- function(m, AvgTau, shape) {
    match.arg(shape, c("constant", "linear pm 0.1", "linear pm 0.2",
                       "sin pm 0.1", "sin pm 0.2"))
    if (shape == "constant") {
        taut <- rep(AvgTau, m)
    } else if (shape == "linear pm 0.1") {
        taut <- seq(from = AvgTau + 0.1, to = AvgTau - 0.1, length.out = m)
    } else if (shape == "linear pm 0.2") {
        taut <- seq(from = AvgTau + 0.2, to = AvgTau - 0.2, length.out = m)
    } else if (shape == "sin pm 0.1") {
        taut <- AvgTau + 0.1 * sin(1:m)
    } else if (shape == "sin pm 0.2") {
        taut <- AvgTau + 0.2 * sin(1:m)
    }
    stopifnot(all(taut >= 0 & taut <= 1))
    return(taut)
}


construct_taut_theta <- function(m, AvgTau, shape, theta = 0) {
    match.arg(shape, c("constant", "linear", "sine"))
    stopifnot(AvgTau + theta <= 1 & AvgTau + theta >= 0 &
                  AvgTau - theta <= 1 & AvgTau - theta >= 0)
    if (shape == "constant") {
        taut <- rep(AvgTau, m)
    } else if (shape == "linear") {
        taut <- seq(from = AvgTau + theta, to = AvgTau - theta, length.out = m)
    } else if (shape == "sine") {
        taut <- AvgTau + theta * sin(1:m)
    }
    stopifnot(all(taut >= 0 & taut <= 1))
    return(taut)
}


construct_pt <- function(m, shape) {
    match.arg(shape, c("constant 0.5", "constant 0.6", "linear 0.5 pm 0.1"))
    if (shape == "constant 0.5") {
        return(rep(0.5, m))
    } else if (shape == "constant 0.6") {
        return(rep(0.6, m))
    } else if (shape == "linear 0.5 pm 0.1") {
        return(seq(from = 0.6, to = 0.4, length.out = m))
    }
}


##### for use in appending ft and gt to data set #####
add_analysis_vars <- function(dta, n, ft, gt) {
    stopifnot(nrow(ft) * n == nrow(dta))
    stopifnot(nrow(gt) * n == nrow(dta))
    
    stopifnot(all(ft[, 1] == 1))
    stopifnot(all(gt[, 1] == 1))
    
    if (ncol(ft) >= 2) {
        for (ivar in 1:(ncol(ft) - 1)) {
            dta[, paste0("mod", ivar)] <- rep(ft[, ivar + 1], n)
        }
    }
    if (ncol(gt) >= 2) {
        for (ivar in 1:(ncol(gt) - 1)) {
            dta[, paste0("ctr", ivar)] <- rep(gt[, ivar + 1], n)
        }
    }
    
    return(dta)
}

get_control_varname <- function(gt) {
    if (ncol(gt) == 1) {
        return(NULL)
    } else {
        return(paste0("ctr", 1:(ncol(gt) - 1)))
    }
}

get_moderator_varname <- function(ft) {
    if (ncol(ft) == 1) {
        return(NULL)
    } else {
        return(paste0("mod", 1:(ncol(ft) - 1)))
    }
}