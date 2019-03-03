#' Spread-Adjusted Linear Pool (SLP) of two Normals
#' @export
#' @description Computes the weight and scale parameters of the
#' SLP of two normals from a rolling training period.
#' @param x A vector of observations of a weather quantity.
#' @param par_one A list with two vectors \code{mean} and \code{sd} of
#' the same length as \code{x} respectively, providing parameters of the
#' normal predictive distribution.
#' @param par_two A list with two vectors \code{mean} and \code{sd}
#' of the same length as \code{x} respectively, providing parameters of the second
#' normal predictive distribution.
#' @param weight_grid The possible values of the SLP weight.
#' @param scale_grid The possible values of the SLP scale parameter.
#' @param train The length of the training period.
#' @details For each forecast date (all dates except for the first
#' \code{train} dates), the optimal combination of SLP weight and scale
#' parameters is found on a grid determined by \code{weight_grid} and
#' \code{scale_grid}, such that the corresponding
#' (predictive distribution, observation)-pair minimizes the average
#' CRPS with respect to the rolling training period of
#' length \code{train}.
#' @return A data frame with the observation, predictive mean and
#' standard deviation of the two normals, and SLP weight and
#' scale parameter for the forecast period.
#' @references Gneiting T.,  Ranjan R. 2013. Combining predictive
#' distributions. \emph{Electronic Journal of Statistics},
#' \strong{7},  1747--1782.
#'
#' @author J. Gross, A. Moeller.
#' @examples
#' combineSLP(16.1, list(mean = 15, sd = 0.8), list(mean = 18, sd = 1), 0)
#'
#'
combineSLP <- function(x, par_one = list(mean = NULL, sd = NULL), par_two = par_one, 
    train = 90, weight_grid = seq(0, 1, 0.1), scale_grid = seq(0.6, 1.4, 
        0.1)) {
    n <- length(x)
    m <- train
    if (n <= m) 
        stop("too less observations")
    if (!.is_weight(weight_grid)) 
        stop("weight grid not between 0 and 1")
    if (any(scale_grid < 0)) 
        stop("scale grid not positive")
    # input moments
    mu1 <- par_one$mean
    sd1 <- par_one$sd
    mu2 <- par_two$mean
    sd2 <- par_two$sd
    # weight-scale grid
    wt_grid <- weight_grid
    ct_grid <- scale_grid
    allgrid <- expand.grid(weight = wt_grid, scale = ct_grid)
    # verification period
    v_period <- (m + 1):n
    
    roll_combine <- function(v_day) {
        # av. crps for each (two-dimensional) grid point
        train_period <- v_day - (m:1)
        crps_combine <- function(allgrid_row) {
            crps_out <- crps_norm(x = x[train_period], mu1 = mu1[train_period], 
                sd1 = (sd1[train_period] * allgrid_row[2]), mu2 = mu2[train_period], 
                sd2 = (sd2[train_period] * allgrid_row[2]), w1 = allgrid_row[1])
            crps.comb <- mean(crps_out)
            return(crps.comb)
        }
        crps_allgrid <- apply(allgrid, 1, crps_combine)
        allgrid_crps <- cbind(allgrid, crps = crps_allgrid)
        # mimimum av. crps
        grid_out <- allgrid_crps[which.min(allgrid_crps$crps), ]
        grid_out <- unlist(grid_out)
        grid_out
    }
    weights_out <- vapply(v_period, roll_combine, numeric(3))
    weights_out <- (as.data.frame(t(weights_out)))[, 1:2]
    # weights and scale for verification period
    comb_out <- cbind(obs = x[v_period], mu1 = mu1[v_period], sd1 = sd1[v_period], 
        mu2 = mu2[v_period], sd2 = sd2[v_period], weights_out)
    comb_out
}
#' Moments of the SLP of two Normals
#' @export
#' @description Computes the mean and standard
#' deviation of the SLP (spread-adjust linear pool) of two normals.
#' @param par_one A list with two vectors \code{mean} and \code{sd} of
#' the same length as \code{x} respectively, providing parameters of the
#' normal distribution.
#' @param par_two A list with two vectors \code{mean} and \code{sd}
#' of the same length as \code{x} respectively, providing parameters of the second
#' normal distribution.
#' @param weight_one A vector of weights corresponding to \code{par_one}.
#' @param scale A vector of scale values.
#' @return A list with two elements \code{mean} and
#' \code{sd} of the SLP combination.
#' @examples
#' slp_moments(list(mean = 15, sd = 0.8), list(mean = 18, sd = 1), 0.7, 1.2)
#' @author J. Gross, A. Moeller.
#'
slp_moments <- function(par_one = list(mean = NULL, sd = NULL), par_two = par_one, 
    weight_one = rep(0.5, length(par_one$mean)), scale = rep(1, length(par_one$sd))) {
    w1 <- weight_one
    if (any((w1 < 0) | (w1 > 1))) 
        stop("weights must lie between 0 and 1")
    w2 <- 1 - weight_one
    s <- scale
    if (any(s < 0)) 
        stop("scale must be positive")
    mu1 <- par_one$mean
    sd1 <- (par_one$sd * s)
    mu2 <- par_two$mean
    sd2 <- (par_two$sd * s)
    mu_slp <- w1 * mu1 + w2 * mu2
    sd_slp <- sqrt(w1 * (mu1^2 + sd1^2) + w2 * (mu2^2 + sd2^2) - mu_slp^2)
    out <- list(mean = mu_slp, sd = sd_slp)
    out
}
# --
