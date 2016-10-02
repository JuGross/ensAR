#' Spread-Adjusted Linear Pool (SLP) of two Normals
#' @export
#' @description Computes the weight and scale parameters of the
#' SLP of two normals from a rolling training period.
#' @param x A vector of observations of a weather quantity.
#' @param pred_one A list with two vectors \code{mean} and \code{sd} of
#' the same length as \code{x} respectively, providing parameters of the
#' normal predictive distribution.
#' @param pred_two A list with two vectors \code{mean} and \code{sd}
#' of the same length as \code{x} respectively, providing parameters of the second
#' normal predictive distribution.
#' @param weight_grid The possible values of the SLP weight.
#' @param scale_grid The possible vlaues of the SLP scale parameter.
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
#'
#'
combineSLP <- function(x, pred_one = list(mean = NULL, sd = NULL),
                       pred_two = pred_one, train = 90,
                       weight_grid = seq(0, 1, 0.1),
                       scale_grid = seq(0.6, 1.4, 0.1)) {
    if (!(length(x) > train)) stop("too less observations")
    mu1 <- pred_one$mean
    sd1 <- pred_one$sd
    mu2 <- pred_two$mean
    sd2 <- pred_two$sd
    train_length <- train
    wt_grid <- weight_grid
    ct_grid <- scale_grid
    allgrid <- expand.grid(weight = wt_grid, scale = ct_grid)

    verification_period <- (train_length + 1):length(mu1)

    roll_combine <- function(verification_day) {
        train_period <- (verification_day - (train_length:1))

        crps_combine <- function(allgrid_row) {
            r_w1 <- allgrid_row[1]
            r_w2 <- 1 - r_w1
            r_ct <- allgrid_row[2]

            r_obs <- x[train_period]
            r_mu1 <- mu1[train_period]
            r_sd1 <- sd1[train_period] * r_ct
            r_mu2 <- mu2[train_period]
            r_sd2 <- sd2[train_period] * r_ct

            crps_out <- crps_norm_mixture(x = r_obs, mu1 = r_mu1,
                                          sd1 = r_sd1, mu2 = r_mu2,
                                          sd2 = r_sd2, w1 = r_w1)

            crps.comb <- mean(crps_out)
            return(crps.comb)
        }
        crps_allgrid <- apply(allgrid, 1, crps_combine)
        allgrid_crps <- cbind(allgrid, crps = crps_allgrid)
        # --
        grid_out <- allgrid_crps[which.min(allgrid_crps$crps), ]
        grid_out <- unlist(grid_out)
        grid_out
    }
    weights_out <- vapply(verification_period, roll_combine, numeric(3))
    weights_out <- (as.data.frame(t(weights_out)))[, 1:2]
    comb_out <- cbind(obs = x[verification_period],
                      mu1 = mu1[verification_period],
        sd1 = sd1[verification_period],
        mu2 = mu2[verification_period],
        sd2 = sd2[verification_period],
        weights_out)
    comb_out
}
#' Predictive Moments of the SLP of two Normals
#' @export
#' @description Computes the predictive mean and predictive standard
#' deviation of the SLP of two normals.
#' @param pred_one A list with two vectors \code{mean} and \code{sd} of
#' the same length as \code{x} respectively, providing parameters of the
#' normal predictive distribution.
#' @param pred_two A list with two vectors \code{mean} and \code{sd}
#' of the same length as \code{x} respectively, providing parameters of the second
#' normal predictive distribution.
#' @param weight_one A vector of weights corresponding to \code{pred_one}.
#' @param scale A vector of scale values.
#' @return A list with two elements \code{mu} (predictive mean) and
#' \code{sd} (predictive standard deviation) of the SLP combination.
#'
#' @author J. Gross, A. Moeller.
#'
slp_moments <- function(pred_one = list(mean = NULL, sd = NULL),
                        pred_two = pred_one,
                        weight_one = rep(0.5, length(pred_one$mean)),
                        scale = rep(1, length(pred_one$sd))) {
    w1 <- weight_one
    w2 <- 1 - weight_one
    s <- scale
    mu1 <- pred_one$mean
    sd1 <- pred_one$sd * s
    mu2 <- pred_two$mean
    sd2 <- pred_two$sd * s
    mu_slp <- w1 * mu1 + w2 * mu2
    var_slp <- w1 * (mu1^2 + sd1^2) + w2 * (mu2^2 + sd2^2) - mu_slp^2
    out <- list(mu = mu_slp, sd = sqrt(var_slp))
    out
}
# --
