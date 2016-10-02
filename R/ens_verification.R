#' CRPS for a Normal or a Mixture of two Normals
#' @export
#' @importFrom stats pnorm dnorm
#' @description Computes the continuous ranked probability score
#' (CRPS) of given observations when the predictive distribution is
#' the normal distribution or a mixture of two normals.
#'
#' @param x A vector of observations for which the CRPS is to be computed.
#' @param mu1 A vector of expectations of the first normal. Must be of the same length as \code{x}.
#' @param sd1 A vector of standard deviations of the first normal. Must be of the same length as \code{x}.
#' @param mu2 A vector of expectations of the second normal. Must be of the same length as \code{x}.
#' @param sd2 A vector of standard deviations of the second normal. Must be of the same length as \code{x}.
#' @param w1 A vector of weights between 0 and 1 associated with the first normal.
#' Must be either of length one, or of the same length as \code{x}.
#' @return A vector of CRPS values.
#' @details Formula (5) from Grimit et al. (2006) is applied
#' for the special case of two normals.
#' @author J. Gross, A. Moeller.
#' @references Grimit E.P., Gneiting T., Berrocal V., Johnson N.A.
#'2006. The continuous ranked probability score for circular
#'variables and its application to mesoscale forecast ensemble
#'verification. \emph{Quarterly Journal of the Royal Meteorological
#'Society}, \strong{132}, 2925--2942.
#'@examples
#'crps_norm(0.5, -1, 1, 2, 1, 0.2)
crps_norm <- function(x, mu1, sd1, mu2 = mu1, sd2 = sd1, w1 = 1) {
    A <- function(m, s) {
        z <- m/s
        2 * s * dnorm(z) + m * (2 * pnorm(z) - 1)
    }
    w2 <- 1 - w1
    val_1 <- w1 * A((x - mu1), sd1)
    if (w2 == 0) {val <- val_1} else {
    val <- val_1 + w2 * A((x - mu2), sd2) -
      0.5 * ((w1^2) * A(0, sqrt(2 * sd1^2)) +
               w1 * w2 * A((mu1 - mu2), sqrt(sd1^2 + sd2^2)) +
               w2 * w1 * A((mu2 - mu1), sqrt(sd1^2 + sd2^2)) +
               (w2^2) * A(0, sqrt(2 * sd2^2)))
    }
    val
}
#' Verification Statistics
#' @export
#' @importFrom stats pnorm dnorm var
#' @description Computes different verification statistics for
#' one predictive normal distribution, or for the spread-adjusted
#' linear pool (SLP) of two predictive normal distributions.
#' @param x A vector of observations of a weather quantity.
#' @param pred_one A list with two vectors \code{mean} and \code{sd} of
#' the same length as \code{x} respectively, providing parameters of the
#' normal predictive distribution.
#' @param pred_two A list with two vectors \code{mean} and \code{sd}
#' of the same length as \code{x} respectively, providing parameters of the second
#' normal predictive distribution. Only required if \code{weight_one}
#' is not \code{NULL} in which case the two parameter vectors must
#' be of the same length as \code{x}.
#' @param weight_one A vector of the same length as \code{x} providing
#' the weights (corresponding to the first normal) used for the spread-adjusted linear pool.
#' @param scale A vector of the same length as \code{x} providing the common scale
#' parameter used for the spread-adjusted linear pool.
#' @details If \code{weight_one} is \code{NULL} verification statistics
#' are computed for a normal distribution with parameters provided by
#' \code{pred_one}. Otherwise the predictive distribution is the
#' spread-adjusted linear pool (SLP) of two normals with parameters
#' provided by \code{pred_one} and \code{pred_two} and
#' SLP parameters \code{weight_one} (weights
#' corresponding to \code{pred_one}) and common scale
#' parameter \code{scale}.
#' @return A list with elements
#' \describe{
#' \item{\code{pit}}{A vector of the same length as \code{x}
#' giving probability integral transform (PIT) values.}
#' \item{\code{crps}}{A vector of the same length as \code{x} giving
#' the continuous ranked probability score (CRPS) values.}
#'\item{\code{dss}}{A vector of the same length as \code{x} giving
#' the Dawid and Sebastiani score (DSS) values.}
#' \item{\code{pitvar}}{The sample variance of the PIT values.}
#' \item{\code{rmv}}{The square root of the mean of the predictive variance.}
#' }
#' @examples
#' veri_stats(17.5, pred_one = list(mean = 15, sd = 1))
#' @author J. Gross, A. Moeller.
veri_stats <- function(x, pred_one = list(mean = NULL, sd = NULL),
                       pred_two = pred_one, weight_one = NULL,
                       scale = NULL) {
    if (is.null(weight_one)) {
        mu <- pred_one$mean
        sd <- pred_one$sd
        pit <- pnorm(x, mu, sd)
        z <- (x - mu)/sd
        dss <- z^2 + 2 * log(sd)
        crps <- sd * (z * (2 * pnorm(z) - 1) + 2 * dnorm(z) - 1/sqrt(pi))
        predvariance <- sd^2
    } else {
        # slp combination
        w1 <- weight_one
        w2 <- 1 - weight_one
        if (is.null(scale)) {
            s <- rep(1, length = length(w1))
        } else {
            s <- scale
        }
        mu1 <- pred_one$mean
        sd1 <- pred_one$sd * s
        mu2 <- pred_two$mean
        sd2 <- pred_two$sd * s
        pit <- w1 * pnorm(x, mu1, sd1) + w2 * pnorm(x, mu2, sd2)
        mu_comb <- w1 * mu1 + w2 * mu2
        var_comb <- w1 * (mu1^2 + sd1^2) + w2 * (mu2^2 + sd2^2) - mu_comb^2
        dss <- ((x - mu_comb)^2)/var_comb + 2 * log(sqrt(var_comb))
        predvariance <- var_comb
        crps <- crps_norm(x = x, mu1 = mu1, sd1 = sd1, mu2 = mu2,
            sd2 = sd2, w1 = w1)
    }
    pit_var <- var(pit, na.rm = TRUE)
    predvariance_rmv <- sqrt(mean(predvariance, na.rm = TRUE))
    out <- list(pit = pit, crps = crps, dss = dss, pitvar = pit_var,
        rmv = predvariance_rmv)
    class(out) <- "veri_stats"
    out
}
#' Synopsis of Verification Statistics
#' @export
#' @description Computes a synopsis for an object of
#' class 'veri_stats'.
#' @param x An object of class \code{'veri_stats'}.
#' @param method_name The name of the applied forecast method.
#' @param digits An integer specifying the digits for \code{round}. If
#' set to \code{NULL}, the values are not rounded.
#' @return A data frame with one row and corresponding row name given by \code{method_name}.
#' @examples
#' x <- veri_stats(17.5, pred_one = list(mean = 15,sd = 1.5))
#' veri_synop(x, method = "No method", digits = 3)
#' #' @author J. Gross, A. Moeller.
veri_synop <- function(x, method_name = NULL, digits = 4) {
    if (!(class(x) == "veri_stats")) stop("input must be of class 'veri_stats'")
    out <- data.frame(Length = length(x$crps),
                            Av.CRPS = mean(x$crps),
                            Av.DSS = mean(x$dss),
                            PIT.var = x$pitvar,
                            RMV = x$rmv)
    rownames(out) <- method_name
    if (!is.null(digits)) out <- round(x= out, digits = digits)
    out
}
