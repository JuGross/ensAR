#' Variance of an AR process
#' @export
#' @importFrom stats ARMAacf
#' @description Comptutes the variance of an autoregressive (AR) process.
#'
#' @param ar A vector of autoregression coefficients.
#' @param i_var The innovations variance.
#' @return The variance of the process.
#' @details The variance of an AR process of order \eqn{p} is
#' \deqn{\gamma_{0} = \frac{\sigma^2}{1 - \phi_{1} \rho_{1} - \phi_{2} \rho_{2} - \cdots
#'  - \phi_{p} \rho_{p}}}
#' where \eqn{\sigma^2} is the variance of the innovations and \eqn{\phi_{j}}
#' are the autoregression coefficients. The autocorrelations
#' \eqn{\rho_{j}} are computed by \code{\link{ARMAacf}}.
#'
#' @author J. Gross, A. Moeller.
#' @examples
#' ## x <- arima.sim(list(ar = 0.7), sd = 0.5, 100000)
#' ## var(x)
#' var_ar(ar = 0.7, i_var= 0.5^2)
var_ar <- function(ar = numeric(), i_var = 1) {
    p <- length(ar)
    if (p == 0)
        var_out <- i_var
    if (p == 1)
        var_out <- i_var/(1 - ar^2)
    if (p > 1)
        var_out <- i_var/sum(ARMAacf(ar = ar) * c(1, -ar))
    var_out
}
