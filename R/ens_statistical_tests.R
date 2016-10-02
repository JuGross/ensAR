#' Diebold-Mariano Test
#'
#' @export
#' @importFrom stats complete.cases acf pnorm
#' @description Performs the Diebold-Mariano test for
#' equal predictive performance of two methods with respect to
#' a scoring rule.
#'
#' @param s1 Score from method 1.
#' @param s2 Score from method 2.
#' @param alternative A character string specifying the alternative
#' hypothesis, must be one of 'two.sided' (default), 'greater' or
#' 'less'.
#' @param h An integer specifying the smallest lag not
#' used in the autocovariance function.
#' @return A object of class \code{'htest'}.
#' @details The null hypothesis is that the difference \code{s1 - s2} has zero mean.
#' The alternative \code{'less'} is that \code{s1 - s2} has negative mean.
#'
#' The test statistic requires values of the autocovariance function
#' of  \code{s1 - s2} with lags smaller than \code{h}, i.e.
#' \code{lag.max = h - 1} (the truncation lag)
#' in \code{\link[stats]{acf}}.
#'
#'  The difference \code{s1 - s2}  may contain missing values,
#'  in which case complete cases are used and a warning is given.
#'
#' @author J. Gross, A. Moeller.
#' @seealso \code{\link[forecast]{dm.test}} in package \pkg{forecast}.
#' @note The function \code{dm_test} is inspired by the function \code{dm.test} from
#' package \pkg{forecast}. It is a 'simpler' version in the sense that it computes
#' p-values from the normal instead of the t distribution.
#' @references
#' Diebold F.X, Mariano R.S. 1995. Comparing predictive accuracy.
#' \emph{Journal of Business & Economic Statistics}, \strong{13}, 253--263.
#'
#' Gneiting T., Katzfuss  M. 2014. Probabilistic forecasting.
#' \emph{Annual Review of Statistics and Its Application},
#' \strong{1}, 125--151.
#'
#' @examples
#' set.seed(1)
#' s1 <- arima.sim(list(ar = 0.7), sd = 0.5, 100)
#' s2 <- arima.sim(list(ar = 0.7), sd = 0.5, 100) - 0.1
#' dm_test(s1, s2)
dm_test <- function(s1, s2,
                    alternative = c("two.sided", "less", "greater"),
                    h = 1) {
    if (length(s1) != length(s2)) {
        stop("imput vectors must have same length")
    }
    alternative <- match.arg(alternative)
    dname <- paste(deparse(substitute(s1)), deparse(substitute(s2)))
    d <- s1 - s2
    if (any(is.na(d))) {
        warning("missig values: autocovariance estimate may not be valid")
    }
    d <- d[complete.cases(d)]
    n_d <- length(d)
    acf_est <- acf(d, type = "covariance", lag.max = h - 1, plot = FALSE)
    d_acf <- acf_est$acf[, , 1]
    d_var <- sum(c(d_acf[1], 2 * d_acf[-1]))/n_d
    if (d_var < 0) {
        S <- NA
        pval <- 0
    } else {
        S <- mean(d)/sqrt(d_var)
        if (alternative == "two.sided")
            pval <- 2 * pnorm(-abs(S)) else if (alternative == "less")
            pval <- pnorm(S) else if (alternative == "greater")
            pval <- pnorm(S, lower.tail = FALSE)
    }
    para <- h - 1
    names(para) = c("truncation lag")
    RVAL <- list(statistic = c(DW = S), parameter = para, p.value = pval,
        alternative = alternative, method = "Diebold-Mariano test", data.name = dname)
    class(RVAL) <- "htest"
    return(RVAL)
}
