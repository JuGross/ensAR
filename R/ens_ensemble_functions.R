#' Rank of Observation
#' @export
#' @importFrom stats complete.cases
#' @description Computes the rank of an observation
#' among the corresponding ensemble forecasts.
#' @param ens A data frame with one observation column and at least
#' one forecast column.
#' @param obs_col The observation column.
#' @param mem_col The column(s) of the forecast members(s).
#' @return A vector of the same length as rows of \code{ens}.
#' @details Ranks are computed with \code{\link[base]{rank}} and
#' \code{ties.method = 'first'}.
#' @author J. Gross, A. Moeller.
#' @examples
#' x <- rank_ensemble(Magdeburg, obs_col = 6, mem_col = 7:56)
#' table(x)

rank_ensemble <- function(ens, obs_col, mem_col) {
    r_ens <- cbind(ens[, obs_col], ens[, mem_col])
    give_rank <- function(x) {
        if (any(!complete.cases(x))) {
            out <- NA
        } else {
            out <- rank(x, ties.method = "first")[1]
        }
        out
    }
    apply(r_ens, 1, give_rank)
}
#' AR Ensemble
#' @export
#' @importFrom stats ar
#' @importFrom stats predict
#' @importFrom zoo na.spline
#' @description For each member of a forecast ensemble,
#' a corresponding autoregressive modified member is generated.
#'
#' @param ens A data frame with one observation column, at least
#' one forecast column, and at least one additional column (e.g. date).
#' @param obs_col The observation column.
#' @param mem_col The column(s) of the forecast members(s).
#' @param train The length of the rolling training period used for
#' fitting an autoregressive process.
#' @param skip A number corresponding to the forecast ahead time (0 for ahead times
#'  not greater than 24 hours, 1 for ahead times greater than 24 hours and not greater
#'  than 48 hours, and so on).
#' @return A list with four elements:
#' \describe{
#' \item{\code{observation}}{A numeric vector containing the
#' observations for the forecast period (original period except for the
#' first \code{train} dates.}
#' \item{\code{forecast}}{A data frame containing the AR modified forecasts.}
#' \item{\code{variance}}{A data frame containing variance estimates
#' corresponding to the forecasts.}
#' \item{\code{additional}}{A data frame containg the variables (columns) of
#' \code{ens} which were
#' not specified by \code{obs_col} and \code{mem_col}.}
#' }
#'
#' @note It is assumed that in each row of the data frame \code{ens}
#' the forecast matches the observation,
#' i.e. for each row the difference between forecast and observation
#' is the forecast error.
#'
#' @author J. Gross, A. Moeller.
#' @examples
#' ar_ensemble(ens = Magdeburg[1:(90 + 1), -c(57,58)],
#'     obs_col = 6, mem_col = 7:56)
#' ar_ensemble(ens = Magdeburg48[1:(90 + 1), -c(57,58)],
#'     obs_col = 6, mem_col = 7:56, skip = 1)
#'
ar_ensemble <- function(ens, obs_col, mem_col, train = 90, skip = 0) {
    if (!is.data.frame(ens)) 
        stop("ensemble input must be a data frame")
    if (any(!complete.cases(ens[, c(obs_col, mem_col)]))) {
        warning("occurring NA's were replaced by spline interpolation")
        ens_approx <- apply(ens[, c(obs_col, mem_col)], 2, zoo::na.spline, 
            na.rm = FALSE)
        ens[, obs_col] <- ens_approx[, 1]
        ens[, mem_col] <- ens_approx[, -1]
    }
    train_length <- train
    forecast_period <- (train_length + 1):nrow(ens)
    y <- ens[, obs_col]  # observation for whole period
    #--
    ar_member <- function(member) {
        # to be applied to ensemble members
        ar_modify <- function(forecast_day) {
            # to be applied to forecast period
            train_period <- (forecast_day - (train_length:1))
            z <- (y - member)[train_period]  # forecast error series
            # --
            if (skip > 0) {
                zs <- z[1:(length(z) - skip)]
                zs_mod <- ar(x = zs, aic = TRUE, order.max = NULL)
                zs_pred <- predict(object = zs_mod, n.ahead = skip, se.fit = FALSE)
                z <- c(zs, zs_pred)
            }
            # --
            z_mod <- ar(x = z, aic = TRUE, order.max = NULL)
            p <- z_mod$order
            mu <- z_mod$x.mean
            a <- z_mod$ar
            # one-step ahead prediction of obs for time k = forecast_day based on
            # times k-1, ..., k-p:
            out_forecast_day <- member[forecast_day] + mu + sum(a * (z[(train_length - 
                (1:p))] - mu))
            # estimated variance of z from model fit:
            out_var_day <- var_ar(ar = a, i_var = z_mod$var.pred)
            return(c(forecast = out_forecast_day, variance = out_var_day))
        }
        out <- vapply(forecast_period, ar_modify, numeric(2))
        return(list(forecast = out[1, ], variance = out[2, ]))
    }
    ret_ens <- apply(ens[, mem_col, drop = FALSE], 2, ar_member)
    # -
    ens_forecast <- as.data.frame(lapply(ret_ens, function(x) x[[1]]))
    ens_var <- as.data.frame(lapply(ret_ens, function(x) x[[2]]))
    # 
    out <- list(observation = ens[forecast_period, obs_col], forecast = ens_forecast, 
        variance = ens_var, additional = ens[forecast_period, -c(obs_col, 
            mem_col)])
    class(out) <- "ar_ens"
    out
    
}
#' Predictive Moments from an AR Ensemble
#' @export
#' @importFrom stats pnorm dnorm optim
#' @description Computes the predictive mean and predictive standard
#' deviation based on a
#' model designed to handle an AR modified ensemble by \code{\link{ar_ensemble}}.
#' @details The predictive mean \code{mu} is the usual mean of
#' the AR modified ensemble members (data frame \code{forecast} from \code{ar_ens}).
#' The predictive standard deviation \code{sd} is the weighted mean of two
#' standard deviations. The first one is the square root of the mean
#' of AR variances (data frame \code{variance} from \code{ar_ens}). The second one
#' is the sample standard deviation (with denominator \eqn{n}, not \eqn{n-1})
#' of the AR modified ensemble members.
#' The weights are chosen in order to minimize the average CRPS computed from a rolling
#' training period of length \code{train} and assuming a predictive Gaussian distribution.
#' @param ar_ens A list generated by \code{\link{ar_ensemble}}.
#' @param train The length of the training period.
#' @return A data frame containing
#' \itemize{
#' \item{\code{obs}:}{ the observation}
#' \item{\code{mu}:}{ predictive mean}
#' \item{\code{sd}:}{ predictive standard deviation}
#' \item{\code{w}:}{ the weight corresponding to the first employed standard deviation}}
#' and additional columns as given by \code{additional} when invoking \code{\link{ar_ensemble}}.
#' @author J. Gross, A. Moeller.
#' @examples
#' mod <- ar_ensemble(ens = Magdeburg[1:(90 + 30 + 1), -c(57,58)],
#'     obs_col = 6, mem_col = 7:56)
#' ar_preddistr(mod) # data frame of one row
ar_preddistr <- function(ar_ens, train = 30) {
    if (!(class(ar_ens) == "ar_ens")) 
        stop("input must be of class 'ar_ens'")
    y <- ar_ens$observation
    n <- length(y)
    m <- train
    if (n <= m) 
        stop("too less observations")
    sample_sd <- function(x) sqrt(mean(scale(x, scale = FALSE)^2))
    # length n objects
    mu_out <- apply(ar_ens$forecast, 1, mean)
    sd_out_1 <- sqrt(apply(ar_ens$variance, 1, mean))
    sd_out_2 <- apply(ar_ens$forecast, 1, sample_sd)
    # length n - m objects
    v_period <- (m + 1):n
    w_out <- rep(1, (n - m))
    sd_out <- sd_out_1[v_period]
    # recomputation of w_out and sd_out in case m > 0
    if (m > 0) {
        roll_preddistr <- function(v_day) {
            train_period <- (v_day - (m:1))
            m_crps <- function(par) {
                m_obs <- y[train_period]
                m_mu <- mu_out[train_period]
                m_sd <- par[1] * sd_out_1[train_period] + (1 - par[1]) * 
                  sd_out_2[train_period]
                # vals <- crps_norm(m_obs, m_mu, m_sd)
                z <- (m_obs - m_mu)/m_sd
                vals <- m_sd * (z * (2 * pnorm(z) - 1) + 2 * dnorm(z) - 
                  1/sqrt(pi))
                mean(vals)
            }
            res <- optim(par = c(w = 0.5), fn = m_crps, method = "L-BFGS-B", 
                lower = 0, upper = 1)
            res$par["w"]
        }
        w_out <- vapply(v_period, roll_preddistr, numeric(1))
        # sd_out <- w_out * (sd_out_1[v_period]) + (1 - w_out) *
        # sd_out_2[v_period]
        sd_out <- .cx_comb(sd_out_1[v_period], sd_out_2[v_period], w_out)
    }
    out <- cbind(obs = y[v_period], mu = mu_out[v_period], sd = sd_out, 
        w = w_out, ar_ens$additional[v_period, ])
    out
}
#' AR Ensemble and Predictive Distribution
#' @export
#' @description A wrapper for sequential processing of the functions
#' \code{\link{ar_ensemble}} and \code{\link{ar_preddistr}}.
#' @param ens A data frame with one observation column, at least
#' one forecast column, and at least one additional column (e.g. date).
#' @param obs_col The observation column.
#' @param mem_col The column(s) of the forecast members(s).
#' @param skip A number corresponding to the forecast ahead time (0 for ahead times
#'  not greater than 24 hours, 1 for ahead times greater than 24 hours and not greater
#'  than 48 hours, and so on).
#' @param train_ar The length of the rolling training period used for
#' fitting an autoregressive process.
#' @param train_crps The length of the additional training period used
#' for computing the predictive standard deviation.
#' @return A data frame.
#' @author J. Gross, A. Moeller.
#' @examples
#' ensembleAR(ens = Magdeburg[1:(90 + 30 + 1), -c(57,58)],
#'     obs_col = 6, mem_col = 7:56)
#'
ensembleAR <- function(ens, obs_col, mem_col, skip = 0, train_ar = 90, 
    train_crps = 30) {
    out1 <- ar_ensemble(ens = ens, obs_col = obs_col, mem_col = mem_col, 
        train = train_ar, skip = skip)
    out2 <- ar_preddistr(out1, train = train_crps)
    out2
}


