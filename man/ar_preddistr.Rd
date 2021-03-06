% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ens_ensemble_functions.R
\name{ar_preddistr}
\alias{ar_preddistr}
\title{Predictive Moments from an AR Ensemble}
\usage{
ar_preddistr(ar_ens, train = 30)
}
\arguments{
\item{ar_ens}{A list generated by \code{\link{ar_ensemble}}.}

\item{train}{The length of the training period.}
}
\value{
A data frame containing
\itemize{
\item{\code{obs}:}{ the observation}
\item{\code{mu}:}{ predictive mean}
\item{\code{sd}:}{ predictive standard deviation}
\item{\code{w}:}{ the weight corresponding to the first employed standard deviation}}
and additional columns as given by \code{additional} when invoking \code{\link{ar_ensemble}}.
}
\description{
Computes the predictive mean and predictive standard
deviation based on a
model designed to handle an AR modified ensemble by \code{\link{ar_ensemble}}.
}
\details{
The predictive mean \code{mu} is the usual mean of
the AR modified ensemble members (data frame \code{forecast} from \code{ar_ens}).
The predictive standard deviation \code{sd} is the weighted mean of two
standard deviations. The first one is the square root of the mean
of AR variances (data frame \code{variance} from \code{ar_ens}). The second one
is the sample standard deviation (with denominator \eqn{n}, not \eqn{n-1})
of the AR modified ensemble members.
The weights are chosen in order to minimize the average CRPS computed from a rolling
training period of length \code{train} and assuming a predictive Gaussian distribution.
}
\examples{
mod <- ar_ensemble(ens = Magdeburg[1:(90 + 30 + 1), -c(57,58)],
    obs_col = 6, mem_col = 7:56)
ar_preddistr(mod) # data frame of one row
}
\author{
J. Gross, A. Moeller.
}
