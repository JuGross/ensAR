#' Finding Gaps of NAs
#' @export
#' @description Computes lengths and frequencies of NA gaps, i.e.
#' sequences of consecutive NAs, in a vector or data frame.
#'
#' @details If the input is a data frame, a missing value is always defined as a
#' not complete case, i.e. each row containing an NA is treated as a
#' missing observation.
#'
#' @param x A vector or data frame.
#' @return A list with three elements of class 'na_gaps'.
#' \describe{
#' \item{\code{table}}{A table of gap lengths and corresponding
#' frequencies.}
#' \item{\code{lengths}}{A vector of gap lengths, actually the names
#' of list element \code{table}.}
#' \item{\code{number}}{The number of incomplete cases.}
#' }
#' @author J. Gross, A. Moeller.
#' @seealso \code{\link[stats]{complete.cases}}
#' @examples
#' na_gaps(c(1, NA, NA, 1, 1, NA, NA, NA, 1, NA, NA, 1, 1, 1, NA, NA))
na_gaps <- function(x) {
    na_ind <- as.numeric(!complete.cases(x))
    na_numb <- sum(na_ind)
    if (na_numb == length(na_ind)) {
        gap_vec <- c(1)
        names(gap_vec) <- as.character(na_numb)
        gap_table <- as.table(gap_vec)
        gap_lengths <- na_numb
    }
    if (na_numb < length(na_ind)) {
        not_na <- 1 - na_ind
        y <- cumsum(not_na) + 1
        y_na <- y * na_ind
        tab1 <- table(y_na)[1]
        tab2 <- table((as.vector(table(y_na)[-1])))
        gap_table <- as.table(c(tab1, tab2))
        gap_lengths <- as.numeric(names(gap_table))
    }
    out <- list(table = gap_table, lengths = gap_lengths, number = na_numb)
    class(out) <- "na_gaps"
    out
}
#--
#' Print NA Gaps
#' @export
#' @description Print method for objects of class 'na_gaps'.
#' @param x An object of class 'na_gaps'
#' @param ... Further arguments passed to S3 method for class 'table'
#' @author J. Gross, A. Moeller.
#' @seealso \code{\link{na_gaps}}
print.na_gaps <- function(x, ...) print.table(x$table, ...)
# --
#' Conversion of Character Vector of Dates
#' @export
#' @description Converts a character vector into a Date object
#' and checks for missing dates.
#' @param x A character vector with dates in format 'yyyymmdd' with no NA's and no duplicates.
#' @return A list with two elements:
#'  \describe{
#' \item{\code{dates}}{Date object in format yyyy-mm-dd, chronologically ordered.}
#' \item{\code{missing}}{Non-occuring dates between chronolocically first and last input date.}
#' }
#' @author J. Gross, A. Moeller.
#' @note If \code{missing} is 'None' and no warning is given, the
#' dates in \code{x} are complete and in chronological order.
#' @note The function is not thoroughly protected against an incorrect input format.
#' @examples x <- c('20020105', '20020102', '20020103')
#' convert_chardate(x)
#'
# --
convert_chardate <- function(x) {
    if (any(is.na(x))) 
        stop("input dates contain NA's")
    if (any(duplicated(x))) 
        stop("input dates contain duplicates")
    x_dates <- as.Date(x, format = "%Y%m%d")
    if (is.unsorted(x_dates)) {
        x_dates <- sort(x_dates)
        warning("input dates are not in chronological order")
    }
    obs_dates <- chron::chron(dates. = as.character(x_dates), 
        format = c(dates = "y-m-d"), out.format = c(dates = "y-m-d"))
    n <- length(obs_dates)
    gen_dates <- chron::seq.dates(obs_dates[1], obs_dates[n], 
        by = "days")
    if (n < length(gen_dates)) {
        mis_dates <- gen_dates[(!(gen_dates %in% obs_dates))]
        mis_dates <- as.Date(mis_dates, format = "%Y%m%d")
    } else {
        mis_dates <- "None"
    }
    out <- list(dates = x_dates, missing = mis_dates)
    out
}
#' Relative Frequency of Previous Equal Elements
#' @export
#' @description Compares each feasible element of a vector to a given number of
#' previous elements and computes the relative frequency of successes,
#' where a success occurs in case of positive equality comparison.
#'
#' @details For each element of \code{x} greater than \code{previous},
#' a success is recorded when it is equal to the last \code{previous} elements,
#' where equality is assessed by \code{\link{all.equal}}.
#'
#' @param x A numeric vector.
#' @param previous The number of previous elements compared to the actual element.
#' @return Relative frequency of possible successes as decribed in the Details section.
#'
#' @author J. Gross, A. Moeller.
#' @examples
#' x <- c(1, 1, 2, 1, 2, 2, 2, 3)
#' previous_equal(x, previous = 1) # 3/7
#' previous_equal(x, previous = 2) # 1/6
#'
previous_equal <- function(x, previous = 10) {
    lower <- previous + 1
    upper <- length(x)
    if (lower >= upper) 
        stop("Not possible")
    range_equal <- function(element) {
        isTRUE(all.equal(x[(element + 1 - lower):element], rep(x[element], 
            lower)))
    }
    out_equal <- vapply(lower:upper, range_equal, logical(1))
    # out_equal
    sum(out_equal)/length(lower:upper)
}
# --

