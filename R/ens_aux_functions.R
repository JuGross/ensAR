#' Check for Weight
#' @description Checks whether all elements of vector lie between 0 and 1.
#' @param x A vector.
#' @return A logical vector of length one.
#'
#' @author J. Gross, A. Moeller.
#' @examples
#' \dontrun{
#' a <- c(0.3, 0.4, 0.9)
#' .is_weight(a)
#' }
#'
.is_weight <- function(x) {
    out <- !any((x < 0) | (x > 1))
    out
}
#'
#' Convex Combination
#' @description Computes a convex combination of two vectors of the same length
#' @param x A vector.
#' @param y A vector of the same length as x.
#' @param w A weight vector of the same length as x.
#' @return A vector of the same length as x and y.
#'
#' @author J. Gross, A. Moeller.
#' @examples \dontrun{
#' x <- c(1, 2, 3)
#' y <- c(3, 2, 1)
#' w <- c(1/2, 1/3, 1/2)
#' .cx_comb(x, y, w)
#' }
#'
.cx_comb <- function(x, y, w) {
    if (!((length(x) == length(y)) & (length(y) == length(w)))) 
        stop("input vectors must have same length")
    if (!.is_weight(w)) 
        stop("weights are not between 0 and 1")
    w * x + (1 - w) * y
}

