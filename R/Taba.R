#
#' Robust Correlation
#'
#' @description Returns the Taba robust linear, Taba rank (monotonic), TabWil, or
#'    TabWil rank correlation coefficient between two numeric vectors.
#' @usage taba(x, y, method = c("taba", "tabarank","tabwil", "tabwilrank"), omega)
#' @param x A numeric vector of length greater than 2 must be same length as y
#' @param y A numeric vector of length greater than 2 must be same length as x
#' @param method A character string of \code{"taba"}, \code{"tabarank"}, \code{"tabwil"}, or
#'    \code{"tabwilrank"} determining if one wants to calculate Taba linear, Taba rank
#'    (monotonic), TabWil, or TabWil rank correlation, respectively. If no method is specified,
#'    the function will output Taba Linear correlation.
#' @param omega Numeric allowing the user to alter the tuning constant. If one is not specified,
#'   the function will default to 0.45 for Taba and Taba rank, and 0.1 for TabWil and TabWil rank.
#'   Range is between 0 and 1.
#' @details This function can be used to compare two non-empty numeric vectors of
#'    length greater than two, or two columns of a data frame or matrix composed
#'    of more than two numeric elements. Missing values in either x or y are
#'    deleted row-wise. The default method is Taba Linear correlation, with the
#'    tuning constant \code{omega}.
#' @return This function returns a the robust linear or monotonic association
#'   between two numeric vectors as a numeric.
#' @seealso
#'   \code{\link{taba.test}} for testing Taba linear or Taba rank (monotonic) correlations
#'   \cr\code{\link{taba.partial}} for partial and semipartial correlations
#'   \cr\code{\link{taba.gpartial}} for generalized partial correlations
#'   \cr\code{\link{taba.matrix}} for calculating correlation, p-value, and distance matricies
#' @references Tabatabai, M., Bailey, S., Bursac, Z. et al. An introduction to new robust linear
#'   and monotonic correlation coefficients. BMC Bioinformatics 22, 170 (2021). https://doi.org/10.1186/s12859-021-04098-4
#'   \cr{\cr{\doi{https://doi.org/10.1186/s12859-021-04098-4}}}
#' @examples
#' x = rnorm(100)
#' y = rnorm(100)
#' taba(x, y)
#' taba(x, y, method = "tabarank", omega = 0.4)
#' taba(x, y, method = "tabwil", omega = 0.22)
#' @import robustbase
#'         stats
#' @export taba

taba = function(x, y, method = c("taba", "tabarank","tabwil", "tabwilrank"), omega) {
  if (missing(method)) {
    method <- "taba"
  }
  na.method <- pmatch(method, c("taba", "tabarank","tabwil", "tabwilrank"))
  if (is.na(na.method)) {
    stop("invalid 'method' argument")
    method <- match.arg(method)
  }
  if (missing(omega)) {
    if (method == "taba" || method == "tabarank") {
      omega <- 0.45
    } else {
      omega <- 0.05
    }
  }
  if (omega > 1 || omega < 0) {
    stop("'omega' must be between 0 and 1")
    omega <- match.arg(omega)
  }
  if (is.data.frame(y) || is.numeric(y)) {
    y <- as.matrix(y)
  }
  if (is.data.frame(x) || is.numeric(x)) {
    x <- as.matrix(x)
  }
  if (!is.matrix(x) && is.null(y)) {
    stop("supply both 'x' and 'y' or a matrix-like 'x'")
  }
  if (!(is.numeric(x) || is.logical(x))) {
    stop("'x' must be numeric")
    stopifnot(is.atomic(x))
  }
  if (!is.null(y)) {
    if (!(is.numeric(y) || is.logical(y)))
      stop("'y' must be numeric")
      stopifnot(is.atomic(y))
  }
  if (sum(is.na(x)) > 0 || sum(is.na(y)) > 0) {
    warning("Missing data included in dataset was removed row-wise. Results may not be accurate.")
    miss = which(complete.cases(x,y) == FALSE)
    x <- x[-miss]
    y <- y[-miss]
  }
  if (length(x) != length(y)) {
    stop("'x' and 'y' must have the same length")
  }
  if (method == "tabarank" || method == "tabwilrank") {
    x <- rank(x)
    y <- rank(y)
  }
  if (Sn(x) == 0 || Sn(y) == 0) {
    s1 <- 1
    s2 <- 1
  } else {
    s1 <- Sn(x)
    s2 <- Sn(y)
  }
  if (method == "taba" || method == "tabarank") {
    medx <- median(x)
    medy <- median(y)
    a <- sum( ((1 / cosh(omega * ((x - medx) / s1))) * ((x - medx) / s1)) *
              ((1 / cosh(omega * ((y - medy) / s2))) * ((y - medy) / s2))    )
    b <- sum( ((1 / cosh(omega * ((x - medx) / s1))) * ((x - medx) / s1))**2 )
    c <- sum( ((1 / cosh(omega * ((y - medy) / s2))) * ((y - medy) / s2))**2 )
    tcor <- a / sqrt(b * c)
  } else {
    u <- (x - median(x))/s1 + (y - median(y))/s2
    v <- (x - median(x))/s1 - (y - median(y))/s2
    a <-  ((1 / cosh(omega * (median(abs(u))**2))) * (median(abs(u))**2)) - ((1 / cosh(omega * (median(abs(v))**2))) * (median(abs(v))**2))
    b <-  ((1 / cosh(omega * (median(abs(u))**2))) * (median(abs(u))**2)) + ((1 / cosh(omega * (median(abs(v))**2))) * (median(abs(v))**2))
    tcor <- a / b
  }
  return(tcor)
}
