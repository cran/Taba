#'
#' Robust Correlation Test
#'
#' @description Tests the association between two numeric vectors using Taba
#'     robust linear or Taba rank (monotonic) correlation.
#' @usage taba.test(x, y, method = c("taba", "tabarank"),
#'           alternative = c("less", "greater", "two.sided"),
#'           omega = 0.45)
#' @param x A numeric vector of length greater than 2 must be same length as y
#' @param y A numeric vector of length greater than 2 must be same length as x
#' @param method A character string of \code{"taba"} or \code{"tabarank"}
#'   determining if one wants to calculate Taba linear or Taba rank (monotonic) correlation,
#'   respectively. If no method is specified, the function will output Taba
#'   linear correlation.
#' @param alternative Character string specifying the alternative hypothesis must be one
#'    of \code{"less"} for negative association, \code{"greater"} for
#'    positive association, or \code{"two.sided"} for difference in association.
#'    If the alternative is not specified, the function will default to a two sided test.
#' @param omega Numeric allowing the user to alter the tuning constant. If one is not specified,
#'   the function will default to 0.45. Range is between 0 and 1.
#' @details This function tests the association of two non-empty numeric vectors of
#'    length greater than two, or two columns of a data frame or matrix composed
#'    of more than two numeric elements. Covariates are combined colomn-wise and can be
#'    numeric vectors, matricies, or data frames with numeric cells. Each column in the
#'    matrix or data frame will be treated as a different covariate, and must have
#'    different names. Missing values in either x or y are deleted row-wise. The two sided
#'    test with the null hypothesis correlation is equal to zero. The default is a two
#'    sided test using Taba Linear correlation, with the tuning constant \code{omega}
#'    equal to 0.45.
#' @return This function returns the robust linear or monotonic association
#'   between two numeric vectors, along with it's respective test statistic, and p-value.
#' @seealso
#'   \code{\link{taba}} for calculating Taba linear or Taba rank (monotonic) correlations
#'   \cr\code{\link{taba.partial}} for partial and semipartial correlations
#'   \cr\code{\link{taba.gpartial}} for generalized partial correlations
#'   \cr\code{\link{taba.matrix}} for calculating correlation, p-value, and distance matricies
#' @references The paper is under review for possible publication.
#' @examples
#' x = rnorm(10)
#' y = rnorm(10)
#' taba.test(x,y)
#' taba.test(x,y,method = "tabarank", alternative = "less")$p.value
#' @import robustbase
#'         stats
#' @export

taba.test = function(x, y, method = c("taba", "tabarank"),
                     alternative = c("less", "greater", "two.sided"),
                     omega = 0.45) {
  if (missing(method)) {
    method <- "taba"
  }
  na.method <- pmatch(method, c("taba", "tabarank"))
  if (is.na(na.method)) {
    stop("invalid 'method' argument")
    method <- match.arg(method)
  }
  if (missing(alternative)) {
    alternative <- "two.sided"
  }
  na.alternative <- pmatch(alternative, c("less","greater","two.sided"))
  if (is.na(na.alternative)) {
    stop("invalid 'alternative' argument")
    alternative <- match.arg(alternative)
  }
  if (missing(omega)) {
    omega <- 0.45
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
    miss <- which(complete.cases(x,y)==FALSE)
    x <- x[-miss]
    y <- y[-miss]
  }
  k <- length(x)
  if (k != length(y)) {
    stop("'x' and 'y' must have the same length")
  }
  if (method == "tabarank") {
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
  medx <- median(x)
  medy <- median(y)
  a <- sum( ((1 / cosh(omega * ((x - medx) / s1))) * ((x - medx) / s1)) *
            ((1 / cosh(omega * ((y - medy) / s2))) * ((y - medy) / s2))    )
  b <- sum( ((1 / cosh(omega * ((x - medx) / s1))) * ((x - medx) / s1))**2 )
  c <- sum( ((1 / cosh(omega * ((y - medy) / s2))) * ((y - medy) / s2))**2 )
  tcor <- a / sqrt(b * c)
  t    <- tcor * sqrt( (k - 2) / (1 - tcor**2) )
  if (alternative == "two.sided") {
    p <- 2*pt(-abs(t), (k - 2))
  }else{
    if (alternative == "greater") {
      p <- pt(-abs(t), (k - 2), lower.tail = TRUE)
    }else{
      p <- pt(-abs(t), (k - 2), lower.tail = FALSE)
    }
  }
  out  <- list(correlation = tcor,
               t.statistic = t,
               p.value     = p )
  return(out)
}