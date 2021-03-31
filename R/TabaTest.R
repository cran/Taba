#'
#' Robust Correlation Test
#'
#' @description Tests the association between two numeric vectors using Taba robust linear,
#'    Taba rank (monotonic), TabWil, or TabWil rank correlation coefficient.
#' @usage taba.test(x, y, method = c("taba", "tabarank", "tabwil", "tabwilrank"),
#'           alternative = c("less", "greater", "two.sided"),
#'           omega, alpha = 0.05)
#' @param x A numeric vector of length greater than 2 must be same length as y
#' @param y A numeric vector of length greater than 2 must be same length as x
#' @param method A character string of \code{"taba"}, \code{"tabarank"}, \code{"tabwil"}, or
#'    \code{"tabwilrank"} determining if one wants to calculate Taba linear, Taba rank
#'    (monotonic), TabWil, or TabWil rank correlation, respectively. If no method is specified,
#'    the function will output Taba Linear correlation.
#' @param alternative Character string specifying the alternative hypothesis must be one
#'    of \code{"less"} for negative association, \code{"greater"} for
#'    positive association, or \code{"two.sided"} for difference in association.
#'    If the alternative is not specified, the function will default to a two sided test.
#' @param omega Numeric allowing the user to alter the tuning constant. If one is not specified,
#'   the function will default to 0.45 for Taba and Taba rank, and 0.1 for TabWil and TabWil rank.
#'   Range is between 0 and 1.
#' @param alpha Type I error rate. Numeric must be between 0 and 1. Default set to 0.05.
#' @details This function tests the association of two non-empty numeric vectors of
#'    length greater than two, or two columns of a data frame or matrix composed
#'    of more than two numeric elements. Covariates are combined colomn-wise and can be
#'    numeric vectors, matricies, or data frames with numeric cells. Each column in the
#'    matrix or data frame will be treated as a different covariate, and must have
#'    different names. Missing values in either x or y are deleted row-wise. The two sided
#'    test with the null hypothesis correlation is equal to zero. The default is a two
#'    sided test using Taba Linear correlation, with tuning constant \code{omega}.
#' @return This function returns the robust linear or monotonic association
#'   between two numeric vectors, along with it's respective test statistic, and p-value.
#' @seealso
#'   \code{\link{taba}} for calculating Taba linear or Taba rank (monotonic) correlations
#'   \cr\code{\link{taba.partial}} for partial and semipartial correlations
#'   \cr\code{\link{taba.gpartial}} for generalized partial correlations
#'   \cr\code{\link{taba.matrix}} for calculating correlation, p-value, and distance matricies
#' @references Tabatabai, M., Bailey, S., Bursac, Z. et al. An introduction to new robust linear
#'   and monotonic correlation coefficients. BMC Bioinformatics 22, 170 (2021). https://doi.org/10.1186/s12859-021-04098-4
#'   \cr{\cr{\doi{https://doi.org/10.1186/s12859-021-04098-4}}}
#' @examples
#' x = rnorm(10)
#' y = rnorm(10)
#' taba.test(x, y)
#' taba.test(x, y, method = "tabarank", alternative = "less")$p.value
#' taba.test(x, y, method = "tabwil", omega = .1)
#' @import robustbase
#'         stats
#' @export taba.test

taba.test = function(x, y, method = c("taba", "tabarank", "tabwil", "tabwilrank"),
                      alternative = c("less", "greater", "two.sided"),
                      omega, alpha = 0.05) {
  if (missing(method)) {
    method <- "taba"
  }
  na.method <- pmatch(method, c("taba", "tabarank","tabwil", "tabwilrank"))
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
    if (method == "taba" || method == "tabarank") {
      omega <- 0.45
    } else {
      omega <- 0.05
    }
  }
  if (missing(alpha)) {
    alpha <- 0.05
  }
  if (alpha > 1 || alpha < 0) {
    stop("'alpha' must be between 0 and 1")
    alpha <- match.arg(alpha)
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
  t    <- tcor * sqrt( (k - 2) / (1 - tcor**2) )
  if (alternative == "two.sided") {
    p <- 2*pt(-abs(t), (k - 2))
    small  <- tanh(atanh(tcor) - qnorm(1 - alpha/2) * sqrt((1 + 0.5*tcor**2) / (k - 3)))
    large  <- tanh(atanh(tcor) + qnorm(1 - alpha/2) * sqrt((1 + 0.5*tcor**2) / (k - 3)))
    lcl <- min(small,large)
    ucl <- max(small,large)
  }else{
    if (alternative == "greater") {
      p <- pt(-abs(t), (k - 2), lower.tail = TRUE)
      small  <- tanh(atanh(tcor) - qnorm(1 - alpha) * sqrt((1 + 0.5*tcor**2) / (k - 3)))
      large  <- tanh(atanh(tcor) + qnorm(1 - alpha) * sqrt((1 + 0.5*tcor**2) / (k - 3)))
    }else{
      p <- pt(-abs(t), (k - 2), lower.tail = FALSE)
      small  <- tanh(atanh(tcor) - qnorm(alpha) * sqrt((1 + 0.5*tcor**2) / (k - 3)))
      large  <- tanh(atanh(tcor) + qnorm(alpha) * sqrt((1 + 0.5*tcor**2) / (k - 3)))
    }
    lcl <- min(small,large)
    ucl <- max(small,large)
  }
  if (tcor == 1) {
    lcl <- 1
    ucl <- 1
    t <- Inf
    p <- 0
  }else if (tcor == -1){
    lcl <- -1
    ucl <- -1
    t <- Inf
    p <- 0
  }
  out  <- list(correlation = tcor,
               LCL = lcl,
               UCL = ucl,
               t.statistic = t,
               p.value     = p )
  return(out)
}
