#'
#' Robust Partial and Semipartial Correlation
#'
#' @description Calculates a partial or semipartial correlation using one of the
#'    specified robust methods Taba linear or Taba rank correlation.
#' @usage taba.partial(x, y, ..., regress, method = c("taba", "tabarank"),
#'              alternative = c("less", "greater", "two.sided"),
#'              semi = c("none", "x", "y"), omega = 0.45)
#' @param x A numeric vector of length greater than 2 must be same length as y and covariates
#'          listed in ...
#' @param y A numeric vector of length greater than 2 must be same length as x and covariates
#'          listed in ...
#' @param ... Numeric vectors used as covariates of length equal to x and y
#' @param regress A string variable "\code{linear}" for linear regression, "\code{logistic}" for binary
#'   logistic regression, and "\code{poisson}" for Poisson regression
#' @param method A character string of \code{"taba"} or \code{"tabarank"} determining
#'   if one wants to calculate Taba linear or Taba rank (monotonic) correlation,
#'   respectively. If no method is specified, the function will output Taba
#'   linear correlation.
#' @param alternative Character string specifying the alternative hypothesis must be one
#'    of \code{"less"} for negative association, \code{"greater"} for
#'    positive association, or \code{"two.sided"} for difference in association.
#'    If the alternative is not specified, the function will default to a two sided test.
#' @param semi A character string specifying which variable (x or y) should be adjusted.
#' @param omega Numeric allowing the user to alter the tuning constant. If one is not specified,
#'   the function will default to 0.45. Range is between 0 and 1.
#' @details This function calculates the partial or semipartial association of two
#'    numeric vectors, or columns of a matrix or data frame composed
#'    of more than two numeric elements, adjusting for covariates of length equal to
#'    x and y. Covariates are combined colomn-wise and can be numeric vectors, matricies,
#'    or data frames with numeric cells. Each column in the matrix or data frame will be
#'    treated as a different covariate, and must have different names from x and y.
#'    Missing values in x, y, or any of the covariates are deleted row-wise.
#'    The default for this function is a two sided test using Taba linear partial
#'    correlation, with the tuning constant \code{omega} equal to 0.45.
#'    Thr variable you are not controlling must be continuous when using semipartial correlation.
#' @return This function returns the robust linear or monotonic association
#'   between two numeric vectors, adjusting for specified covariates. In addition,
#'   this function can provide the semipartial correlation, if specified.
#' @seealso
#'   \code{\link{taba}} for calculating Taba linear or Taba rank (monotonic) correlations
#'   \cr\code{\link{taba.test}} for testing Taba linear or Taba rank (monotonic) correlations
#'   \cr\code{\link{taba.gpartial}} for generalized partial correlations
#'   \cr\code{\link{taba.matrix}} for calculating correlation, p-value, and distance matricies
#' @references The paper is under review for possible publication.
#' @examples
#' x = rnorm(100)
#' y = rnorm(100)
#' z1 = rnorm(100)
#' z2 = rnorm(100)
#' z3 = rnorm(100)
#' taba.partial(x, y, z1, z2, z3, method = "tabarank")
#' taba.partial(x, y, z2, alternative = "less", semi = "x")
#' @import robustbase
#'         stats
#' @export taba.partial

taba.partial = function(x, y, ..., regress, method = c("taba", "tabarank"),
                        alternative = c("less", "greater", "two.sided"),
                        semi = c("none", "x", "y"), omega = 0.45) {
  if (missing(method)) {
    method <- "taba"
  }
  na.method <- pmatch(method, c("taba", "tabarank"))
  if (is.na(na.method)) {
    stop("invalid 'methods' argument")
    method <- match.arg(method)
  }
  if (missing(regress)) {
    regress <- "linear"
  }
  na.regress <- pmatch(regress, c("linear", "logistic", "poisson"))
  if (is.na(na.regress)) {
    stop("invalid 'regress' argument")
    regress <- match.arg(regress)
  }
  if (missing(alternative)) {
    alternative <- "two.sided"
  }
  na.alternative <- pmatch(alternative, c("less","greater","two.sided"))
  if (is.na(na.alternative)) {
    stop("invalid 'alternative' argument")
    alternative <- match.arg(alternative)
  }
  if (missing(semi)) {
    semi <- "none"
  }
  na.semi <- pmatch(semi, c("none", "x", "y"))
  if (is.na(na.semi)) {
    stop("invalid 'semi' argument")
    semi <- match.arg(semi)
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
  Covariates <- cbind.data.frame(...)
  covlen <- length(Covariates)
  if (covlen == 0) {
    stop("No covariates entered")
  }
  if (sum(sapply(Covariates,is.numeric)) != covlen) {
      stop("All covariates must be numeric")
      stopifnot(is.atomic(y))
  }
  if (sum(is.na(x)) > 0 || sum(is.na(y)) > 0 || sum(is.na(Covariates)) > 0) {
    warning("Missing data included in dataset was removed row-wise. Results may not be accurate.")
    miss <- which(complete.cases(x,y,Covariates) == FALSE)
    x <- x[-miss]
    y <- y[-miss]
    Covariates = as.data.frame(Covariates[-miss,])
  }
  k = length(x)
  if (k != length(y)) {
    stop("'x','y', and 'Covariares' must have the same length")
  }
  if (semi != "y") {
    xres <- switch(regress, "linear" = lm(x~., data = Covariates)$residuals,
                   "logistic" = glm(x~., family = binomial(link = "logit"),
                                    data = Covariates)$residuals,
                   "poisson" = glm(x~., family = poisson(link = "log"),
                                   data = Covariates)$residuals)
  }else{
    xres <- x
  }
  if (semi != "x") {
    yres  <- switch(regress, "linear" = lm(y~., data = Covariates)$residuals,
                    "logistic" = glm(y~., family = binomial(link = "logit"),
                                     data = Covariates)$residuals,
                    "poisson" = glm(y~., family = poisson(link = "log"),
                                    data = Covariates)$residuals)

  }else{
    yres <- y
  }
  if (method == "tabarank") {
    xres <- rank(xres)
    yres <- rank(yres)
  }
  if (Sn(xres) == 0 || Sn(yres) == 0) {
    s1 <- 1
    s2 <- 1
  } else {
    s1 <- Sn(xres)
    s2 <- Sn(yres)
  }
  medx <- median(xres)
  medy <- median(yres)
  a <- sum( ((1 / cosh(omega * ((xres - medx) / s1))) * ((xres - medx) / s1)) *
            ((1 / cosh(omega * ((yres - medy) / s2))) * ((yres - medy) / s2))    )
  b <- sum( ((1 / cosh(omega * ((xres - medx) / s1))) * ((xres - medx) / s1))**2 )
  c <- sum( ((1 / cosh(omega * ((yres - medy) / s2))) * ((yres - medy) / s2))**2 )
  tcor <- a / sqrt(b * c)
  tTaba <- ( tcor * sqrt((k - 2) / (1 - tcor**2)) )
  if (alternative == "two.sided") {
    p <- 2*pt(-abs(tTaba), (k - 2))
  }else{
    if (alternative == "greater") {
      p <- pt(-abs(tTaba), (k - 2), lower.tail = TRUE)
    }else{
      p <- pt(-abs(tTaba), (k - 2), lower.tail = FALSE)
    }
  }
  TabaC <- list(correlation = tcor,
                t.statistic = tTaba,
                p.value = p )
  return(TabaC)
}
