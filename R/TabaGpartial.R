#'
#' Generalized Taba Partial and Taba Rank Partial Correlation
#'
#' @description Calculates a generalized partial correlation using one of the
#'    specified robust methods Taba linear or Taba rank correlation.
#' @usage taba.gpartial(x, y, xcov, ycov, regress.x, regress.y,
#'               method = c("taba","tabarank","tabwil","tabwilrank"),
#'               alternative = c("less", "greater", "two.sided"),
#'               omega)
#' @param x A numeric vector of length greater than 2 must be same length as y and covariates
#'          listed in x and ycov
#' @param y A numeric vector of length greater than 2 must be same length as x and covariates
#'          listed in y and xcov
#' @param xcov A data frame, matrix, or numeric vectors combined columnwize used as covariates for x,
#'     which have length equal to x
#' @param ycov A data frame, matrix, or numeric vectors combined columnwize used as covariates for y,
#'     which have length equal to y
#' @param regress.x A string variable "\code{linear}" for linear regression, "\code{logistic}" for binary
#'   logistic regression, and "\code{poisson}" for Poisson regression
#' @param regress.y A string variable "\code{linear}" for linear regression, "\code{logistic}" for binary
#'   logistic regression, and "\code{poisson}" for Poisson regression
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
#' @details This function generalizes the partial correlation. In the event that the controlling
#'    variables for x and y are identical, it reduces to Taba, Taba rank, TabWil, and TabWil
#'    rank partial correlation. Covariates used to control for x
#'    should be represented columnwise in a matrix or data frame as \code{xcov}. Similarly,
#'    covariates used to control for y should be represented columnwise in a matrix or
#'    data frame as \code{ycov}. When controling an outcome variable with one covariate,
#'    a vector will suffice. Because x and y refer to the outcome varibales, names of
#'    covariates (or control variables) must not be named "x" or "y". The user has the
#'    option of using different regression methods when controling each outcome variable.
#'    Missing values in x, y, or any of the covariates are deleted row-wise. All categorical
#'    variables must be converted to type factor prior to using this function. \cr
#'    The default for this function is a two sided test using generalized partial Taba
#'    correlation using a linear regression to obtain residuals, with the tuning
#'    constant \code{omega} equal to 0.45.
#' @return This function returns the robust association
#'   between two numeric vectors, adjusting for specified covariates. In addition,
#'   this function can provide the semipartial correlation, if specified.
#' @seealso
#'   \code{\link{taba}} for calculating Taba linear or Taba rank (monotonic) correlations
#'   \cr\code{\link{taba.test}} for testing Taba linear or Taba rank Monotonic correlations
#'   \cr\code{\link{taba.partial}} for partial and semipartial correlations
#'   \cr\code{\link{taba.matrix}} for calculating correlation, p-value, and distance matricies
#' @references Tabatabai, M., Bailey, S., Bursac, Z. et al. An introduction to new robust linear
#'   and monotonic correlation coefficients. BMC Bioinformatics 22, 170 (2021). https://doi.org/10.1186/s12859-021-04098-4
#'   \cr{\cr{\doi{https://doi.org/10.1186/s12859-021-04098-4}}}
#' @examples
#' x = rnorm(100)
#' y = rnorm(100)
#' z1 = rnorm(100)
#' z2 = rnorm(100)
#' z3 = rnorm(100)
#' w = sample(c(0,1), replace=TRUE, size=100)
#' taba.gpartial(x, y, xcov = cbind(z1, z2), ycov = cbind(z1, z3), method = "tabarank")
#' taba.gpartial(x, y, z2, ycov = cbind(z1, z2), alternative = "less")
#' taba.gpartial(w, y, z1, cbind(z2, z3),regress.x = "logistic")
#' @import robustbase
#'         stats
#' @export taba.gpartial

taba.gpartial = function(x, y, xcov, ycov, regress.x, regress.y,
                         method = c("taba", "tabarank", "tabwil", "tabwilrank"),
                         alternative = c("less", "greater", "two.sided"),
                         omega) {
  if (missing(method)) {
    method <- "taba"
  }
  na.method <- pmatch(method, c("taba", "tabarank", "tabwil", "tabwilrank"))
  if (is.na(na.method)) {
    stop("invalid 'methods' argument")
    method <- match.arg(method)
  }
  if (missing(regress.x)) {
    regress.x <- "linear"
  }
  na.regress.x <- pmatch(regress.x, c("linear", "logistic", "poisson"))
  if (is.na(na.regress.x)) {
    stop("invalid 'regress.x' argument")
    regress.x <- match.arg(regress.x)
  }
  if (missing(regress.y)) {
    regress.y <- "linear"
  }
  na.regress.y <- pmatch(regress.y, c("linear", "logistic", "poisson"))
  if (is.na(na.regress.y)) {
    stop("invalid 'regress.y' argument")
    regress.y <- match.arg(regress.y)
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
  xcov <- cbind.data.frame(xcov)
  ycov <- cbind.data.frame(ycov)
  Covariates <- merge(xcov, ycov, by = intersect(colnames(xcov), colnames(ycov)))
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
    Covariates <- as.data.frame(Covariates[-miss,])
    xcov <- xcov[-miss,]
    ycov <- ycov[-miss,]
  }
  k = length(x)
  if (k != length(y)) {
    stop("'x','y', and 'Covariares' must have the same length")
  }
  xres <- switch(regress.x, "linear" = lm(x~., data = xcov)$residuals,
                 "logistic" = glm(x~., family = binomial(link = "logit"), data = xcov)$residuals,
                 "poisson" = glm(x~., family = poisson(link = "log"), data = xcov)$residuals)
  yres <- switch(regress.y, "linear" = lm(y~., data = ycov)$residuals,
                 "logistic" = glm(y~., family = binomial(link = "logit"), data = ycov)$residuals,
                 "poisson" = glm(y~., family = poisson(link = "log"), data = ycov)$residuals)
  if (method == "tabarank" || method == "tabwilrank") {
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
  if (method == "taba" || method == "tabarank") {
    medx <- median(xres)
    medy <- median(yres)
    a <- sum( ((1 / cosh(omega * ((xres - medx) / s1))) * ((xres - medx) / s1)) *
              ((1 / cosh(omega * ((yres - medy) / s2))) * ((yres - medy) / s2))    )
    b <- sum( ((1 / cosh(omega * ((xres - medx) / s1))) * ((xres - medx) / s1))**2 )
    c <- sum( ((1 / cosh(omega * ((yres - medy) / s2))) * ((yres - medy) / s2))**2 )
    tcor <- a / sqrt(b * c)
  } else {
    u <- (xres - median(xres))/s1 + (yres - median(yres))/s2
    v <- (xres - median(xres))/s1 - (yres - median(yres))/s2
    a <-  ((1 / cosh(omega * (median(abs(u))**2))) * (median(abs(u))**2)) - ((1 / cosh(omega * (median(abs(v))**2))) * (median(abs(v))**2))
    b <-  ((1 / cosh(omega * (median(abs(u))**2))) * (median(abs(u))**2)) + ((1 / cosh(omega * (median(abs(v))**2))) * (median(abs(v))**2))
    tcor <- a / b
  }
  tTaba <- ( tcor * sqrt((k - 2 - covlen) / (1 - tcor**2)) )
  if (alternative == "two.sided") {
    p <- 2*pt(-abs(tTaba), (k - 2 - covlen))
  }else{
    if (alternative == "greater") {
      p <- pt(-abs(tTaba), (k - 2 - covlen), lower.tail = TRUE)
    }else{
      p <- pt(-abs(tTaba), (k - 2 - covlen), lower.tail = FALSE)
    }
  }
  TabaC <- list(correlation = tcor,
                t.statistic = tTaba,
                p.value = p )
  return(TabaC)
}

