#'
#' Robust Correlation Matrix
#'
#' @description Calculates a correlation, distance, and p-value matrix using one of the
#'     specified robust methods Taba linear or Taba rank correlation.
#' @usage taba.matrix(x, y = NULL, ..., method = c("taba", "tabarank"),
#'             alternative = c("less", "greater", "two.sided"),
#'             omega = 0.45)
#' @param x A numeric vector of length greater than 2 must be same length as all other vectors.
#' @param y A numeric vector of length greater than 2 must be same length as all other vectors.
#' @param ... Numeric vector(s) of length equal to x and y. May be of class matrix
#'    or data.frame, whose columns will be compared and whose column's length must be of
#'    equal length to x and y. Not one vector or column name can be "x" or "y."
#' @param method A character string of \code{"taba"} or \code{"tabarank"}
#'    determining if one wants to calculate Taba linear or Taba rank (monotonic) correlation,
#'    respectively. If no method is specified, the function will output Taba
#'    Linear correlation.
#' @param alternative Character string specifying the alternative hypothesis must be one
#'    of \strong{\code{"less"}} for negative association, \code{"greater"} for
#'    positive association, or \code{"two.sided"} for difference in association.
#'    If the alternative is not specified, the function will default to a two sided test.
#' @param omega Numeric allowing the user to alter the tuning constant. If one is not specified,
#'    the function will default to 0.45. Range is between 0 and 1.
#' @details This function uses Taba linear or Taba rank (monotonic) correlation to
#'    calculate the association of two or more numeric vectors. Numeric vectors under \code{...}
#'    are combined colomn-wise with x and y. When inserting a single matrix x, the function will
#'    calculate the correlation matix using the columns of matrix x. \cr
#'    Matricies or data frames with numeric cells can be inserted in \code{...}, whereby
#'    each column in the matrix or data frame will be treated as a different vector
#'    for comparison. Columns must all have different names from each other. No vector
#'    or column should be named "x" or "y," as these refer to the first two vectors respectively,
#'    if inserted as a vector or matrix with no name. Missing values in any of the vectors
#'    are deleted row-wise. \cr
#'    The default for this function is a two sided test using Taba linear partial correlation,
#'    with the tuning constant \code{omega} equal to 0.45.
#' @return This function returns the robust linear or monotonic association
#'   between two or more numeric vectors, as a matrix; the distance matrix, as type dist;
#'   and a p-value matrix corresponding to the correlation matrix.
#' @seealso
#'   \code{\link{taba}} for calculating Taba linear or Taba rank (monotonic) correlations
#'   \cr\code{\link{taba.test}} for testing Taba linear or Taba rank (monotonic) correlations
#'   \cr\code{\link{taba.gpartial}} for generalized partial correlations
#'   \cr\code{\link{taba.partial}} for partial and semipartial correlations
#' @references The paper is under review for possible publication.
#' @examples
#' x = rnorm(100)
#' y = rnorm(100)
#' z1 = rnorm(100)
#' z2 = rnorm(100)
#' z3 = rnorm(100)
#' Z = cbind(z1,z3)
#' colnames(Z) = c("A","B")
#' taba.matrix(x, y, z1, z2, z3, method = "tabarank")
#' taba.matrix(x, y, z2, Z, alternative = "less", omega = 0.4)
#' taba.matrix(Z, method = "tabarank")
#' @import robustbase
#'         stats
#' @export taba.matrix

taba.matrix = function(x, y = NULL, ..., method = c("taba", "tabarank"),
                       alternative = c("less", "greater", "two.sided"),
                       omega = 0.45) {
  if (missing(method)) {
    method <- "taba"
  }
  na.method <- pmatch(method, c("taba", "tabarank"))
  if (is.na(na.method)) {
    stop("invalid 'methods' argument")
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
    x <- as.data.frame(x)
  if (!(is.numeric(x[,1]) || is.logical(x[,1]))) {
    stop("'x' must be numeric")
    stopifnot(is.atomic(x[,1]))
  }
  if (!is.null(y)) { #Inserted
      y <- as.data.frame(y)
    if (!(is.numeric(y[,1]) || is.logical(y[,1])))
      stop("'y' must be numeric")
      stopifnot(is.atomic(y[,1]))
    if (missing(...)) {
      if (sum(is.na(x)) > 0 || sum(is.na(y)) > 0) {
        warning("Missing data included in dataset was removed row-wise. Results may not be accurate.")
        miss <- which(complete.cases(x,y) == FALSE)
        x <- x[-miss,]
        y <- y[-miss,]
      }
      frame <- as.matrix(cbind(x,y))
      n <- ncol(frame)
      Tab.x <- matrix(nrow = n, ncol = n)
      pmatrix <- Tab.x
    }else{
      Vectors <- cbind.data.frame(...)
      if ((length(x[,1]) != length(y[,1])) || (length(x[,1]) != length(Vectors[,1]))) {
        stop("all vectors must have the same length")
      }
      if (sum(sapply(Vectors,is.numeric)) != length(Vectors)) {
        stop("All vectors must be numeric")
        stopifnot(is.atomic(y))
      }
      if (sum(is.na(x[,1])) > 0 || sum(is.na(y[,1])) > 0 || sum(is.na(Vectors)) > 0) {
        warning("Missing data included in dataset was removed row-wise. Results may not be accurate.")
        miss <- which(complete.cases(x,y,Vectors) == FALSE)
        x <- x[-miss,]
        y <- y[-miss,]
        Vectors <- Vectors[-miss,]
      }
      frame <- as.matrix(cbind(x,y,Vectors))
      n <- ncol(frame)
      Tab.x <- matrix(nrow = n, ncol = n)
      pmatrix <- Tab.x
    }
  }else{
    if (sum(is.na(x)) > 0) {
      warning("Missing data included in dataset was removed row-wise. Results may not be accurate.")
      miss <- which(complete.cases(x) == FALSE)
      x <- x[-miss,]
    }
    frame <- as.matrix(x)
    n <- ncol(frame)
    Tab.x <- matrix(nrow = n, ncol = n)
    pmatrix <- Tab.x
  }
  Tab = function(x, y, method, alternative, omega) {
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
    lenx <- length(x)
    tTaba <- tcor * sqrt((lenx - 2)/(1 - tcor**2))
    if (alternative == "two.sided") {
      p <- 2*pt(-abs(tTaba), (lenx - 2))
    }else{
      if (alternative == "greater") {
        p <- pt(-abs(tTaba), (lenx - 2), lower.tail = TRUE)
      }else{
        p <- pt(-abs(tTaba), (lenx - 2), lower.tail = FALSE)
      }
    }
    TabaC <- list(correlation = tcor,
                  t.statistic = tTaba,
                  p.value = p )
    return(TabaC)
  }

  for (j in 1:(n-1)) {
    for (i in (j+1):n) {
      tmp <- Tab(frame[,i],frame[,j], method = method, alternative = alternative, omega = omega)
      Tab.x[i,j] <- tmp$correlation
      pmatrix[i,j] <- tmp$p.value
    }
  }

diag(Tab.x) <- 1
diag(pmatrix) <- 1

Tab.x <- as.data.frame(Tab.x)
pmatrix <- as.data.frame(pmatrix)
colnames(Tab.x) <- colnames(frame)
rownames(Tab.x) <- colnames(frame)
colnames(pmatrix) <- colnames(frame)
rownames(pmatrix) <- colnames(frame)

#delete any rows with missing values
while (sum(is.nan(Tab.x[,1])) != 0) {
  for (i in 1:nrow(Tab.x)) {
    if ((nrow(Tab.x)-1) == sum(is.nan(Tab.x[,1]), na.rm = F)) {
      Tab.x <- Tab.x[-1,]
      Tab.x <- Tab.x[,-1]
      break
    }
    if (is.nan(Tab.x[i,1])) {
      Tab.x <- Tab.x[-i,]
      Tab.x <- Tab.x[,-i]
      break
    }
  }
}
while (sum(is.nan(pmatrix[,1])) != 0) {
  for (i in 1:nrow(pmatrix)) {
    if ((nrow(pmatrix)-1) == sum(is.nan(pmatrix[,1]), na.rm = F)) {
      pmatrix <- pmatrix[-1,]
      pmatrix <- pmatrix[,-1]
      break
    }
    if (is.nan(Tab.x[i,1])) {
      pmatrix <- pmatrix[-i,]
      pmatrix <- pmatrix[,-i]
      break
    }
  }
}
#Distance Matrix
Distance <- (1-as.dist(Tab.x))/2

#complete correlation matrix
  for (j in 1:(n-1)) {
  for (i in (j+1):n) {
    Tab.x[j,i] <- Tab.x[i,j]
    pmatrix[j,i] <- pmatrix[i,j]
  }
}

CorList <- list("cmatrix" = Tab.x, "distance" = Distance, "pmatrix" = pmatrix)
return(CorList)
}
