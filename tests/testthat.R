library(testthat)
library(Taba)

########################################    Taba
test_that("Taba 1", {
  x = rnorm(100)
  y = rnorm(100)
  expect_error(taba(x,y,method = "bghsjxcmv"),"invalid 'method' argument")
})
test_that("Taba 2", {
  x = rnorm(100)
  y = rnorm(100)
  expect_error(taba(x,y,omega = 10),"'omega' must be between 0 and 1")
})
test_that("Taba 3", {
  x = "test"
  y = NULL
  expect_error(taba(x,y),"supply both 'x' and 'y' or a matrix-like 'x'")
})
test_that("Taba 4", {
  x = rnorm(100)
  y = rnorm(100)
  x = factor(x)
  expect_error(taba(x,y),"'x' must be numeric")
})
test_that("Taba 5", {
  x = rnorm(100)
  y = factor(rnorm(100))
  expect_error(taba(x,y),"'y' must be numeric")
})
test_that("Taba 6", {
  x = rnorm(100)
  y = rnorm(100)
  x[5] = NA
  y[72] = NA
  expect_warning(taba(x,y),"Missing data included in dataset was removed row-wise. Results may not be accurate.")
})
test_that("Taba 7", {
  x = rnorm(100)
  y = rnorm(99)
  expect_error(taba(x,y),"'x' and 'y' must have the same length")
})


##########################################          Taba Test
test_that("TabaTest 1", {
  x = rnorm(100)
  y = rnorm(100)
  expect_error(taba.test(x,y,method = "bghsjxcmv"),"invalid 'method' argument")
})
test_that("TabaTest 2", {
  x = rnorm(100)
  y = rnorm(100)
  expect_error(taba.test(x,y,alternative = "sjhdbfurfh" ),"invalid 'alternative' argument")
})
test_that("TabaTest 3", {
  x = rnorm(100)
  y = rnorm(100)
  expect_error(taba.test(x,y,omega = 10),"'omega' must be between 0 and 1")
})
test_that("TabaTest 4", {
  x = "test"
  y = NULL
  expect_error(taba.test(x,y),"supply both 'x' and 'y' or a matrix-like 'x'")
})
test_that("TabaTest 5", {
  x = as.character(rnorm(100))
  y = rnorm(100)
  expect_error(taba.test(x,y),"'x' must be numeric")
})
test_that("TabaTest 6", {
  x = rnorm(100)
  y = factor(rnorm(100))
  expect_error(taba.test(x,y),"'y' must be numeric")
})
test_that("TabaTest 7", {
  x = rnorm(100)
  y = rnorm(100)
  x[5] = NA
  y[72] = NA
  expect_warning(taba.test(x,y),"Missing data included in dataset was removed row-wise. Results may not be accurate.")
})
test_that("TabaTest 8", {
  x = rnorm(100)
  y = rnorm(99)
  expect_error(taba.test(x,y),"'x' and 'y' must have the same length")
})


############################################          Taba Partial
test_that("TabaPartial 1", {
  x = rnorm(100)
  y = rnorm(100)
  z = rnorm(100)
  expect_error(taba.partial(x,y,z,method = "bghsjxcmv"),"invalid 'methods' argument")
})
test_that("Tabapartial 2", {
  x = rnorm(100)
  y = rnorm(100)
  z = rnorm(100)
  expect_error(taba.partial(x,y,z,alternative = "sjhdbfurfh" ),"invalid 'alternative' argument")
})
test_that("Tabapartial 3", {
  x = rnorm(100)
  y = rnorm(100)
  z = rnorm(100)
  expect_error(taba.partial(x,y,z,semi = "sjhdbfurfh" ),"invalid 'semi' argument")
})
test_that("TabaPartial 4", {
  x = rnorm(100)
  y = rnorm(100)
  z = rnorm(100)
  expect_error(taba.partial(x,y,z,omega = 10),"'omega' must be between 0 and 1")
})
test_that("TabaPartial 5", {
  x = "test"
  y = NULL
  z = rnorm(100)
  expect_error(taba.partial(x,y,z),"supply both 'x' and 'y' or a matrix-like 'x'")
})
test_that("TabaPartial 6", {
  x = factor(rnorm(100))
  y = rnorm(100)
  z = rnorm(100)
  expect_error(taba.partial(x,y,z),"'x' must be numeric")
})
test_that("TabaPartial 7", {
  x = rnorm(100)
  y = as.character(rnorm(100))
  z = rnorm(100)
  expect_error(taba.partial(x,y,z),"'y' must be numeric")
})
test_that("TabaPartial 8", {
  x = rnorm(100)
  y = rnorm(100)
  z = rnorm(100)
  x[5] = NA
  y[72] = NA
  expect_warning(taba.partial(x,y,z),"Missing data included in dataset was removed row-wise. Results may not be accurate.")
})
test_that("TabaPartial 9", {
  x = rnorm(100)
  y = rnorm(99)
  z = rnorm(100)
  expect_error(taba.partial(x,y,z),"'x','y', and 'Covariares' must have the same length")
})
test_that("TabaPartial 10", {
  x = rnorm(100)
  y = rnorm(100)
  expect_error(taba.partial(x,y),"No covariates entered")
})
test_that("TabaPartial 10", {
  x = rnorm(100)
  y = rnorm(100)
  z = factor(rnorm(100))
  expect_error(taba.partial(x,y,z),"All covariates must be numeric")
})


############################################          Taba Matrix
test_that("TabaMatrix 1", {
  x = rnorm(100)
  y = rnorm(100)
  z = rnorm(100)
  expect_error(taba.matrix(x,y,z,method = "bghsjxcmv"),"invalid 'methods' argument")
})
test_that("TabaMatrix 2", {
  x = rnorm(100)
  y = rnorm(100)
  z = rnorm(100)
  expect_error(taba.matrix(x,y,z,alternative = "sjhdbfurfh" ),"invalid 'alternative' argument")
})
test_that("TabaMatrix 3", {
  x = rnorm(100)
  y = rnorm(100)
  z = rnorm(100)
  expect_error(taba.matrix(x,y,z,omega = 10),"'omega' must be between 0 and 1")
})
test_that("TabaMatrix 4", {
  x = factor(rnorm(100))
  y = rnorm(100)
  z = rnorm(100)
  expect_error(taba.matrix(x,y,z),"'x' must be numeric")
})
test_that("TabaMatrix 5", {
  x = rnorm(100)
  y = as.character(rnorm(100))
  z = rnorm(100)
  expect_error(taba.matrix(x,y,z),"'y' must be numeric")
})
test_that("TabaMatrix 6", {
  x = rnorm(100)
  y = rnorm(100)
  z = rnorm(100)
  x[5] = NA
  y[72] = NA
  expect_warning(taba.matrix(x,y,z),"Missing data included in dataset was removed row-wise. Results may not be accurate.")
})
test_that("TabaMatrix 7", {
  x = rnorm(100)
  y = rnorm(99)
  z = rnorm(100)
  expect_error(taba.matrix(x,y,z),"all vectors must have the same length")
})
test_that("TabaMatrix 8", {
  x = rnorm(100)
  y = rnorm(100)
  z = factor(rnorm(100))
  expect_error(taba.matrix(x,y,z),"All vectors must be numeric")
})



################################################  Compare between Taba
test_that("Compare 1", {
  x = rnorm(100)
  y = rnorm(100)
  expect_equal(taba(x,y),as.numeric(taba.matrix(y,x)$cmatrix[2,1]))
})
test_that("Compare 2", {
  x = rnorm(100)
  y = rnorm(100)
  expect_equal(taba(x,y),1-as.numeric(taba.matrix(x,y)$distance[1])*2)
})
test_that("Compare 3", {
  x = rnorm(100)
  y = rnorm(100)
  expect_equal(taba(x,y,method = "tabarank", omega = .3),
               taba(y,x,method = "tabarank", omega = .3))
})
test_that("Compare 4", {
  x = rnorm(100)
  y = rnorm(100)
  expect_equal(taba(x,y,method = "tabarank"),
               taba(y,x,method = "tabarank", omega = .45))
})
test_that("Compare 5", {
  x = rnorm(100)
  y = rnorm(100)
  expect_equal(as.numeric(taba.matrix(x,y,method = "tabarank", omega = .3)$cmatrix[2,1]),
               1-as.numeric(taba.matrix(x,y,method = "tabarank", omega = .3)$distance[1])*2)
})

test_that("Compare 6", {
  x = rt(100,2)
  y = rt(100,20)
  z = rnorm(100,10,.5)
  w = rnorm(100,-3,1.4)
  w1 = rt(100,33)
  w2 = rnorm(100,-100,.1)
  e1 = lm(x~z+w)$residuals
  e2 = lm(y~w+w1+w2)$residuals
  expect_equal(taba(e1,e2),
               taba.gpartial(x,y,xcov = cbind(z,w),ycov = cbind(w,w1,w2))$correlation)
})

#  expect_output() prints specified output?
#  expect_message() displays specified message?
#  expect_warning() displays specified warning?
#  expect_error()
