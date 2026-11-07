test_that("makeTSextra replaces control fields rather than duplicating them", {
  rn <- function() 1
  z <- makeTSextra(list(samplingmethod=99, keep=7), 2, rn)
  expect_identical(z$samplingmethod, 2)
  expect_identical(z$rnull, rn)
  expect_identical(z$keep, 7)
  expect_equal(sum(names(z)=="samplingmethod"), 1L)
  expect_equal(sum(names(z)=="rnull"), 1L)
})

test_that("maketypeTS dispatches all supported interfaces", {
  expect_identical(maketypeTS(Continuous=TRUE, WithWeights=FALSE)$typeTS, 1L)
  expect_identical(maketypeTS(Continuous=TRUE, WithWeights=TRUE)$typeTS, 2L)
  expect_identical(maketypeTS(Continuous=FALSE, WithWeights=FALSE)$typeTS, 5L)

  c2 <- function(x,y) c(A=0)
  c3 <- function(x,y,TSextra) c(A=0)
  d3 <- function(x,y,vals) c(A=0)
  d4 <- function(x,y,vals,TSextra) c(A=0)

  expect_identical(maketypeTS(c2, TRUE, FALSE)$typeTS, 3L)
  expect_identical(maketypeTS(c3, TRUE, FALSE)$typeTS, 4L)
  expect_identical(maketypeTS(d3, FALSE, FALSE)$typeTS, 6L)
  expect_identical(maketypeTS(d4, FALSE, FALSE)$typeTS, 7L)

  expect_error(maketypeTS(function(x) 0, TRUE, FALSE), "2 or 3")
  expect_error(maketypeTS(function(x,y) 0, FALSE, FALSE), "3 or 4")
})

test_that("set_seed validates seeds", {
  expect_invisible(set_seed(NULL))
  expect_error(set_seed(NA_real_), "seed")
  expect_error(set_seed(c(1,2)), "seed")
  expect_error(set_seed(Inf), "seed")
})
