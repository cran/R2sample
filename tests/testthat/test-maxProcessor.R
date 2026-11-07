test_that("makemaxProcessor respects forced and explicit single-processor choices", {
  dta <- list(x=1:5, y=6:10)
  TS <- function(x,y) c(A=mean(x)-mean(y))

  expect_identical(
    makemaxProcessor(4, dta, TS, 3L, list(), 100,
                     useSingleProcessor=TRUE),
    1L
  )
  expect_identical(
    makemaxProcessor(2, dta, TS, 3L, list(), 100),
    2L
  )
})

test_that("makemaxProcessor validates explicit maxProcessor", {
  dta <- list(x=1:5, y=6:10)
  TS <- function(x,y) c(A=0)

  expect_error(makemaxProcessor(0, dta, TS, 3L, list(), 100),
               "maxProcessor")
  expect_error(makemaxProcessor(1.5, dta, TS, 3L, list(), 100),
               "maxProcessor")
  expect_error(makemaxProcessor(NA_real_, dta, TS, 3L, list(), 100),
               "maxProcessor")
})

test_that("makemaxProcessor uses timing rule when maxProcessor is missing", {
  dta <- list(x=1:5, y=6:10)
  TS <- function(x,y) c(A=0)

  local_mocked_bindings(
    timecheck=function(...) 0.001,
    .package="R2sample"
  )

  expect_identical(
    makemaxProcessor(dta=dta, TS=TS, typeTS=3L, TSextra=list(),
                     B=100, SuppressMessages=TRUE),
    1L
  )
})

test_that("timecheck supports discrete typeTS through calcTS", {
  dta <- list(x=c(10L,20L), y=c(12L,18L), vals=c(1,2))
  TS <- function(x,y,vals) c(A=sum(abs(x-y)))
  tm <- timecheck(dta, TS, 6L, list())
  expect_true(is.numeric(tm))
  expect_length(tm, 1L)
  expect_true(is.finite(tm))
  expect_gte(tm, 0)
})
