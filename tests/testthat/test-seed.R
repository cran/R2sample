test_that("twosample_test is reproducible with seed", {
  set.seed(501)
  x <- rnorm(35)
  y <- rnorm(40, 0.3)

  a <- twosample_test(x, y, B=19, seed=12345, maxProcessor=1,
                      SuppressMessages=TRUE, doMethods=c("KS","AD"))
  runif(20)
  b <- twosample_test(x, y, B=19, seed=12345, maxProcessor=1,
                      SuppressMessages=TRUE, doMethods=c("KS","AD"))
  expect_identical(a,b)
})

test_that("twosample_power is reproducible with seed", {
  f <- function(mu) list(x=rnorm(30), y=rnorm(35,mu))
  a <- twosample_power(f, c(0,0.5), B=20, seed=23456,
                       maxProcessor=1, SuppressMessages=TRUE)
  runif(20)
  b <- twosample_power(f, c(0,0.5), B=20, seed=23456,
                       maxProcessor=1, SuppressMessages=TRUE)
  expect_identical(a,b)
})

test_that("twosample_test_adjusted_pvalue is reproducible with seed", {
  set.seed(503)
  x <- rnorm(30)
  y <- rnorm(30,0.25)

  get_messages <- function()
    capture_messages(
      twosample_test_adjusted_pvalue(
        x, y, B=c(19,19), seed=34567,
        SuppressMessages=TRUE, doMethods=c("KS","AD")
      )
    )

  a <- get_messages()
  runif(20)
  b <- get_messages()

  expect_identical(a, b)
})

test_that("public simulation functions reject invalid seeds", {
  x <- 1:10
  y <- 2:11
  f <- function(mu) list(x=rnorm(10),y=rnorm(10,mu))

  expect_error(twosample_test(x,y,B=5,seed=NA_real_,maxProcessor=1),"seed")
  expect_error(twosample_power(f,0,B=5,seed=Inf,maxProcessor=1),"seed")
  expect_error(
    twosample_test_adjusted_pvalue(
      x,y,B=c(5,5),seed=c(1,2),doMethods="KS"
    ),
    "seed"
  )
})
