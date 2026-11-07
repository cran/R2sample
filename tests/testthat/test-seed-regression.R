test_that("bincounter does not modify its input", {
  x <- c(4, 1, 3, 2)
  x0 <- x
  R2sample:::bincounter(x, c(1, 2, 3, 4))
  expect_identical(x, x0)
})

test_that("chi_test followed by permutation is reproducible", {
  set.seed(503)
  x <- rnorm(30)
  y <- rnorm(30, 0.25)

  f <- function() {
    R2sample:::set_seed(34567)
    dta <- list(x=x, y=y)
    extra <- R2sample:::makeTSextra(samplingmethod=1)

    R2sample:::chi_test(dta, nbins=c(50,10),
                        minexpcount=5, typeTS=1)

    replicate(5, R2sample:::gen_sim_data(dta, extra),
              simplify=FALSE)
  }

  a <- f()
  runif(20)
  b <- f()
  expect_identical(a, b)
})

test_that("adjusted p value is reproducible with non-chi methods", {
  set.seed(503)
  x <- rnorm(30)
  y <- rnorm(30, 0.25)

  get_messages <- function()
    capture_messages(
      R2sample::twosample_test_adjusted_pvalue(
        x, y, B=c(19,19), seed=34567,
        SuppressMessages=TRUE, doMethods=c("KS","AD")
      )
    )

  a <- get_messages()
  runif(20)
  b <- get_messages()
  expect_identical(a,b)
})

test_that("adjusted p value is reproducible when chi methods are included", {
  set.seed(504)
  x <- rnorm(30)
  y <- rnorm(30, 0.25)

  get_messages <- function()
    capture_messages(
      R2sample::twosample_test_adjusted_pvalue(
        x, y, B=c(19,19), seed=45678,
        SuppressMessages=TRUE,
        doMethods=c("KS","ES small")
      )
    )

  a <- get_messages()
  runif(20)
  b <- get_messages()
  expect_identical(a,b)
})
