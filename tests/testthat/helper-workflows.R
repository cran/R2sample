# Helpers adapted from tests_twosample_routines.Rmd.

expect_valid_test_result <- function(x) {
  expect_type(x, "list")
  expect_named(x, c("statistics", "p.values"))
  expect_true(is.numeric(x$statistics))
  expect_true(is.numeric(x$p.values))
  expect_equal(length(x$statistics), length(x$p.values))
  expect_true(all(is.finite(x$statistics)))
  expect_true(all(is.finite(x$p.values)))
  expect_true(all(x$p.values >= 0 & x$p.values <= 1))
}

expect_valid_power_result <- function(x, allow_na = FALSE) {
  expect_true(is.matrix(x))
  expect_true(is.numeric(x))
  expect_gt(length(x), 0L)

  if (allow_na) {
    observed <- x[!is.na(x)]
    expect_gt(length(observed), 0L)
    expect_true(all(is.finite(observed)))
    expect_true(all(observed >= 0 & observed <= 1))
  } else {
    expect_true(all(is.finite(x)))
    expect_true(all(x >= 0 & x <= 1))
  }
}

myTSc1 <- function(x, y) {
  out <- abs(mean(x) - mean(y))
  names(out) <- "Mean"
  out
}

myTSc2 <- function(x, y, TSextra) {
  xtmp <- diff(TSextra$quant(x))
  ytmp <- diff(TSextra$quant(y))
  out <- c(abs(median(x) - median(y)), abs(xtmp - ytmp))
  names(out) <- c("Median", "IQR")
  out
}

TSextra_cont <- list(quant = function(x) quantile(x, c(0.25, 0.75)))

hybrid_cont <- function(mu) {
  x <- rnorm(40, mu)
  list(x = x, y = rnorm(50, mean(x), sd(x)))
}

rnull_cont <- function(dta) {
  list(x = dta$x, y = rnorm(length(dta$y), mean(dta$x), sd(dta$x)))
}

myTSd1 <- function(x, y, vals) {
  out <- abs(mean(x * vals) - mean(y * vals))
  names(out) <- "D"
  out
}

myTSd2 <- function(x, y, vals, TSextra) {
  xtmp <- TSextra$quant(rep(vals, x))
  ytmp <- TSextra$quant(rep(vals, y))
  out <- c(abs(median(x * vals) - median(y * vals)), abs(xtmp - ytmp))
  names(out) <- c("Median", "IQR")
  out
}

TSextra_disc <- list(quant = function(x) diff(quantile(x, c(0.25, 0.75))))

hybrid_disc <- function(p) {
  x <- table(c(0:5, rbinom(120, 5, p))) - 1
  phat <- mean(rep(0:5, x)) / 5
  list(
    x = x,
    y = table(c(0:5, rbinom(140, 5, phat))) - 1,
    vals = 0:5
  )
}

f_disc <- function(p) {
  vals <- 0:5
  x <- table(c(0:5, rbinom(120, 5, 0.5))) - 1
  y <- table(c(0:5, rbinom(140, 5, p))) - 1
  list(x = x, y = y, vals = vals)
}

rnull_disc <- function(dta) {
  phat <- mean(rep(0:5, dta$x)) / 5
  list(
    x = dta$x,
    y = table(c(0:5, rbinom(sum(dta$y), 5, phat))) - 1,
    vals = 0:5
  )
}
