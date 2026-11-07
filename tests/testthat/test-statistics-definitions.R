# Checks of the built-in test statistics against their published definitions.

ref_stats <- function(x, y) {
  x <- sort(x); y <- sort(y); n <- length(x); m <- length(y); N <- n + m
  z <- sort(c(x, y)); D <- ecdf(x)(z) - ecdf(y)(z); j <- 1:(N - 1)
  c(KS = max(abs(D)),
    Kuiper = max(D, 0) + max(-D, 0),
    CvM = n * m / N^2 * sum(D^2),                       # Anderson (1962)
    AD = n * m * sum(D[-N]^2 / (j * (N - j))),          # Pettitt (1976)
    Wassp1 = sum(abs(D[-N]) * diff(z)))
}

ref_zhang <- function(x, y) {
  x <- sort(x); y <- sort(y); n <- length(x); m <- length(y); N <- n + m
  z <- sort(c(x, y)); k <- 1:N
  Fn <- (ecdf(x)(z) * n - 0.5 * (z %in% x)) / n
  Gm <- (ecdf(y)(z) * m - 0.5 * (z %in% y)) / m
  H <- (k - 0.5) / N
  e <- function(p) ifelse(p > 0, p * log(p), 0) + ifelse(p < 1, (1 - p) * log(1 - p), 0)
  rk <- rank(c(x, y)); r <- sort(rk[1:n]); s <- sort(rk[-(1:n)])
  c(ZA = -sum((n * e(Fn) + m * e(Gm)) / ((k - 0.5) * (N - k + 0.5))),
    ZK = max(n * e(Fn) + m * e(Gm) - N * e(H)),
    ZC = (sum(log(n / (1:n - 0.5) - 1) * log(N / (r - 0.5) - 1)) +
          sum(log(m / (1:m - 0.5) - 1) * log(N / (s - 0.5) - 1))) / N)
}

test_that("continuous statistics agree with published definitions", {
  set.seed(1)
  for (i in 1:50) {
    x <- rnorm(sample(2:30, 1)); y <- rnorm(sample(2:30, 1), 0.3)
    a <- R2sample:::TS_cont(x, y)
    b <- ref_stats(x, y)
    expect_equal(unname(a[names(b)]), unname(b), tolerance = 1e-8)
    zh <- ref_zhang(x, y)
    # ZA and ZC are returned with the sign changed
    expect_equal(unname(a[c("ZA", "ZK", "ZC")]),
                 unname(c(-zh["ZA"], zh["ZK"], -zh["ZC"])), tolerance = 1e-6)
  }
  # example where the old Kuiper initialisation was wrong
  expect_equal(unname(R2sample:::TS_cont(c(0, 2, 3), c(1, 1.5, 4, 5))["Kuiper"]), 2/3)
  expect_false("LR" %in% names(R2sample:::TS_cont(1:5, 2:7)))
})

test_that("discrete statistics equal continuous ones on the same data", {
  set.seed(2)
  for (i in 1:30) {
    vals <- sort(unique(round(rnorm(8), 2)))
    k <- length(vals)
    x <- sample(0:5, k, replace = TRUE); y <- sample(1:5, k, replace = TRUE)
    if (sum(x) == 0) x[1] <- 1
    dta <- list(x = x, y = y, vals = vals)
    d <- R2sample:::TS_disc(as.integer(x), as.integer(y), vals, R2sample:::weights(dta))
    c0 <- R2sample:::TS_cont(rep(vals, x), rep(vals, y))
    expect_equal(unname(d), unname(c0[names(d)]), tolerance = 1e-8)
  }
})

test_that("weighted statistics reduce to unweighted ones and use the weights", {
  set.seed(3)
  x <- rnorm(20); y <- rnorm(25)
  a <- R2sample:::TSw_cont(x, y, rep(1, 20), rep(1, 25))
  b <- R2sample:::TS_cont(x, y)
  expect_equal(unname(a), unname(b[names(a)]), tolerance = 1e-10)
  wx <- runif(20, 0.5, 2); wy <- runif(25, 0.5, 2)
  dta <- list(x = x, y = y, wx = wx, wy = wy)
  expect_equal(R2sample:::calcTS(dta, R2sample:::TSw_cont, 2L, list()),
               R2sample:::TSw_cont(x, y, wx, wy))
})

test_that("discrete permutation sampling is valid", {
  set.seed(4)
  x <- c(1, 1, 2, 1, 30, 1, 1, 2, 1, 1); y <- c(1, 0, 1, 1, 25, 1, 2, 0, 1, 1)
  for (method in 1:2) {
    for (i in 1:200) {
      s <- R2sample:::gen_disc(x, y, 1:10, list(samplingmethod = method))
      expect_true(all(s$x >= 0 & s$y >= 0))
      expect_equal(s$x + s$y, x + y)
      expect_equal(sum(s$x), sum(x))
    }
  }
})

test_that("large sample p values are valid probabilities", {
  for (z in c(0.05, 0.3473, 0.46136, 0.74346, 1.2, 5)) {
    p <- R2sample:::asymptotic_pvalues(c(0.01, 0.02, z, 1), 1e4, 1e4)
    expect_true(all(p >= 0 & p <= 1))
  }
  expect_equal(unname(R2sample:::asymptotic_pvalues(c(0.01, 0.02, 0.46136, 1), 1e4, 1e4)["CvM"]),
               0.05, tolerance = 1e-4)
})

test_that("power_newtest uses all parameter values", {
  set.seed(5)
  f <- function(mu) list(x = rnorm(30), y = rnorm(30, mu))
  TSp <- function(x, y) {o <- ks.test(x, y)$p.value; names(o) <- "KSp"; o}
  pw <- R2sample:::power_newtest(TSp, f, xparam = c(0, 2), yparam = c(0, 0),
                                 TSextra = list(), B = 50)
  expect_gt(pw[2, 1], 0.9)
})
