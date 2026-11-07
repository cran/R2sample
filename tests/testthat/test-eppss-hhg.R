# Checks of the Epps-Singleton and HHG tests

test_that("Epps-Singleton agrees with scipy.stats.epps_singleton_2samp", {
  x <- round(sin((1:60) * 1.7) * 2 + (1:60)/30, 6)
  y <- round(cos((1:80) * 1.3) * 2.5, 6)
  # reference values from scipy 1.18.1
  expect_equal(unname(R2sample:::EppsS(x, y)["statistic"]), 35.02613354940631,
               tolerance = 1e-10)
  # small samples use the small sample correction
  expect_equal(unname(R2sample:::EppsS(x[1:20], y[1:15])["statistic"]),
               3.0085139027478904, tolerance = 1e-10)
  expect_equal(unname(R2sample:::TS_cont(x, y)["EppsS"]),
               unname(R2sample:::EppsS(x, y)["statistic"]), tolerance = 1e-12)
})

test_that("discrete Epps-Singleton equals continuous one on the expanded data", {
  vals <- c(1, 2, 3, 5, 8); x <- c(3, 0, 4, 2, 5); y <- c(1, 6, 2, 3, 1)
  d <- R2sample:::EppsS(x, y, vals)
  c0 <- R2sample:::EppsS(rep(vals, x), rep(vals, y))
  expect_equal(unname(d), unname(c0), tolerance = 1e-10)
})

test_that("HHG statistics are averages over all partitions", {
  lr <- function(X, Y, nx, ny) {
    T <- X + Y; N <- nx + ny
    s <- 0
    if (X > 0) s <- s + X * log(X * N/(nx * T))
    if (Y > 0) s <- s + Y * log(Y * N/(ny * T))
    2 * s
  }
  brute <- function(x, y, mmax) {
    A <- length(x); nx <- sum(x); ny <- sum(y)
    sapply(2:mmax, function(m) {
      cuts <- utils::combn(1:(A - 1), m - 1, simplify = FALSE)
      mean(sapply(cuts, function(cc) {
        b <- c(0, cc, A)
        sum(sapply(1:m, function(j) {
          I <- (b[j] + 1):b[j + 1]
          lr(sum(x[I]), sum(y[I]), nx, ny)
        }))
      }))
    })
  }
  set.seed(4)
  for (i in 1:5) {
    x <- sample(1:6, 7, replace = TRUE); y <- sample(1:6, 7, replace = TRUE)
    expect_equal(as.numeric(R2sample:::HHG_stats_disc(x, y, 5)),
                 brute(x, y, 5), tolerance = 1e-10)
  }
})

test_that("large sample p values of EppsS and HHG are valid probabilities", {
  set.seed(5)
  out <- R2sample::twosample_test(rnorm(300), rnorm(300), UseLargeSample = TRUE,
                                  SuppressMessages = TRUE)
  p <- out$p.values[c("EppsS", "HHG")]
  expect_true(all(p >= 0 & p <= 1))
})
