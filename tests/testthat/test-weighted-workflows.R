test_that("weighted continuous workflows run", {
  set.seed(303)
  w <- function(x) dnorm(x) / dt(x, 3)

  x <- rnorm(40)
  wx <- rep(1, length(x))
  y <- rt(50, 3)
  wy <- w(y)

  out <- twosample_test(x, y, wx = wx, wy = wy, B = 19,
                        maxProcessor = 1, SuppressMessages = TRUE)
  expect_valid_test_result(out)

  f <- function(m) {
    x <- rnorm(40, m)
    wx <- rep(1, length(x))
    pool <- rt(500, 3)
    pool <- pool[abs(pool) < 5]
    y <- pool[seq_len(50)]
    wy <- dnorm(y) / dt(y, 3)
    list(x = x, y = y, wx = wx, wy = wy)
  }

  pwr <- twosample_power(f, c(0, 0.5), B = 20,
                         maxProcessor = 1, SuppressMessages = TRUE)
  expect_valid_power_result(pwr)
  expect_equal(nrow(pwr), 2L)
})
