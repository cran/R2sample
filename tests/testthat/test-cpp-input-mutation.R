test_that("continuous built-in test statistics do not modify inputs", {
  x <- c(3, 1, 4, 2)
  y <- c(8, 6, 7, 5)
  x0 <- x
  y0 <- y

  TS_cont(x, y)

  expect_identical(x, x0)
  expect_identical(y, y0)
})


test_that("weighted continuous built-in test statistics do not modify inputs", {
  x <- c(3, 1, 4, 2)
  y <- c(8, 6, 7, 5)
  wx <- c(0.1, 0.2, 0.3, 0.4)
  wy <- c(0.4, 0.3, 0.2, 0.1)

  x0 <- x
  y0 <- y
  wx0 <- wx
  wy0 <- wy

  TSw_cont(x, y, wx, wy)

  expect_identical(x, x0)
  expect_identical(y, y0)
  expect_identical(wx, wx0)
  expect_identical(wy, wy0)
})


test_that("discrete built-in test statistics do not modify inputs", {
  x <- c(4L, 2L, 5L, 3L)
  y <- c(3L, 6L, 2L, 4L)
  vals <- c(1, 2, 3, 4)
  adw <- rep(1, 4)

  x0 <- x
  y0 <- y
  vals0 <- vals
  adw0 <- adw

  TS_disc(x, y, vals, adw)

  expect_identical(x, x0)
  expect_identical(y, y0)
  expect_identical(vals, vals0)
  expect_identical(adw, adw0)
})
