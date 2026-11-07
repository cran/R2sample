test_that("run.studies handles user-supplied tests", {
  set.seed(404)
  extra_cont <- list(nbins = 5, statistic = FALSE)
  cont <- suppressMessages(
    capture.output(
      ans_cont <- run.studies(
        1:2, 1:2, TS = myTS_cont, TSextra = extra_cont,
        data_type = "cont", With.p.value = TRUE,
        SuppressMessages = TRUE, B = 20
      )
    )
  )
  expect_true(is.matrix(ans_cont))
  expect_true(is.numeric(ans_cont))
  expect_gt(length(ans_cont), 0L)
  expect_true(all(is.finite(ans_cont)))
  expect_true(all(ans_cont >= 0 & ans_cont <= 1))
  cat("\nMARK 6\n"); flush.console()
  extra_disc <- list(nbins = 5, statistic = TRUE)
  disc <- suppressMessages(
    capture.output(
      ans_disc <- run.studies(
        1:2, 1:2, TS = myTS_disc, TSextra = extra_disc,
        data_type = "disc", SuppressMessages = TRUE, B = 20
      )
    )
  )
  expect_true(is.matrix(ans_disc))
  cat("\nMARK 7\n"); flush.console()
  expect_true(is.numeric(ans_disc))
  cat("\nMARK 8\n"); flush.console()
  expect_gt(length(ans_disc), 0L)
  cat("\nMARK 9\n"); flush.console()
  expect_true(all(is.finite(ans_disc)))
  cat("\nMARK 10\n"); flush.console()
  expect_true(all(ans_disc >= 0 & ans_disc <= 1))
})