test_that("discrete two-sample workflows run", {
  set.seed(202)
  dta <- f_disc(0.5)
  x <- dta$x
  y <- dta$y
  vals <- dta$vals

  out_list <- twosample_test(dta, B = 19, maxProcessor = 1,
                             SuppressMessages = TRUE, doMethods = "default")
  expect_valid_test_result(out_list)

  out_xy <- twosample_test(x, y, vals, B = 19, maxProcessor = 1,
                           SuppressMessages = TRUE, doMethods = c("KS", "AD"))
  expect_valid_test_result(out_xy)
  expect_named(out_xy$p.values, c("KS", "AD"))

  out_null <- twosample_test(x, y, vals, rnull = rnull_disc, B = 19,
                             maxProcessor = 1, SuppressMessages = TRUE,
                             doMethods = c("KS", "AD"))
  expect_valid_test_result(out_null)

  out_custom1 <- twosample_test(x, y, vals, TS = myTSd1, B = 19,
                                maxProcessor = 1, SuppressMessages = TRUE)
  expect_valid_test_result(out_custom1)
  expect_named(out_custom1$statistics, "D")

  out_custom2 <- twosample_test(x, y, vals, TS = myTSd2,
                                TSextra = TSextra_disc, B = 19,
                                maxProcessor = 1, SuppressMessages = TRUE)
  expect_valid_test_result(out_custom2)
  expect_named(out_custom2$statistics, c("Median", "IQR"))

  adjusted_messages <- capture_messages(
    adjusted <- twosample_test_adjusted_pvalue(
      dta, B = c(19, 19), SuppressMessages = TRUE,
      doMethods = c("KS", "AD")
    )
  )
  expect_type(adjusted, "list")
  expect_true(adjusted$adjusted.p.value >= 0 && adjusted$adjusted.p.value <= 1)
  expect_true(any(grepl("adjusted p value", adjusted_messages,
                        ignore.case = TRUE)))

  # Keep power checks deliberately small.  The R Markdown validation file
  # performs the comprehensive power studies; testthat only needs to verify
  # that the discrete power workflows execute and return probabilities.
  pwr_custom <- twosample_power(f_disc, 0.6, TS = myTSd1,
                                B = 10, maxProcessor = 1,
                                SuppressMessages = TRUE)
  expect_valid_power_result(pwr_custom)

  pwr_hybrid <- twosample_power(hybrid_disc, 0.6, TS = myTSd1,
                                rnull = rnull_disc, B = 10,
                                maxProcessor = 1, SuppressMessages = TRUE)
  expect_valid_power_result(pwr_hybrid, allow_na = TRUE)
})
