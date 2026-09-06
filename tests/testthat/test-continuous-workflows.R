test_that("continuous two-sample workflows run", {
  set.seed(101)
  dta <- hybrid_cont(0)
  x <- dta$x
  y <- dta$y

  out_list <- twosample_test(dta, B = 19, maxProcessor = 1,
                             SuppressMessages = TRUE, doMethods = "default")
  expect_valid_test_result(out_list)

  out_xy <- twosample_test(x, y, B = 19, maxProcessor = 1,
                           SuppressMessages = TRUE, doMethods = c("KS", "AD"))
  expect_valid_test_result(out_xy)
  expect_named(out_xy$p.values, c("KS", "AD"))

  out_null <- twosample_test(x, y, rnull = rnull_cont, B = 19,
                             maxProcessor = 1, SuppressMessages = TRUE,
                             doMethods = c("KS", "AD"))
  expect_valid_test_result(out_null)

  out_custom1 <- twosample_test(x, y, TS = myTSc1, B = 19,
                                maxProcessor = 1, SuppressMessages = TRUE)
  expect_valid_test_result(out_custom1)
  expect_named(out_custom1$statistics, "Mean")

  out_custom2 <- twosample_test(x, y, TS = myTSc2, TSextra = TSextra_cont,
                                B = 19, maxProcessor = 1,
                                SuppressMessages = TRUE)
  expect_valid_test_result(out_custom2)
  expect_named(out_custom2$statistics, c("Median", "IQR"))

  adjusted_messages <- capture_messages(
    adjusted <- twosample_test_adjusted_pvalue(
      dta, B = c(19, 19), SuppressMessages = TRUE,
      doMethods = c("KS", "AD")
    )
  )
  expect_null(adjusted)
  expect_true(any(grepl("adjusted p value", adjusted_messages,
                        ignore.case = TRUE)))

  f <- function(m) list(x = rnorm(40), y = rnorm(50, m))
  pwr <- twosample_power(f, c(0, 0.5), B = 20, maxProcessor = 1,
                         SuppressMessages = TRUE)
  expect_valid_power_result(pwr)
  expect_equal(nrow(pwr), 2L)

  pwr_custom <- twosample_power(f, c(0, 0.5), TS = myTSc1,
                                B = 20, maxProcessor = 1,
                                SuppressMessages = TRUE)
  expect_valid_power_result(pwr_custom)

  pwr_hybrid <- twosample_power(hybrid_cont, c(0, 0.5), rnull = rnull_cont,
                                TS = myTSc2, TSextra = TSextra_cont,
                                B = 20, maxProcessor = 1,
                                SuppressMessages = TRUE)
  expect_valid_power_result(pwr_hybrid)
  
 
})
