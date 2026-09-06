test_that("use_parallel chooses the appropriate number of processors", {
  
  # Short computation: parallel overhead is not worthwhile
  expect_equal(
    use_parallel(0.001, 2, 10, 2),
    1L
  )
  
  # Long computation: use requested processors
  expect_equal(
    use_parallel(1, 2, 10, 2),
    2L
  )
  
  # Respect explicit request for one processor
  expect_equal(
    use_parallel(1, 2, 10, 1),
    1L
  )
})