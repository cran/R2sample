#' Test computation time
#'
#' Estimate the time required for one calculation of the test statistic.
#'
#' @param dta data set
#' @param TS test-statistic routine
#' @param typeTS format of TS
#' @param TSextra additional information passed to TS
#' @return Mean computation time in seconds
#' @export
timecheck <- function(dta, TS, typeTS, TSextra) {
  f <- function() calcTS(dta, TS, typeTS, TSextra)
  a <- microbenchmark::microbenchmark(f(), unit="seconds", times=10)
  as.numeric(summary(a)[["mean"]])
}
