#' test function
#' @param  dta data set
#' @param  TS test statistics
#' @param  typeTS format of TS
#' @param  TSextra additional info TS
#' @return Mean computation time
#' @export
timecheck=function(dta, TS, typeTS, TSextra) {
  if(typeTS==1) f=function() TS(dta$x, dta$y)
  if(typeTS==2) f=function() TS(dta$x, dta$y, dta$wx, dta$wy)
  if(typeTS==3) f=function() TS(dta$x, dta$y)
  if(typeTS==4) f=function() TS(dta$x, dta$y, TSextra)
  a=microbenchmark::microbenchmark(f(), 
              unit="seconds", times=10)
  as.numeric(summary(a)["mean"])
}
