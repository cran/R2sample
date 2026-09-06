#' New test for continuous data
#' 
#' This function runs a test for continuous data
#' @param x first data set
#' @param y second data set
#' @param TSextra list with extra info
#' @return a number
#' @export
myTS_cont = function(x, y, TSextra) {
  z=c(x, y)
  nbins=TSextra$nbins
  bins=quantile(z, probs = (0:nbins)/nbins)
  Ox=hist(x, bins, plot = FALSE)$counts
  Oy=hist(y, bins, plot = FALSE)$counts
  tmp=sqrt(sum(Ox)/sum(Oy))
  chi = sum((Ox/tmp-Oy*tmp)^2/(Ox+Oy))
  p = 1-stats::pchisq(chi, length(Ox)-1)
  if(TSextra$statistic) {out=chi}
  else out=p
  names(out)="Chi"
  out
}
#' This function runs a test for discrete data
#' @param x counts for first data set
#' @param y counts for second data set
#' @param vals values of discrete data
#' @param TSextra list with extra info
#' @return a number
#' @export
myTS_disc = function(x, y, vals, TSextra) {
  tmp=sqrt(sum(x)/sum(y))
  chi = sum((x/tmp-y*tmp)^2/(x+y))
  p = 1-stats::pchisq(chi, length(x)-1)
  if(TSextra$statistic) {out=chi}
  else out=p
  names(out)="Chi"
  out
} 