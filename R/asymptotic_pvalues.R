#' This function finds the p values of several tests based on large sample theory
#'
#' KS: Kolmogorov distribution. Kuiper: asymptotic distribution with the
#' 1/sqrt(n) correction term (Kuiper 1960, Stephens 1965), using the effective
#' sample size nm/(n+m). CvM: limiting distribution of the two-sample statistic
#' (Anderson and Darling 1952 series). AD: Marsaglia and Marsaglia (2004)
#' approximation of the limiting distribution.
#' @param  x a vector of test statistics
#' @param  n size of sample 1
#' @param  m size of sample 2
#' @keywords internal
#' @return A vector of p values.
#' @noRd

asymptotic_pvalues = function(x, n, m) {
  neff = n*m/(n+m)
  KSp = function(x) {
    z = x*sqrt(neff)
    if(z<=0) return(1)
    k = 1:100
    min(1, max(0, 2*sum((-1)^(k-1)*exp(-2*k^2*z^2))))
  }
  Kuiperp = function(kp) {
    z = kp*sqrt(neff)
    if(z<0.4) return(1)
    j = 1:100
    s1 = sum(2*(4*j^2*z^2-1)*exp(-2*j^2*z^2))
    s2 = sum(j^2*(4*j^2*z^2-3)*exp(-2*j^2*z^2))
    min(1, max(0, s1 - 8*z/(3*sqrt(neff))*s2))
  }
  CvMp = function(x) {
    if(x<=0) return(1)
    j = 0:50
    u = (4*j+1)^2/(16*x)
    cf = exp(lgamma(j+0.5)-lgamma(0.5)-lgamma(j+1))
    cdf = sum(cf*sqrt(4*j+1)*besselK(u, 0.25, expon.scaled=TRUE)*exp(-2*u))/(pi*sqrt(x))
    min(1, max(0, 1-cdf))
  }
  ADp <- function(x) {
    if(x<=0) return(1)
    if(x < 2)
      x=exp(-1.2337141/x)/sqrt(x)*(2.00012+(.247105-
              (.0649821-(.0347962-(.011672-.00168691*x)*x)*x)*x)*x)
    else
      x=exp(-exp(1.0776-(2.30695-(.43424-
                (.082433-(.008056 -.0003146*x)*x)*x)*x)*x))
    min(1, max(0, 1-x))
  }
  out=c(KSp(x[1]), Kuiperp(x[2]), CvMp(x[3]), ADp(x[4]))
  names(out) = c("KS", "Kuiper", "CvM", "AD")
  out
}
