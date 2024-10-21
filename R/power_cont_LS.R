#' Find the power of various discrete tests via permutation.
#' @param rxy a function that generates x and y data.
#' @param alpha A numeric constant
#' @param B Number of simulation runs.
#' @param xparam  arguments for r1.
#' @param yparam  arguments for r2.
#' @keywords internal
#' @return A numeric matrix of powers

power_cont_LS = function(rxy, alpha=0.05, B=1000, xparam=0.0, yparam=0.0) { 
# Find out how many tests are to be done, sample sizes etc. 
  dta=rxy(xparam[1],yparam[1])
  x = dta["x"]
  y = dta["y"]
  if(length(dta)==4) {
    wx = dta["wx"]
    wy = dta["wy"]
  }  
  tmp=twosample_test(dta$x, dta$y, UseLargeSample = TRUE)$p.value
  pwr=matrix(0, length(xparam), length(tmp))
#  l loop over values in xparam, yparam 
  for(m in 1:B) {
#  m loop over simulation runs 
    for(l in 1:length(xparam)) {
        dta=rxy(xparam[l],yparam[l])
        if(length(dta)==2)
          tmp=twosample_test(dta$x, dta$y, UseLargeSample = TRUE)
        else
          tmp=twosample_test(dta$x, dta$y,
                             wx=dta$wx, wy=dta$wy, UseLargeSample = TRUE)
        pwr[l, ]=pwr[l, ]+ifelse(tmp$p.values<alpha, 1, 0)
    }  # end of loop of power simulation
  } 
  pwr/B
}

