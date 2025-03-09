#' This function estimates the power of test routines that calculate p value(s)
#' @param TS routine to calculate test statistics.
#' @param f routine that generates data.
#' @param param_alt values of parameter under the alternative hypothesis. 
#' @param TSextra list passed to TS.
#' @param alpha =0.05  type I error.
#' @param B = 1000 number of simulation runs to estimate the power.
#' @return A matrix of power values

power_newtest <- function(TS, f, param_alt, TSextra, alpha=0.05, B=1000) {
     Continuous=ifelse("vals" %in% names(formals(TS)), FALSE, TRUE)  
     dta=f(param_alt[1])
     if(Continuous) {
        if("aaa"%in% names(TSextra)) tmp=TS(dta$x, dta$y)
        else tmp=TS(dta$x, dta$y, TSextra)
     }
     else {
        if("aaa"%in% names(TSextra)) tmp=TS(dta$x, dta$y, dta$vals)
        else tmp=TS(dta$x, dta$y, dta$vals, TSextra)
     } 
     pwr=matrix(0, length(param_alt), length(TS))
     rownames(pwr)=param_alt
     colnames(pwr)=names(tmp)
     for(j in 1:length(param_alt)) {
        for(i in 1:B) {
          dta=f(param_alt[j])
          if(Continuous) {
             if("aaa"%in% names(TSextra)) tmp=TS(dta$x, dta$y)
             else tmp=TS(dta$x, dta$y, TSextra)
          }
          else {
            if("aaa"%in% names(TSextra)) tmp=TS(dta$x, dta$y, dta$vals)
            else tmp=TS(dta$x, dta$y, dta$vals, TSextra)
          }
          pwr[j, ]=pwr[j, ]+ifelse(tmp<alpha,1,0) 
        }
     }
     pwr/B
}
