#' This function runs the chi-square test for continuous or discrete data
#' @param  rxy  a function to generate data
#' @param  alpha =0.05 type I error probability of test
#' @param  B =1000 number of simulation runs
#' @param  xparam vector of parameter values
#' @param  yparam vector of parameter values
#' @param  nbins =c(50, 10) number of desired bins
#' @param  minexpcount =5 smallest number of counts required in each bin
#' @param  typeTS  type of problem, continuous/discrete, with/without weights
#' @return A matrix of power values

chi_power = function(rxy, alpha=0.05, B=1000, xparam, yparam, 
                     nbins=c(50, 10), minexpcount=5, typeTS) {
   dta = rxy(xparam[1], yparam[1])
   if(typeTS<5) Continuous = TRUE
   else Continuous = FALSE
   if(Continuous) {
       pwr=matrix(0, length(xparam), 4)
       colnames(pwr) = c("ES large", "ES small", "EP large", "EP small")
   }     
   else {
       pwr=matrix(0, length(xparam), 2)
       colnames(pwr) = c("large", "small")
   }
   for(i in 1:length(xparam)) {
     for(j in 1:B) {
         dta=rxy(xparam[i], yparam[i])
         tmp = chi_test(dta, nbins, minexpcount,typeTS)$p.values 
         pwr[i, ] = pwr[i, ] + ifelse(tmp<alpha, 1, 0)
     }
     pwr[i, ] = pwr[i, ]/B
   }
   pwr
}
