#'  twosample_test
#'
#' This function runs a number of two sample tests using Rcpp and parallel computing.
#' @param  x  a vector of numbers or of counts  
#' @param  y a vector of numbers or of counts  
#' @param  vals a vector of numbers, the values of a discrete vector random variable. If it is missing, continuous data is assumed.
#' @param  B number of simulation runs for permutation test
#' @param  nbins Number of bins for chi square tests.
#' @param  maxProcessor Maximum number of cores to use. If maxProcessor=1 no parallel computing is used.
#' @param  discretize  Should continuous data be binned?
#' @param  doMethod Which methods should be included? If missing default methods are used.
#' @return A list of two numeric vectors, the test statistics and the p values. 
#' @export 
#' @examples
#'  twosample_test(rnorm(1000), rt(1000, 4), B=1000, maxProcessor = 1)
#'  vals=1:5
#'  x=table(sample(vals, size=100, replace=TRUE))
#'  y=table(sample(vals, size=100, replace=TRUE, prob=c(1,1,2,1,1)))
#'  twosample_test(x, y, vals, maxProcessor = 1)

twosample_test=function(x, y, vals, B=5000, nbins=c(100,10), maxProcessor=10, 
                        discretize=FALSE, doMethod) {
    default.methods = list(cont=c("chi small", "ZA", "ZK", "Wassp1"), 
                           disc=c("chi small", "ZA", "Kuiper", "Wassp1"))                             
    m=parallel::detectCores()  
    if(m==1) maxProcessor=1
    allmethods=c("chi large", "chi small", "t test", "KS", "Kuiper", "CvM", "AD", "LR", 
                   "ZA", "ZK", "ZC", "Wassp1")               
    if(missing(doMethod)) {
       if(missing(vals) && !discretize) doMethod = default.methods$cont
       else doMethod = default.methods$disc
    }
    else if(doMethod[1]=="all") doMethod=allmethods
# Check whether continuous data has many ties
    if(missing(vals)) {
       if(length(unique(c(x,y)))*3<length(c(x,y)))
         message("There appear to be many ties. If data is actually discrete, run discrete 
         version of twosample_test. For details see Help\n")
    }
# if data is continuous and discretize=TRUE, bin x and y data using nbins[1] 
# equal-probability bins and rerun twosample_test.
    if(missing(vals) & discretize) {
         
         bins=stats::quantile(c(x,y),0:nbins[1]/nbins[1])
         newx=graphics::hist(x, bins, plot=FALSE)$counts
         newy=graphics::hist(y, bins, plot=FALSE)$counts
         vals=(bins[-1]+bins[-(nbins[1]+1)])/2
         return(twosample_test(newx, newy, vals=vals, B=B, nbins=nbins, 
             maxProcessor=maxProcessor, doMethod=doMethod))
    }  
# if data is discrete, check that all vectors have equal lengths, that there are no empty classes
# and ZK method is not run
    if(!missing(vals)) {
       if(length(x)!=length(vals) | length(y)!=length(vals)) {
           message("x, y and vals have to have equal lengths. (x or y can be 0)!\n")
           return(NULL)       
       }
       if(min(x+y)==0) {
           message("Bins with 0 counts of both x and y are not allowed!\n")
           return(NULL)
       }
       doMethod=doMethod[doMethod!="ZK"]
    }      

# is data is continuous, do a quick timing run and if full run takes more than 
# a few seconds print an estimate of the run time. 
    if(missing(vals)) {
        tm=microbenchmark::microbenchmark(
          perm_test_cpp(x=x, y=y, B=100, nbins=nbins, doMethod=doMethod),
          times=1)$time/1e9/100*B/max(1,m-1)
        unit="seconds"  
        if(tm>10) {
          if(tm>120) {tm=round(tm/60);unit="minutes"}  
          else tm=round(tm,-1)
          message("Estimated computation time :", tm, " - ", 
                  tm+ifelse(unit=="seconds",10,1), " ", unit, "\n") 
       }                             
    }
    
# if either only one core is present, B=0 or maxProcessor==1, run perm_test_cpp. 
# If B=0, return test statistics only.    
    if(B==0 | maxProcessor==1) {
       if(missing(vals)) out=perm_test_cpp(x=x, y=y, B=B, nbins=nbins, doMethod=doMethod)
       else out=perm_test_cpp(x=x, y=y, vals=vals, B=B, nbins=nbins, doMethod=doMethod) 
       names(out$statistics)=allmethods  
       if(B==0) return(out$statistics[doMethod]) 
       names(out$p.values)=allmethods  
       out$statistics=out$statistics[doMethod]
       out$p.values=out$p.values[doMethod]
       return(out)
    }
# run perm_test_cpp in parallel. Use one less core than is present, at most maxProcessor.    
    m=min(maxProcessor+1, m)-1
    cl=parallel::makeCluster(m)
    if(missing(vals)) z=parallel::clusterCall(cl, perm_test_cpp, x=x, y=y, 
                                  B=B/m, nbins=nbins, doMethod=doMethod)
    else z=parallel::clusterCall(cl, perm_test_cpp, x=x, y=y, vals=vals, 
                                  B=B/m, nbins=nbins, doMethod=doMethod)
    parallel::stopCluster(cl)
      
# average p values over cores    
    statistics=z[[1]]$statistics
    names(statistics)=allmethods  
    p=z[[1]]$p.values
    for(i in 2:m) p=p+z[[i]]$p.values
    names(p)=allmethods    
    out=list(statistics=statistics[doMethod], p.values=round(p[doMethod]/m, 4))
    out
}
