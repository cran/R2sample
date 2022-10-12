#' twosample_power
#'
#' Find the power of various two sample tests using Rcpp and parallel computing.
#' @param  rxy  function to generate a list with data sets x, y and (optional) vals
#' @param  alpha =0.05, level of hypothesis test 
#' @param  B =1000, number of simulation runs for permutation test and power.
#' @param  avals parameter for function that generates x and y.
#' @param  bvals parameter for function that generates x and y.
#' @param  nbins =c(100,10), number of bins for chi large and chi small.
#' @param  maxProcessor =10, maximum number of cores to use. If maxProcessor=1 no parallel computing is used.
#' @param  doMethod ="all", which methods should be included?
#' @return A numeric vector of power values.
#' @export 
#' @examples
#'  rxy=function(a,b) list(x=rnorm(50), y=rnorm(50,0.5))
#'  twosample_power(rxy, B=100, maxProcessor = 1)
#'  rxy=function(a,b) list(x=table(sample(1:10, size=1000, replace=TRUE)), 
#'        y=table(sample(1:10, size=1200, replace=TRUE)), vals=1:10)
#'  twosample_power(rxy, B=100, maxProcessor = 1)

twosample_power=function(rxy, alpha=0.05, B=1000, avals=0, bvals=0, 
                       nbins=c(100,10), maxProcessor=10, doMethod="all") {

# which methods should be run?   
    if(doMethod[1]=="all") 
        doMethod=c("chi large", "chi small", "t test", "KS", "Kuiper", 
                    "CvM", "AD", "LR",  "ZA", "ZK", "ZC", "Wassp1")
# If data is discrete, make sure ZK method is eliminated
    if(length(rxy(avals[1],bvals[1]))==3)  #Data is discrete
        doMethod=doMethod[doMethod!="ZK"]  
       
# check that avals and bvals have the same length. 
# If they do, matrix of powers is returned without row names.
# If one of them is a scalar, make it the same length as the other and use those
# values as row names
    if(length(avals)!=length(bvals)) {
       if(min(c(length(avals),length(bvals)))>1) {
          message("lengths of parameter vectors not compatible!\n")
          return(NULL)
       }
       if(length(avals)==1) {
           avals=rep(avals, length(bvals))
           rnames=bvals   
       }    
       else {
           bvals=rep(bvals, length(avals))
           rnames=avals
       }    
    }
    else rnames=1:length(avals)
                        
    m=parallel::detectCores()
    if(m==1 | maxProcessor==1) {
       if(m==1) message("No multiple cores found, using one core\n")    
       tmp=power_cpp(rxy=rxy, alpha=alpha, B=B, 
                  xparam=avals, yparam=bvals, nbins=nbins, doMethod=doMethod)[,doMethod]
       if(length(avals)>1 & length(doMethod)>1) rownames(tmp)=rnames           
       return(tmp)
    }  
# Use one less processor than is present, or at most maxProcessor    
    m=min(m, maxProcessor+1)-1
    cl <- parallel::makeCluster(m)
    z=parallel::clusterCall(cl, power_cpp, rxy=rxy, alpha=alpha, B=round(B/m), 
                     xparam=avals, yparam=bvals, nbins=nbins, doMethod=doMethod)
    parallel::stopCluster(cl)  
# Average power of cores    
    out=0*z[[1]]
    for(i in 1:m) out=out+z[[i]]
    rownames(out)=rnames           
    out[,doMethod]/m
}
