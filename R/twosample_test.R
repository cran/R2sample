#' This function runs a number of two sample tests using Rcpp and parallel computing.
#' @param  x  a vector of numbers if data is continuous or of counts  if data is discrete.
#' @param  y a vector of numbers if data is continuous or of counts  if data is discrete.
#' @param  vals =NA, a vector of numbers, the values of a discrete random variable. NA if data is continuous data.
#' @param  TS routine to calculate test statistics for non-chi-square tests
#' @param  TSextra additional info passed to TS, if necessary
#' @param  wx A numeric vector of weights of x.
#' @param  wy A numeric vector of weights of y.
#' @param  B =5000, number of simulation runs for permutation test
#' @param  nbins =c(50,10), number of bins for chi square tests.
#' @param  maxProcessor maximum number of cores to use. If missing (the default) no parallel processing is used.
#' @param  UseLargeSample should p values be found via large sample theory if n,m>10000?
#' @param  samplingmethod ="independence" or "MCMC" for discrete data
#' @param  doMethods ="all" Which methods should be included? If missing all methods are used.
#' @return A list of two numeric vectors, the test statistics and the p values. 
#' @export 
#' @examples
#'  R2sample::twosample_test(rnorm(1000), rt(1000, 4), B=1000)
#'  myTS=function(x,y) {z=c(mean(x)-mean(y),sd(x)-sd(y));names(z)=c("M","S");z}
#'  R2sample::twosample_test(rnorm(1000), rt(1000, 4), TS=myTS, B=1000)
#'  vals=1:5
#'  x=table(sample(vals, size=100, replace=TRUE))
#'  y=table(sample(vals, size=100, replace=TRUE, prob=c(1,1,3,1,1)))
#'  R2sample::twosample_test(x, y, vals)

twosample_test=function(x, y, vals=NA, TS, TSextra, wx=rep(1, length(x)),
                        wy=rep(1, length(y)), B=5000, nbins=c(50,10),
                        maxProcessor,  UseLargeSample, 
                        samplingmethod="independence", doMethods="all") {
    default.methods = list(cont=c("EP small", "ZA", "ZK", "Wassp1"), 
                           disc=c("small", "ZA", "Kuiper", "Wassp1"))
    samplingmethod=ifelse(samplingmethod=="independence", 1, 2)
    if(missing(maxProcessor)) maxProcessor=1
    if(missing(TSextra)) TSextra = list(aaaa=0)
    if(any(is.na(vals))) Continuous=TRUE
    else Continuous=FALSE
    if(missing(UseLargeSample))
      UseLargeSample=ifelse(min(length(x), length(y))<1e4, FALSE, TRUE)  
    outchi = list(statistics=NULL, p.values=NULL)
# what methods are to be run?  
    CustomTS=TRUE
    if(missing(TS)) { # do included methods
      CustomTS=FALSE
      if(Continuous) {
          if(all( abs(c(wx,wy)-round(c(wx,wy)))<1e-10 )) { # No weights
             outchi = chi_test(list(x=x, y=y), nbins, 5, typeTS=1)
             typeTS = 1
             TS=TS_cont
             tmp=TS(x,y)
          }
          else { # weighted data
            outchi = chi_test(list(x=x, y=y, wx=wx, wy=wy), nbins, 5, typeTS=2)
            typeTS = 2
            TS=TSw_cont
            tmp=TS(x, y, wx, wy)
            doMethods=c("KS", "Kuiper", "CvM", "AD", names(outchi[[1]]))
            if(min(wx,wy)<0.01 | max(wx,wy)>5)
               message("Some of the weights are either exceptionally small or large. In either case the tests can be unreliable!")
          }
          
      }
      else {
          if(all( abs(c(x,y)-round(c(x,y)))<1e-10 )) { # no weights
            typeTS=5
            TS=TS_disc
            tmp=TS(x, y, vals, rep(1, length(x)))
            dta=list(x=x, y=y, vals=vals)
            outchi = chi_test(dta, nbins, 5, typeTS=5)        
          }
          else {
            typeTS=6
            TS=TSw_disc
            tmp=TS(x, y)
            if("aaaa" %in% names(TSextra)) {
               message("Weighted discrete data, but no TSextra with sum of squared weights, so we can't do the chi square tests")
            }
            else {
               dta=list(x=x, y=y, vals=vals, TSextra=TSextra)
               outchi = chi_test(dta, nbins, 5, typeTS=6)        
               
            }  
            doMethods=c("KS", "Kuiper", "CvM", "AD")
          }
      }  
      chimethods = names(outchi[[1]])
    }  
    else { # do user-supplied tests
      if(substr(deparse(TS)[2], 1, 5)==".Call") {
         if(maxProcessor>1) {
            message("Parallel Programming is not possible if custom TS is written in C++. Switching to single processor")  
            maxProcessor=1
         }  
      }
      if(Continuous) {
        if(length(formals(TS))==2) {
          typeTS=3
          tmp=TS(x, y)
        }
        else {
          typeTS=4
          tmp=TS(x, y, TSextra)
        }
      }
      else {
        if(length(formals(TS))==3) {
          typeTS=7
          TSextra = list("aaaa"=0)
          tmp=TS(x, y, vals)
        }  
        else {
          typeTS=8
          tmp=TS(x, y, vals, TSextra)
        }
      }
      if(is.null(names(tmp))) {
         message("output of TS routine has to be a named vector!")
        return(NULL)
      }  
    }

#  set number of processors for parallel programming
    m=parallel::detectCores()  
    if(maxProcessor>1) maxProcessor=min(maxProcessor, m-1)
    
# Check whether continuous data has many ties
    if(Continuous) {
       if(length(unique(c(x,y)))*3<length(c(x,y)))
         message("There appear to be many ties. If data is actually discrete, run discrete 
         version of twosample_test. For details see Help\n")
    }
  
# if data is discrete, check that all vectors have equal lengths, that there are no empty classes
# and ZK method is not run
    if(!Continuous) {
       if(length(x)!=length(vals) | length(y)!=length(vals)) {
           message("x, y and vals have to have equal lengths. (x or y can be 0)!\n")
           return(NULL)       
       }
       if(min(x+y)==0) {
           message("Bins with 0 counts of both x and y are not allowed!\n")
           return(NULL)
       }
    }      

# if data is continuous and min(n,m)>10000, use asymptotic formulas for p value
    if(UseLargeSample) {
       statistics = tmp[c("KS", "Kuiper", "CvM", "AD")]
       p.values=asymptotic_pvalues(statistics, length(x), length(y))  
       s = c(statistics, outchi$statistics)
       p = c(p.values, outchi$p.values)
       return(list(statistics=s, p.values=p))
    }  
 
# if either only one core is present, B=0 or maxProcessor=1, run perm_test. 
    if(B==0 | maxProcessor==1) {
       if(Continuous) outTS = perm_test_cont(x, y, TS, typeTS, TSextra, wx=wx, wy=wy,B=B)
       else outTS = perm_test_disc(x, y, vals, TS, typeTS, TSextra, 
                                   samplingmethod=samplingmethod, B=B)
    }
    else {
# run perm_test_cpp in parallel. Use one less core than is present, at most maxProcessor.    
      m=min(maxProcessor+1, m)-1
      cl=parallel::makeCluster(m)
      if(Continuous) z=parallel::clusterCall(cl, perm_test_cont, x=x, y=y, 
                            TS=TS, typeTS=typeTS, TSextra=TSextra, wx=wx, wy=wy,B=B/m)
      else z=parallel::clusterCall(cl, perm_test_disc, x=x, y=y, vals=vals, 
                            TS=TS, typeTS=typeTS, TSextra=TSextra,
                            samplingmethod=samplingmethod, B=B/m)
      parallel::stopCluster(cl)
# average p values over cores    
      p=z[[1]]$p.values
      for(i in 2:m) p=p+z[[i]]$p.values
      p = round(p/m, 4)  
      outTS = list(statistics=z[[1]]$statistics, p.values=p)
    }
    if(CustomTS) return(signif.digits(outTS))
    s = c(outTS$statistics, outchi$statistics)
    p = c(outTS$p.values, outchi$p.values)
    out = list(statistics=s, p.values=p)
    out=signif.digits(out)
    if(doMethods[1]=="all") return(out)
    if(typeTS%in%c(2,6)) return(out)
    if(doMethods[1]=="default") {
      if(Continuous) doMethods = default.methods$cont
      else doMethods = default.methods$disc
    }
    out = list(statistics=s[doMethods], p.values=p[doMethods])
    out
}
