#' This function runs a number of two sample tests using Rcpp and parallel computing.
#' @param  x  a vector of numbers if data is continuous or of counts  if data is discrete or a list with the data
#' @param  y a vector of numbers if data is continuous or of counts  if data is discrete.
#' @param  vals =NA, a vector of numbers, the values of a discrete random variable. NA if data is continuous data.
#' @param  TS routine to calculate test statistics for non-chi-square tests
#' @param  TSextra additional info passed to TS, if necessary
#' @param  wx A numeric vector of weights of x.
#' @param  wy A numeric vector of weights of y.
#' @param  B =5000, number of simulation runs for permutation test
#' @param  nbins =c(50,10), number of bins for chi square tests.
#' @param  minexpcount =5, minimum required expected counts for chi-square tests.
#' @param  maxProcessor maximum number of cores to use. If missing (the default) no parallel processing is used.
#' @param  UseLargeSample should p values be found via large sample theory if n,m>10000?
#' @param  samplingmethod ="independence" or "MCMC" for discrete data
#' @param  rnull a function that generates data from a model, possibly with parameter estimation.
#' @param  doMethods ="all" Which methods should be included? If missing all methods are used.
#' @param  SuppressMessages = FALSE print informative messages?
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
                        minexpcount=5, maxProcessor,  UseLargeSample, 
                        samplingmethod="independence", rnull, 
                        SuppressMessages = FALSE, doMethods="all") {
    default.methods = list(cont=c("EP small", "ZA", "ZK", "Wassp1"), 
                           disc=c("small", "ZA", "Kuiper", "Wassp1"))
    if(is.list(x)) {
        y=x$y
        if("vals"%in%names(x)) vals=x$vals
        if("wx"%in%names(x)) wx=x$wx
        if("wy"%in%names(x)) wy=x$wy
        x=x$x
    }
    if(missing(maxProcessor)) {
      if(min(c(length(x), length(y)))<1e5) {
        maxProcessor=1
      }   
      else {
        maxProcessor=max(1, parallel::detectCores(logical = FALSE)-1)   
        if(!SuppressMessages) message(paste("Using",maxProcessor," cores"))  
      }  
    }
    samplingmethod=ifelse(samplingmethod=="independence", 1, 2)
    if(missing(TSextra)) TSextra = list(samplingmethod=samplingmethod)
    else TSextra=c(TSextra, samplingmethod=samplingmethod)
    if(!missing(rnull)) TSextra=c(TSextra, rnull=rnull)
    Continuous=ifelse(any(is.na(vals)), TRUE, FALSE)
    CustomTS=ifelse(missing(TS), FALSE, TRUE)
    WithWeights=ifelse(all(c(wx,wy)==1), FALSE, TRUE)
    if(missing(UseLargeSample))
       UseLargeSample=ifelse(min(length(x), length(y))<1e4, FALSE, TRUE)  
    outchi = list(statistics=NULL, p.values=NULL)
# what methods are to be run?  
    if(!CustomTS) { # do included methods
      if(Continuous) {
          if(!WithWeights) { # No weights
             typeTS = 1
             TS=TS_cont
             dta=list(x=x, y=y)
             tmp=TS(x,y)
             outchi = chi_test(dta, nbins, minexpcount, typeTS=1)
          }
          else { # weighted data
            typeTS = 2
            TS=TSw_cont
            dta=list(x=x, y=y, wx=wx, wy=wy)
            tmp=TS(x, y, wx, wy)
            outchi = chi_test(dta, nbins, minexpcount, typeTS=2)
            doMethods=c("KS", "Kuiper", "CvM", "AD", names(outchi[[1]]))
            if(min(wx,wy)<0.01 | max(wx,wy)>5)
               if(!SuppressMessages) message("Some of the weights are either exceptionally small or large. In either case the tests can be unreliable!")
          }
          
      }
      else {
          typeTS=5
          TS=TS_disc
          dta=list(x=x, y=y, vals=vals)
          TSextra=c(TSextra, list(adw=weights(dta)))
          tmp=TS(x, y, vals, TSextra$adw)
          outchi = chi_test(dta, nbins, minexpcount, typeTS=5)        
      }  
      chimethods = names(outchi[[1]])
    }  
    else { # do user-supplied tests
      if(substr(deparse(TS)[2], 1, 5)==".Call") {
         if(maxProcessor>1) {
            if(!SuppressMessages) message("Parallel Programming is not possible if custom TS is written in C++. Switching to single processor")  
            maxProcessor=1
         }  
      }
      if(Continuous) {
        if(length(formals(TS))==2) {
          typeTS=3
          dta=list(x=x, y=y)
          tmp=TS(x, y)
        }
        else {
          typeTS=4
          dta=list(x=x, y=y)
          tmp=TS(x, y, TSextra)
        }
      }
      else {
        if(length(formals(TS))==3) {
          typeTS=6
          dta=list(x=x, y=y, vals=vals)
          tmp=TS(x, y, vals)
        }  
        else {
          typeTS=7
          dta=list(x=x, y=y, vals=vals)
          tmp=TS(x, y, vals, TSextra)
        }
      }
      if(is.null(names(tmp))) {
         message("output of TS routine has to be a named vector!")
        return(NULL)
      }  
    }

#  set number of processors for parallel programming
    m=parallel::detectCores(logical = FALSE) 
    if(maxProcessor>1) maxProcessor=min(maxProcessor, m-1)
    
# Check whether continuous data has many ties
    if(Continuous) {
       if(length(unique(c(x,y)))*3<length(c(x,y)))
         if(!SuppressMessages) message("There appear to be many ties. If data is actually discrete, run discrete 
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
       if(!missing(rnull)) {
          message("Large sample methods are not implemented for GOF-Twosample hybrid problem")
          return(NULL)
       }
       statistics = tmp[c("KS", "Kuiper", "CvM", "AD")]
       p.values=asymptotic_pvalues(statistics, length(x), length(y))  
       s = c(statistics, outchi$statistics)
       p = c(p.values, outchi$p.values)
       return(list(statistics=s, p.values=p))
    }  
    
# if either only one core is present, B=0 or maxProcessor=1, run testC.
    if(B==0 | maxProcessor==1) {
       outTS = testC(dta, TS, typeTS, TSextra,B=B)
    }
    else {
# run testC in parallel. Use one less core than is present, at most maxProcessor.    
      cl=parallel::makeCluster(maxProcessor)
      z=parallel::clusterCall(cl, testC, dta=dta, 
                      TS=TS, typeTS=typeTS, TSextra=TSextra, B= round(B/m))
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
