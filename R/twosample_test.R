#' Tests for the univariate two-sample problem
#' 
#' This function runs a number of two sample tests using Rcpp and parallel computing.
#' 
#' For details consult vignette("R2sample","R2sample")
#' 
#' @param  x  a vector of numbers if data is continuous or of counts  if data is discrete or a list with the data
#' @param  y a vector of numbers if data is continuous or of counts  if data is discrete.
#' @param  vals =NA, a vector of numbers, the values of a discrete random variable. NA if data is continuous data.
#' @param  TS routine to calculate test statistics for non-chi-square tests
#' @param  TSextra additional info passed to TS, if necessary
#' @param  wx A numeric vector of weights of x.
#' @param  wy A numeric vector of weights of y.
#' @param  B =5000, number of simulation runs for permutation test
#' @param  seed =NULL, optional seed for reproducibility
#' @param  nbins =c(50,10), number of bins for chi square tests.
#' @param  minexpcount =5, minimum required expected counts for chi-square tests.
#' @param  maxProcessor maximum number of cores to use. If missing (the default) no parallel processing is used.
#' @param  UseLargeSample should p values be found via large sample theory if n,m>10000?
#' @param  samplingmethod ="Binomial" or "independence" for discrete data
#' @param  rnull a function that generates data from a model, possibly with parameter estimation.
#' @param  doMethods ="all" a vector of codes for the methods to include. If "all", all methods are used.
#'         For continuous data without weights and for discrete data this includes the HHG test 
#'         (Heller, Heller, Kaufman, Brill and Gorfine 2016), see the vignette.
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
                        wy=rep(1, length(y)), B=5000, 
                        seed=NULL, nbins=c(50,10),
                        minexpcount=5, maxProcessor,  UseLargeSample, 
                        samplingmethod="Binomial", rnull, 
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
    set_seed(seed)
    if(!is.numeric(samplingmethod))
        samplingmethod=ifelse(samplingmethod=="independence", 1, 2)
    if(missing(rnull))
      TSextra <- makeTSextra(TSextra, samplingmethod)
    else
      TSextra <- makeTSextra(TSextra, samplingmethod, rnull)
    Continuous=ifelse(any(is.na(vals)), TRUE, FALSE)
    CustomTS=ifelse(missing(TS), FALSE, TRUE)
    WithWeights=ifelse(all(c(wx,wy)==1), FALSE, TRUE)
    if(CustomTS) tsinfo <- maketypeTS(TS, Continuous, WithWeights)
    else tsinfo <- maketypeTS(Continuous=Continuous, WithWeights=WithWeights)
    TS <- tsinfo$TS
    typeTS <- tsinfo$typeTS
    if(!Continuous) {
       if(length(x)!=length(vals) | length(y)!=length(vals)) {
           message("x, y and vals have to have equal lengths. (x or y can be 0)!\n")
           return(NULL)
       }
       if(any(duplicated(vals))) {
           message("vals can not have repeated values!\n")
           return(NULL)
       }
       if(is.unsorted(vals)) { # all methods assume vals are increasing
           I=order(vals)
           x=x[I]
           y=y[I]
           vals=vals[I]
       }
    }
# add the HHG test to the built-in tests (continuous data without weights, discrete data)
    if(!CustomTS && !WithWeights) TS <- add_HHG(TS, x, y, Continuous)
    if(missing(UseLargeSample))
       UseLargeSample=ifelse(Continuous && !WithWeights && missing(rnull) &&
                               min(length(x), length(y))>=1e4, TRUE, FALSE)
    if(UseLargeSample && (!Continuous || WithWeights)) {
       if(!SuppressMessages) message("Large sample methods are only available for continuous data without weights. Using permutation test.")
       UseLargeSample=FALSE
    }
    outchi = list(statistics=NULL, p.values=NULL)
    if(test_methods(doMethods,Continuous, UseLargeSample, WithWeights))
       return(NULL)
# what methods are to be run?  
    if(!CustomTS) {
      if(Continuous) {
        if(WithWeights) dta=list(x=x,y=y,wx=wx,wy=wy)
        else dta=list(x=x,y=y)
        tmp=calcTS(dta, TS, typeTS, TSextra)
        outchi=chi_test(dta, nbins, minexpcount, typeTS=typeTS)
        if(WithWeights) {
          doMethods=c("KS","Kuiper","CvM","AD",names(outchi[[1]]))
          if(min(wx,wy)<0.01 | max(wx,wy)>5)
            if(!SuppressMessages) message("Some of the weights are either exceptionally small or large. In either case the tests can be unreliable!")
        }
      } else {
        dta=list(x=x,y=y,vals=vals)
        TSextra$adw=weights(dta)
        tmp=calcTS(dta, TS, typeTS, TSextra)
        outchi=chi_test(dta, nbins, minexpcount, typeTS=typeTS)
      }
      chimethods=names(outchi[[1]])
    } else {
      if(Continuous) dta=list(x=x,y=y)
      else dta=list(x=x,y=y,vals=vals)
      tmp=calcTS(dta, TS, typeTS, TSextra)
      if(is.null(names(tmp))) {
        message("output of TS routine has to be a named vector!")
        return(NULL)
      }
    }

# Choose number of processors after dta, TS and typeTS are known.
    if(missing(maxProcessor))
      maxProcessor <- makemaxProcessor(
        dta=dta, TS=TS, typeTS=typeTS, TSextra=TSextra, B=B,
        SuppressMessages=SuppressMessages,
        useSingleProcessor=tsinfo$useSingleProcessor
      )
    else
      maxProcessor <- makemaxProcessor(
        maxProcessor=maxProcessor,
        dta=dta, TS=TS, typeTS=typeTS, TSextra=TSextra, B=B,
        SuppressMessages=SuppressMessages,
        useSingleProcessor=tsinfo$useSingleProcessor
      )

# Assure that B is a multiple of maxProcessor.
    B=ceiling(B/maxProcessor)*maxProcessor
    
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
# Epps-Singleton: asymptotic chi-square distribution; HHG: distribution-free,
# p value from the simulated null distribution of its statistic
       es=EppsS(x, y)
       statistics=c(statistics, EppsS=unname(es["statistic"]), HHG=unname(tmp["HHG"]))
       p.values=c(p.values, 
                  EppsS=stats::pchisq(es[["statistic"]], es[["df"]], lower.tail=FALSE),
                  HHG=hhg_pvalue_cont(x, y))
       s = c(statistics, outchi$statistics)
       p = c(p.values, outchi$p.values)
       return(list(statistics=round(s,4), p.values=round(p,4)))
    }  
    
# if either only one core is present, B=0 or maxProcessor=1, run testC.
    if(B==0 | maxProcessor==1) {
       outTS = testC(dta, TS, typeTS, TSextra,B=B)
    }
    else {
# run testC in parallel. Use one less core than is present, at most maxProcessor.    
      cl=parallel::makeCluster(maxProcessor)
      on.exit(parallel::stopCluster(cl), add=TRUE)
      if(!is.null(seed))
        parallel::clusterSetRNGStream(cl, iseed=seed)
      z=parallel::clusterCall(cl, testC, dta=dta, 
                      TS=TS, typeTS=typeTS, TSextra=TSextra, 
                      B= round(B/maxProcessor))
# average p values over cores    
      p=z[[1]]$p.values
      for(i in 2:length(z)) p=p+z[[i]]$p.values
      p = round(p/length(z), 4)  
      outTS = list(statistics=z[[1]]$statistics, p.values=p)
    }
    if(CustomTS) return(signif_digits(outTS))
    s = c(outTS$statistics, outchi$statistics)
    p = c(outTS$p.values, outchi$p.values)
    out = list(statistics=s, p.values=p)
    out=signif_digits(out)
    if(doMethods[1]=="all") return(out)
    if(typeTS%in%c(2,6)) return(out)
    if(doMethods[1]=="default") {
      if(Continuous) doMethods = default.methods$cont
      else doMethods = default.methods$disc
    }
    out = list(statistics=s[doMethods], p.values=p[doMethods])
    out
}
