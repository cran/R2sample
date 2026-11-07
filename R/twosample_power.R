#' Power estimation for two-sample methods
#' 
#' Find the power of various two sample tests using Rcpp and parallel computing.
#' 
#' For details consult vignette("R2sample","R2sample")
#' 
#' This routine runs a number of different two-sample tests for univariate data,
#' either discrete or continuous. The user can also provide their own test method.
#'  
#' @param  f  function to generate a list with data sets x, y and (optional) vals, weights
#' @param  ... additional arguments passed to f, up to 2
#' @param  TS routine to calculate test statistics for non-chi-square tests
#' @param  TSextra additional info passed to TS, if necessary
#' @param  With.p.value =FALSE does user supplied routine return p values?
#' @param  alpha =0.05, the level of the hypothesis test 
#' @param  B =1000, number of simulation runs.
#' @param  seed =NULL, optional seed for reproducibility
#' @param  nbins =c(50,10), number of bins for chi large and chi small.
#' @param  minexpcount =5 minimum required count for chi square tests
#' @param  UseLargeSample should p values be found via large sample theory if n,m>10000?
#' @param  samplingmethod ="Binomial" or independence in discrete data case
#' @param  rnull a function that generates data from a model, possibly with parameter estimation.
#' @param  SuppressMessages = FALSE print informative messages?
#' @param  maxProcessor maximum number of cores to use. If maxProcessor=1 no parallel computing is used.
#' @return A numeric vector of power values.
#' @export 
#' @examples
#'  # Power of standard normal vs. normal with mean mu.
#'  f1=function(mu) list(x=rnorm(25), y=rnorm(25, mu))
#'  #Power of uniform discrete distribution vs. with different probabilities.
#'  twosample_power(f1, mu=c(0,2), B=100, maxProcessor = 1)
#'  f2=function(n, p) list(x=table(sample(1:5, size=1000, replace=TRUE)), 
#'        y=table(sample(1:5, size=n, replace=TRUE, 
#'        prob=c(1, 1, 1, 1, p))), vals=1:5)
#'  twosample_power(f2, n=c(1000, 2000), p=c(1, 1.5), B=100, maxProcessor = 1)
#'  # Compare power of a new test with those in package:
#'  myTS=function(x,y) {z=c(mean(x)-mean(y),sd(x)-sd(y));names(z)=c("M","S");z}
#'  cbind(twosample_power(f1, mu=c(0,2), TS=myTS,B=100, maxProcessor = 1),
#'        twosample_power(f1, mu=c(0,2), B=100, maxProcessor = 1))
#'  # Power estimation if routine returns a p value, using two parameters for 
#'  #  the alternative
#'  f3=function(mu, n) list(x=rnorm(25), y=rnorm(n, mu))
#'  myTS2=function(x, y) {out=ks.test(x,y)$p.value; names(out)="KSp"; out}      
#'  twosample_power(f3, c(0,1), c(25, 30), TS=myTS2, With.p.value = TRUE,  B=100)
#'  
twosample_power=function(f, ..., TS, TSextra, With.p.value=FALSE, 
            alpha=0.05, B=1000, seed=NULL,
            nbins=c(50,10), minexpcount=5, UseLargeSample, 
            samplingmethod="Binomial", rnull, 
            SuppressMessages=FALSE, maxProcessor) {

  if(!missing(rnull)) {
    if(!missing(UseLargeSample) && isTRUE(UseLargeSample)) {
      message("Large sample methods are not implemented for GOF-Twosample hybrid problem")
      return(NULL)
    }
    UseLargeSample=FALSE
  }
# create function rxy which generates data, with two arguments 
    dots <- list(...)
    if(length(formals(f))==0) { # f has 0 arguments
       if(length(dots)>0)
         message("f has no arguments! Provided arguments will be ignored")
       rxy=function(a=0, b=0) f()
       avals=0
       bvals=0
    }
    if(length(formals(f))==1) { # f has 1 argument
       if(length(dots)!=1) 
         stop("f has to have one argument!", call. = FALSE)
       rxy=function(a=0, b=0) f(a)
       avals=dots[[1]]
       bvals=0
    }
    if(length(formals(f))==2) { # f has 2 arguments
       if(length(dots)!=2)
         stop("f has to have two arguments!", call. = FALSE)
       rxy=function(a=0, b=0) f(a,b)
       avals=dots[[1]]
       bvals=dots[[2]]
    }
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
    }    
    else {
      bvals=rep(bvals, length(avals))
    }    
  }

# generate one data set as an example, do some setup 
  dta = rxy(avals[1], bvals[1])
  if(length(dta)<2 | length(dta)>4) {
    message("routine f has to create a list with 2 (continous data without weights), 3 (discrete data) or 4 (continous data with weights)")
    return(NULL)
  }
  if( (!("x"%in%names(dta))) | (!("y"%in%names(dta))) ) {
    message("routine f has to create a list with elements named x and y, vals (optional, for discrete data) , wx and wy (optional, for continuous weighted data")
    return(NULL)
  }
  x=dta$x
  y=dta$y
  if("vals"%in%names(dta)) vals=dta$vals
  if("wx"%in%names(dta)) wx=dta$wx
  else wx=rep(1, length(x))
  if("wy"%in%names(dta)) wy=dta$wy
  else wy=rep(1, length(y))
  set_seed(seed)
  if(!is.numeric(samplingmethod))
     samplingmethod=ifelse(samplingmethod=="independence", 1, 2)
  Continuous = ifelse("vals" %in% names(dta), FALSE, TRUE)
  if(missing(rnull))
    TSextra <- makeTSextra(TSextra, samplingmethod)
  else
    TSextra <- makeTSextra(TSextra, samplingmethod, rnull)
  CustomTS=ifelse(missing(TS), FALSE, TRUE)
  WithWeights=ifelse(all(c(wx,wy)==1), FALSE, TRUE)
  if(min(wx, wy)<0.01 | max(wx, wy)>5)
      if(!SuppressMessages) message("Some of the weights are either exceptionally small or large. In either case the tests can be unreliable!")
  if(missing(UseLargeSample)) {
      UseLargeSample=ifelse(max(length(x), length(y))<1e4, FALSE, TRUE)
      if(!Continuous) UseLargeSample=FALSE
  }  
  if(!Continuous) vals=dta$vals
  if(CustomTS) tsinfo <- maketypeTS(TS, Continuous, WithWeights)
  else tsinfo <- maketypeTS(Continuous=Continuous, WithWeights=WithWeights)
  TS <- tsinfo$TS
  typeTS <- tsinfo$typeTS
# add the HHG test to the built-in tests (continuous data without weights, discrete data).
# For discrete data its null distribution is found from the first data set.
  if(!CustomTS && !WithWeights) TS <- add_HHG(TS, x, y, Continuous)

  if(!CustomTS) {
    if(!Continuous) TSextra$adw=weights(dta)
    tmp=calcTS(dta, TS, typeTS, TSextra)
    if(Continuous && WithWeights) {
      doMethods=c("KS","Kuiper","CvM","AD")
      if(min(wx,wy)<0.01 | max(wx,wy)>5)
        if(!SuppressMessages) message("Some of the weights are either exceptionally small or large. In either case the tests can be unreliable!")
    }
  } else {
    tmp=calcTS(dta, TS, typeTS, TSextra)
    if(is.list(tmp)) {
      message("TS should return a named vector of either test statistic(s) or p value(s), but not both.")
      return(NULL)
    }
    if(is.null(names(tmp))) {
      message("output of TS routine has to be a named vector!")
      return(NULL)
    }
  }
  methodnames=names(tmp)

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

# do chi square tests if built-in tests are used
  pwrchi=NULL
  if( (typeTS %in% c(1,2,5)) & (!UseLargeSample) ) {
      pwrchi = chi_power(rxy, alpha, B[1], avals, bvals, 
                           nbins, minexpcount, typeTS)
  }   

# Run tests for power
  if(UseLargeSample) {
       if(!SuppressMessages) message("Using methods with large sample theories")
       if(maxProcessor==1) { # no parallel processing 
           pwr=power_cont_LS(rxy, alpha, B[1], avals, bvals)
       }
       else {
          cl <- parallel::makeCluster(maxProcessor)
          on.exit(parallel::stopCluster(cl), add=TRUE)
          if(!is.null(seed))
            parallel::clusterSetRNGStream(cl, iseed=seed)
          z=parallel::clusterCall(cl, power_cont_LS, rxy=rxy, 
              alpha=alpha, B=round(B/maxProcessor) , xparam=avals, yparam=bvals)
          
          # Average power over cores  
          pwr=z[[1]]
          for(i in 2:maxProcessor) pwr=pwr+z[[i]]
          pwr = pwr/maxProcessor
       }
  }    
  else {
    if(With.p.value) {
       pwr=power_newtest(TS, rxy, xparam=avals, yparam=bvals ,
                         TSextra, alpha, B)
    }
    else {
       pwr=powerR(rxy=rxy, xparam=avals, yparam=bvals, TS=TS, typeTS, 
               TSextra, alpha=alpha, B=B,
               SuppressMessages = SuppressMessages,
               maxProcessor=maxProcessor, seed=seed)
    }
  }  
  if(!UseLargeSample) colnames(pwr) = methodnames
  if(!CustomTS) pwr = cbind(pwr, pwrchi)
  if(length(list(...))==0) rownames(pwr)=NULL
  if(length(list(...))==1) rownames(pwr)=avals
  if(length(list(...))==2) rownames(pwr)=paste0(avals,"|",bvals)  
  round(pwr, 3)
}
