#' Find the power of various two sample tests using Rcpp and parallel computing.
#' @param  f  function to generate a list with data sets x, y and (optional) vals, weights
#' @param  ... additional arguments passed to f, up to 2
#' @param  TS routine to calculate test statistics for non-chi-square tests
#' @param  TSextra additional info passed to TS, if necessary
#' @param  alpha =0.05, the level of the hypothesis test 
#' @param  B =1000, number of simulation runs.
#' @param  nbins =c(50,10), number of bins for chi large and chi small.
#' @param  minexpcount =5 minimum required count for chi square tests
#' @param  UseLargeSample should p values be found via large sample theory if n,m>10000?
#' @param  samplingmethod =independence or MCMC in discrete data case
#' @param  rnull a function that generates data from a model, possibly with parameter estimation.
#' @param  SuppressMessages = FALSE print informative messages?
#' @param  maxProcessor maximum number of cores to use. If maxProcessor=1 no parallel computing is used.
#' @return A numeric vector of power values.
#' @export 
#' @examples
#'  f=function(mu) list(x=rnorm(25), y=rnorm(25, mu))
#'  twosample_power(f, mu=c(0,2), B=100, maxProcessor = 1)
#'  f=function(n, p) list(x=table(sample(1:5, size=1000, replace=TRUE)), 
#'        y=table(sample(1:5, size=n, replace=TRUE, 
#'        prob=c(1, 1, 1, 1, p))), vals=1:5)
#'  twosample_power(f, n=c(1000, 2000), p=c(1, 1.5), B=100, maxProcessor = 1)

twosample_power=function(f, ..., TS, TSextra, alpha=0.05, B=1000, 
            nbins=c(50,10), minexpcount=5, UseLargeSample, 
            samplingmethod="independence", rnull, 
            SuppressMessages=FALSE, maxProcessor) {

  if(!missing(UseLargeSample) & !missing(rnull)) {
      message("Large sample methods are not implemented for GOF-Twosample hybrid problem")
      return(NULL)
  }  

# create function rxy which generates data, with two arguments                       
    if(length(list(...))==0) { # f has 0 arguments
       rxy=function(a=0, b=0) f()
       avals=0
       bvals=0
    }
    if(length(list(...))==1) { # f has 1 argument
       rxy=function(a=0, b=0) f(a)
       avals=list(...)[[1]]
       bvals=0
    }
    if(length(list(...))==2) { # f has 2 arguments
       rxy=function(a=0, b=0) f(a,b)
       avals=list(...)[[1]]
       bvals=list(...)[[2]]
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
      rnames=bvals   
    }    
    else {
      bvals=rep(bvals, length(avals))
      rnames=avals
    }    
  }
  else rnames=1:length(avals)
# generate one data set as an example, do some setup 
  dta = rxy(avals[1], bvals[1])
  if(length(dta)<2 | length(dta)>4) {
    message("routine f has to create a list with 2 (continous data without weights), 3 (discrete data) or 4 (continous data with weights)")
    return(NULL)
  }
  if( (!("x"%in%names(dta))) | (!("x"%in%names(dta))) ) {
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

  samplingmethod=ifelse(samplingmethod=="independence", 1, 2)
  Continuous = ifelse("vals" %in% names(dta), FALSE, TRUE)
  if(missing(TSextra)) TSextra = list(samplingmethod=samplingmethod)
  else TSextra=c(TSextra, samplingmethod=samplingmethod)
  if(!missing(rnull)) TSextra=c(TSextra, rnull=rnull)
  CustomTS=ifelse(missing(TS), FALSE, TRUE)
  WithWeights=ifelse(all(c(wx,wy)==1), FALSE, TRUE)
  if(!Continuous & missing(maxProcessor)) {
     maxProcessor=1
     if(!SuppressMessages) message("For discrete data only a single processor is used if maxProcessor is not specified")
  }
  if(min(wx, wy)<0.01 | max(wx, wy)>5)
      if(!SuppressMessages) message("Some of the weights are either exceptionally small or large. In either case the tests can be unreliable!")
  if(missing(UseLargeSample)) {
      UseLargeSample=ifelse(max(length(x), length(y))<1e4, FALSE, TRUE)
      if(!Continuous) UseLargeSample=FALSE
  }  
  if(!Continuous) vals=dta$vals
  if(!CustomTS) { # do included methods
    if(Continuous) {
      if(!WithWeights) { # No weights
        typeTS = 1
        TS=TS_cont
        tmp=TS(x,y)
      }
      else { # weighted data
        typeTS = 2
        TS=TSw_cont
        tmp=TS(x, y, wx, wy)
        doMethods=c("KS", "Kuiper", "CvM", "AD")
        if(min(wx,wy)<0.01 | max(wx,wy)>5)
          if(!SuppressMessages) message("Some of the weights are either exceptionally small or large. In either case the tests can be unreliable!")
        }
      }
      else {
        typeTS=5
        TS=TS_disc
        TSextra=c(TSextra, list(adw=weights(dta)))
        tmp=TS(x, y, vals, TSextra$adw)
      }  
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
        tmp=TS(x, y)
      }
      else {
        typeTS=4
        tmp=TS(x, y, TSextra)
      }
    }
    else {
      if(length(formals(TS))==3) {
        typeTS=6
        tmp=TS(x, y, vals)
      }  
      else {
        typeTS=7
        tmp=TS(x, y, vals, TSextra)
      }
    }
    if(is.null(names(tmp))) {
      message("output of TS routine has to be a named vector!")
      return(NULL)
    }  
  }
  methodnames=names(tmp)
  
# do chi square tests if built-in tests are used
  pwrchi=NULL
  if( (typeTS %in% c(1,2,5)) & (!UseLargeSample) ) {
      pwrchi = chi_power(rxy, alpha, B[1], avals, bvals, 
                           nbins, minexpcount, typeTS)
  }   

# Run tests for power
  if(missing(maxProcessor) & Continuous) {
       maxProcessor=max(1, parallel::detectCores(logical = FALSE)-1)
  }   
  if(UseLargeSample) {
       if(!SuppressMessages) message("Using methods with large sample theories")
       if(maxProcessor==1) { # no parallel processing 
           pwr=power_cont_LS(rxy, alpha, B[1], avals, bvals)
       }
       else {
          cl <- parallel::makeCluster(maxProcessor)
          z=parallel::clusterCall(cl, power_cont_LS, rxy=rxy, 
              alpha=alpha, B=round(B/maxProcessor) , xparam=avals, yparam=bvals)
          parallel::stopCluster(cl)  
          # Average power over cores  
          pwr=z[[1]]
          for(i in 2:maxProcessor) pwr=pwr+z[[i]]
          pwr = pwr/maxProcessor
       }
  }    
  else {
    pwr=powerR(rxy=rxy, xparam=avals, yparam=bvals, TS=TS, typeTS, 
               TSextra, alpha=alpha, B=B,
               SuppressMessages = SuppressMessages, maxProcessor=maxProcessor)
  }  
  if(UseLargeSample) {
    colnames(pwr) = c("KS", "Kuiper", "CvM", "AD", 
                     "ES large", "ES small" , "EP large", "EP small")
  }  
  else colnames(pwr) = methodnames
  rownames(pwr) = rnames
  if(!CustomTS) pwr = cbind(pwr, pwrchi)
  if(nrow(pwr)==1) {
    nm=colnames(pwr) 
    pwr=c(pwr)
    names(pwr)=nm
  }  
  pwr
}
