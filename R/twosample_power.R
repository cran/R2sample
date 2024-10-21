#' Find the power of various two sample tests using Rcpp and parallel computing.
#' @param  f  function to generate a list with data sets x, y and (optional) vals, weights
#' @param  ... additional arguments passed to f, up to 2
#' @param  TS routine to calculate test statistics for non-chi-square tests
#' @param  TSextra additional info passed to TS, if necessary
#' @param  alpha =0.05, the level of the hypothesis test 
#' @param  B =c(1000, 2000), number of simulation runs for power and permutation test.
#' @param  nbins =c(50,10), number of bins for chi large and chi small.
#' @param  minexpcount =5 minimum required count for chi square tests
#' @param  UseLargeSample should p values be found via large sample theory if n,m>10000?
#' @param  samplingmethod =independence or MCMC in discrete data case
#' @param  maxProcessor =10, maximum number of cores to use. If maxProcessor=1 no parallel computing is used.
#' @return A numeric vector of power values.
#' @export 
#' @examples
#'  f=function(mu) list(x=rnorm(25), y=rnorm(25, mu))
#'  twosample_power(f, mu=c(0,2), B=c(100, 100), maxProcessor = 1)
#'  f=function(n, p) list(x=table(sample(1:5, size=1000, replace=TRUE)), 
#'        y=table(sample(1:5, size=n, replace=TRUE, 
#'        prob=c(1, 1, 1, 1, p))), vals=1:5)
#'  twosample_power(f, n=c(1000, 2000), p=c(1, 1.5), B=c(100, 100), maxProcessor = 1)

twosample_power=function(f, ..., TS, TSextra, alpha=0.05, B=c(1000, 1000), 
            nbins=c(50,10), minexpcount=5, UseLargeSample, 
            samplingmethod="independence", maxProcessor=10) {
       
  samplingmethod=ifelse(samplingmethod=="independence", 1, 2)
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
  tmp = rxy(avals[1], bvals[1])
  x=tmp$x
  y=tmp$y
  Continuous = ifelse("vals" %in% names(tmp), FALSE, TRUE)
  Weights=FALSE
  if(Continuous & length(tmp)==4) {
    Weights=TRUE
    wx=tmp$wx
    wy=tmp$wy
  }  
  if(!Continuous & any( abs(c(x,y)-round(c(x,y)))>1e-10 )) {
    Weights=TRUE
    if(min(x,y)<0.01 | max(x,y)>5)
      message("Some of the weights are either exceptionally small or large. In either case the tests can be unreliable!")
  } 
  if(missing(UseLargeSample))
    UseLargeSample=ifelse(min(length(x), length(y))<1e3, FALSE, TRUE)
  if(!Continuous) vals=tmp$vals
  if(missing(TS)) { # do included methods
    CustomTS=FALSE
    TSextra = list(aaaa=0)
    if(Continuous) { # Continuous Data
      if(!Weights) {
        typeTS = 1
        TS=TS_cont
        tmp=TS(x, y)
      }
      else {
        typeTS = 2
        TS=TSw_cont
        tmp=TS(x, y, wx, wy)
      }
    }
    else { # Discrete Data
      if(all(abs(c(x,y)-round(c(x,y)))<1e-10)) { # no weights
        typeTS=5
        TS=TS_disc
        tmp=TS(x, y, vals, rep(1, length(x)))
      }
      else { # with weights
        typeTS=6
        TS=TSw_disc
        tmp=TS(x, y)
        doMethods=c("KS", "Kuiper", "CvM", "AD")
      }
    }  
#   already do chi square tests if built-in tests are used
    if(typeTS %in% c(1,2,5,6)) {
       z = rxy(avals[1], bvals[1])
       pwrchi=NULL
       if(typeTS==6 & !("TSextra"%in%names(z))) pwrchi=NULL
       else {
         if(!UseLargeSample)
             pwrchi = chi_power(rxy, alpha, B[1], avals, bvals, 
                               nbins, minexpcount, typeTS)
       }   
    }   
    else pwrchi = NULL
  }
  else { # do user-supplied tests
    UseLargeSample=FALSE
    CustomTS=TRUE
    outchi=NULL # no chi square tests
    if(substr(deparse(TS)[2], 1, 5)==".Call") {
      if(maxProcessor>1) {
        message("Parallel Programming is not possible if custom TS is written in C++. Switching to single processor")  
        maxProcessor=1
      }  
    }
    if(Continuous) {
      if(length(formals(TS))==2) {
        typeTS=3
        TSextra = list(aaaa=0) # nonsense TSextra, so item exists
        tmp=TS(x, y)  # TS routine without TSextra
      }
      else {
        typeTS=4
        tmp=TS(x, y, TSextra) # TS routine with TSextra
      }
    }
    else {
      if(length(formals(TS))==3) {
        typeTS=7
        TSextra = list(aaaa=0)
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
  methodnames=names(tmp)

# Run tests for power
    m=parallel::detectCores()
    if(m==1 | maxProcessor==1) { # no parallel processing 
       if(m==1) message("No multiple cores found, using one core\n") 
       if(Continuous) {
          if(UseLargeSample) 
            pwr=power_cont_LS(rxy, alpha, B[1], avals, bvals)
          else                  
            pwr=power_cont(rxy=rxy, TS=TS, typeTS, TSextra, alpha=alpha, B=B, 
                xparam=avals, yparam=bvals)
          
       }    
       else
          pwr=power_disc(rxy=rxy, TS=TS, typeTS, TSextra, alpha=alpha, 
                samplingmethod=samplingmethod, B=B, xparam=avals, yparam=bvals)
    }  
    else { # Use one less processor than is present, or at most maxProcessor
      m=min(m, maxProcessor+1)-1
      cl <- parallel::makeCluster(m)
      if(Continuous) {
        if(UseLargeSample) 
          z=parallel::clusterCall(cl, power_cont_LS, rxy=rxy, 
                    alpha=alpha, B=round(B[1]/m), xparam=avals, yparam=bvals)
        else  
           z=parallel::clusterCall(cl, power_cont, rxy=rxy, TS=TS, typeTS, TSextra, alpha=alpha, 
                                   B=c(round(B[1]/m), B[2]), xparam=avals, yparam=bvals)
      }  
      else
        z=parallel::clusterCall(cl, power_disc, rxy=rxy, TS=TS, typeTS, TSextra, alpha=alpha, 
              samplingmethod=samplingmethod, B=c(round(B[1]/m), B[2]), xparam=avals, yparam=bvals)
      parallel::stopCluster(cl)  
      # Average power of cores  
      pwr=z[[1]]
      for(i in 2:m) pwr=pwr+z[[i]]
      pwr = pwr/m
    }
    if(UseLargeSample) 
      colnames(pwr) = c("KS", "Kuiper", "CvM", "AD", 
                        "ES large", "ES small" , "EP large", "EP small")
    else colnames(pwr) = methodnames
    rownames(pwr) = rnames
    if(!CustomTS) pwr = cbind(pwr, pwrchi)
    pwr
}
