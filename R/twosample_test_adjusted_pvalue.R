#' Adjusted p values for simultaneous testing in the two-sample problem.
#'
#' This function runs a number of two sample tests using Rcpp and parallel computing and then finds the correct p value for the combined tests.
#'
#' For details consult vignette("R2sample","R2sample")
#' 
#' @param  x  a vector of numbers if data is continuous or of counts  if data is discrete, or a list with the data.
#' @param  y a vector of numbers if data is continuous or of counts  if data is discrete.
#' @param  vals =NA, a vector of numbers, the values of a discrete random variable. NA if data is continuous data.
#' @param  TS routine to calculate test statistics for non-chi-square tests
#' @param  TSextra additional info passed to TS, if necessary
#' @param  wx A numeric vector of weights of x.
#' @param  wy A numeric vector of weights of y.
#' @param  B =c(5000, 1000), number of simulation runs for permutation test
#' @param  seed =NULL, optional seed for reproducibility
#' @param  nbins =c(50,10), number of bins for chi square tests.
#' @param  minexpcount = 5, minimum required expected counts for chi-square tests
#' @param  samplingmethod ="independence" or "Binomial" for discrete data
#' @param  rnull routine for parametric bootstrap
#' @param  SuppressMessages = FALSE print informative messages?
#' @param  doMethods ="all" a vector of codes for the methods to include. If "all", all methods are used.
#' @return (invisibly) a list with the p values of the individual tests and the adjusted p value of the combined test. These are also printed as messages.
#' @examples
#'  x=rnorm(100)
#'  y=rt(200, 4)
#'  R2sample::twosample_test_adjusted_pvalue(x, y, B=c(500, 500))
#'  vals=1:5
#'  x=table(c(1:5, sample(1:5, size=100, replace=TRUE)))-1
#'  y=table(c(1:5, sample(1:5, size=100, replace=TRUE, prob=c(1,1,3,1,1))))-1
#'  R2sample::twosample_test_adjusted_pvalue(x, y, vals, B=c(500, 500))
#' @export
twosample_test_adjusted_pvalue=function(x, y, vals=NA, TS, TSextra, wx=rep(1, length(x)),
                        wy=rep(1, length(y)), B=c(5000, 1000), 
                        seed=NULL, nbins=c(50,10),
                        minexpcount=5, samplingmethod="independence", 
                        rnull, SuppressMessages=FALSE, doMethods) {
    default.methods = list(cont=c("ES small", "ZA", "ZK", "Wassp1","Kuiper"), 
                           disc=c("small", "ZA", "Kuiper", "Wassp1"))
    all.methods = list(cont=c("KS","Kuiper","CvM","AD","ZA","ZK","ZC","Wassp1","EppsS","HHG",
                             "ES large","ES small","EP large","EP small"),
                       disc=c("KS","Kuiper","CvM","AD","ZA","Wassp1","EppsS","HHG","large","small"))                                          
    if(is.list(x)) {
      y=x$y
      if("vals"%in%names(x)) vals=x$vals
      if("wx"%in%names(x)) wx=x$wx
      if("wy"%in%names(x)) wy=x$wy
      x=x$x
    }
    if(!is.numeric(samplingmethod))
      samplingmethod=ifelse(samplingmethod=="independence", 1, 2)
    if(!any(is.na(vals)) && is.unsorted(vals)) {
      I=order(vals)
      x=x[I]
      y=y[I]
      vals=vals[I]
    }
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
# add the HHG test to the built-in tests (continuous data without weights, discrete data)
    if(!CustomTS && !WithWeights) TS <- add_HHG(TS, x, y, Continuous)
    if(length(B)==1) B=c(B, B)
    if(missing(doMethods)) {
      if(Continuous) doMethods=default.methods$cont
      else doMethods=default.methods$disc
    }
    if(doMethods[1]=="all") {
      if(Continuous) doMethods=all.methods$cont
      else doMethods=all.methods$disc
    }
    if(test_methods(doMethods,Continuous, FALSE, WithWeights))
      return(NULL)
# what methods are to be run?  
    if(!CustomTS) {
      if(Continuous) {
        if(WithWeights) dta=list(x=x,y=y,wx=wx,wy=wy)
        else dta=list(x=x,y=y)
        tmp=calcTS(dta, TS, typeTS, TSextra)
        if(WithWeights && (min(wx,wy)<0.01 | max(wx,wy)>5))
          if(!SuppressMessages) message("Some of the weights are either exceptionally small or large. In either case the tests can be unreliable!")
      } else {
        dta=list(x=x,y=y,vals=vals)
        TSextra$adw=weights(dta)
        tmp=calcTS(dta, TS, typeTS, TSextra)
      }
    } else {
      if(Continuous) dta=list(x=x,y=y)
      else dta=list(x=x,y=y,vals=vals)
      tmp=calcTS(dta, TS, typeTS, TSextra)
      if(is.null(names(tmp))) {
        message("output of TS routine has to be a named vector!")
        return(NULL)
      }
    }

#   Start the first reproducible Monte Carlo stream only after all setup and
#   calculation of the observed statistic have been completed.
    set_seed(seed)

#   Find matrix of values of test statistics for simulated data
    num_tests=length(tmp)
    A=matrix(0, B[1], num_tests)
    for(i in 1:B[1]) {
       simdta=gen_sim_data(dta, TSextra)
       A[i, ]=calcTS(simdta, TS, typeTS, TSextra)
    }  
    pvalsTS=matrix(0, B[2]+1, num_tests)
    colnames(pvalsTS)=names(tmp)
    chi.methods <- if(Continuous)
      c("ES large","ES small","EP large","EP small")
    else
      c("large","small")
    needCHI <- !CustomTS && any(doMethods %in% chi.methods)

    if(needCHI) {
      pvalsCHI=matrix(0, B[2]+1, length(chi.methods))
      colnames(pvalsCHI)=chi.methods
    } else {
      pvalsCHI=NULL
    }

    for(i in 1:(B[2]+1)) {
      if(i==1) {
        tmp=calcTS(dta, TS, typeTS, TSextra)
        if(needCHI)
          pvalsCHI[i, ]=chi_test(dta, nbins=nbins, minexpcount=minexpcount, typeTS=typeTS)$p.values
      }  
      else {
        simdta=gen_sim_data(dta, TSextra)
        tmp=calcTS(simdta, TS, typeTS, TSextra)
        if(needCHI)
          pvalsCHI[i, ]=chi_test(simdta, nbins=nbins, minexpcount=minexpcount, typeTS=typeTS)$p.values
      }  
# p value = proportion of simulated statistics at least as large as the observed one
      for(j in 1:num_tests) 
        pvalsTS[i,j]=sum(A[,j]>=tmp[j]-1e-9*abs(tmp[j]))/B[1]    
    }
    pvals=cbind(pvalsTS, pvalsCHI) 
    pvals=pvals[,doMethods,drop=FALSE]
    minp_x=min(pvals[1, ])
    minp_sim=apply(pvals[-1, ,drop=FALSE], 1, min)
# adjusted p value = P(smallest p value <= observed smallest p value | H0)
    minp_adj=round(mean(minp_sim<=minp_x), 4)
    message("p values of individual tests:")
    for(i in 1:ncol(pvals)) message(paste(colnames(pvals)[i],": ", pvals[1,i]))
    message(paste0("adjusted p value of combined tests: ", minp_adj))
    invisible(list(p.values=pvals[1, ], adjusted.p.value=minp_adj))
    
}
