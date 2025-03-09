#' This function runs a number of two sample tests using Rcpp and parallel computing and then finds the correct p value for the combined tests.
#' @param  x  a vector of numbers if data is continuous or of counts  if data is discrete.
#' @param  y a vector of numbers if data is continuous or of counts  if data is discrete.
#' @param  vals =NA, a vector of numbers, the values of a discrete random variable. NA if data is continuous data.
#' @param  TS routine to calculate test statistics for non-chi-square tests
#' @param  TSextra additional info passed to TS, if necessary
#' @param  wx A numeric vector of weights of x.
#' @param  wy A numeric vector of weights of y.
#' @param  B =c(5000, 1000), number of simulation runs for permutation test
#' @param  nbins =c(50,10), number of bins for chi square tests.
#' @param  minexpcount = 5, minimum required expected counts for chi-square tests
#' @param  samplingmethod ="independence" or "MCMC" for discrete data
#' @param  doMethods  Which methods should be included? 
#' @return A list of two numeric vectors, the test statistics and the p values. 
#' @export 
#' @examples
#'  x=rnorm(100)
#'  y=rt(200, 4)
#'  R2sample::twosample_test_adjusted_pvalue(x, y, B=c(500, 500))
#'  vals=1:5
#'  x=table(c(1:5, sample(1:5, size=100, replace=TRUE)))-1
#'  y=table(c(1:5, sample(1:5, size=100, replace=TRUE, prob=c(1,1,3,1,1))))-1
#'  R2sample::twosample_test_adjusted_pvalue(x, y, vals, B=c(500, 500))

twosample_test_adjusted_pvalue=function(x, y, vals=NA, TS, TSextra, wx=rep(1, length(x)),
                        wy=rep(1, length(y)), B=c(5000, 1000), nbins=c(50,10),
                        minexpcount=5, samplingmethod="independence", doMethods) {
    default.methods = list(cont=c("ES small", "ZA", "ZK", "Wassp1","Kuiper"), 
                           disc=c("small", "ZA", "Kuiper", "Wassp1"))
    all.methods = list(cont=c("KS","Kuiper","CvM","AD","LR","ZA","ZK","ZC","Wassp1",
                             "ES large","ES small","EP large","EP small"),
                       disc=c("KS","Kuiper","CvM","AD","LR","ZA","Wassp1","large","small"))                                          
    samplingmethod=ifelse(samplingmethod=="independence", 1, 2)
    if(length(B)==1) B=c(B, B)
    if(missing(TSextra)) TSextra = list(aaaa=0)
    if(any(is.na(vals))) Continuous=TRUE
    else Continuous=FALSE
    if(missing(doMethods)) {
      if(Continuous) doMethods=default.methods$cont
      else doMethods=default.methods$disc
    }
    if(doMethods[1]=="all") {
      if(Continuous) doMethods=all.methods$cont
      else doMethods=all.methods$disc
    }
# what methods are to be run?  
    CustomTS=TRUE
    if(missing(TS)) { # do included methods
      CustomTS=FALSE
      if(Continuous) {
          if(all( abs(c(wx,wy)-round(c(wx,wy)))<1e-10 )) { # No weights
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
               message("Some of the weights are either exceptionally small or large. In either case the tests can be unreliable!")
          }
          
      }
      else {
          dta=list(x=x, y=y, vals=vals)
          if(all( abs(c(x,y)-round(c(x,y)))<1e-10 )) { # no weights
            typeTS=5
            TS=TS_disc
            adw=weights(dta)
            tmp=TS(x, y, vals, adw)
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
            }  
            doMethods=c("KS", "Kuiper", "CvM", "AD")
          }
      }  
    }  
    else { # do user-supplied tests
  
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
        dta=list(x=x, y=y, vals=vals)
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
      doMethods=names(tmp)
    }
#   Find matrix of values of test statistics after permutation
    if(Continuous) {n=length(x);m=length(y);xy=c(x,y);Ixy=1:length(xy);wxy=c(wx,wy)}
    num_tests=length(tmp)
    A=matrix(0, B[1], num_tests)
    for(i in 1:B[1]) {
         if(typeTS<=4) {
              I=sample(Ixy)
              if(typeTS==1) A[i, ]=TS(xy[I[1:n]], xy[I[(n+1):(n+m)]]) 
              if(typeTS==2) A[i, ]=TS(xy[I[1:n]], xy[I[(n+1):(n+m)]], wxy[I[1:n]], wxy[I[(n+1):(n+m)]])
              if(typeTS==3) A[i, ]=TS(xy[I[1:n]], xy[I[(n+1):(n+m)]])
              if(typeTS==4) A[i, ]=TS(xy[I[1:n]], xy[I[(n+1):(n+m)]], TSextra)
         
         }
         else {
              dta=permute_disc(dta, samplingmethod)
              if(typeTS==5) A[i, ]=TS(dta$x, dta$y, vals, adw)
              if(typeTS==6) A[i, ]=TS(dta$x, dta$y)    
              if(typeTS==7) A[i, ]=TS(dta$x, dta$y, vals)
              if(typeTS==8) A[i, ]=TS(dta$x, dta$y, vals, TSextra)
         
         }
    } 
    pvalsTS=matrix(0, B[2]+1, num_tests)
    colnames(pvalsTS)=names(tmp)
    pvalsCHI=matrix(0, B[2]+1, ifelse(Continuous, 4, 2))
    if(Continuous) colnames(pvalsCHI)=c("ES large","ES small","EP large","EP small")
    else colnames(pvalsCHI)=c("large","small")
    for(i in 1:(B[2]+1)) {
          if(Continuous) {
              if(i==1) I=Ixy #For data
              else I=sample(Ixy)
              if(typeTS==1) tmp=TS(xy[I[1:n]], xy[I[(n+1):(n+m)]]) 
              if(typeTS==2) tmp=TS(xy[I[1:n]], xy[I[(n+1):(n+m)]], wxy[I[1:n]], wxy[I[(n+1):(n+m)]])
              if(typeTS==3) tmp=TS(xy[I[1:n]], xy[I[(n+1):(n+m)]])
              if(typeTS==4) tmp=TS(xy[I[1:n]], xy[I[(n+1):(n+m)]], TSextra)
         }  
         else {
              if(i==1) dta=list(x=x,y=y,vals=vals)
              else dta=permute_disc(dta, samplingmethod)
              if(typeTS==5) tmp=TS(dta$x, dta$y, vals, adw)
              if(typeTS==6) tmp=TS(dta$x, dta$y)    
              if(typeTS==7) tmp=TS(dta$x, dta$y, vals)
              if(typeTS==8) tmp=TS(dta$x, dta$y, vals, TSextra)
         
         }         
         for(j in 1:num_tests) pvalsTS[i,j]=pvalsTS[i,j]+sum(tmp[j]>A[,j])/B[2]    
         if(!CustomTS) {
            if(Continuous) {
               if(i==1) dta=list(x=x,y=y)
               else dta=list(x=xy[I[1:n]], y=xy[I[(n+1):(n+m)]])
            }   
            pvalsCHI[i, ]=chi_test(dta, nbins=nbins, minexpcount=minexpcount, typeTS=typeTS)$p.values
         }
         else pvalsCHI=NULL
    }
    pvals=cbind(pvalsTS, pvalsCHI) 
    pvals=pvals[,doMethods,drop=FALSE]
    minp_x=min(pvals[1, ])
    minp_sim=apply(pvals[-1, ,drop=FALSE], 1, min)
    z=seq(0, 1, length=250)
    y=z
    for(i in 1:250) y[i]=sum(minp_sim<=z[i])/B[2]
    I=c(1:250)[z>minp_x][1]-1
    slope=(y[I+1]-y[I])/(z[I+1]-z[I])
    minp_adj=round(y[I]+slope*(minp_x-z[I]),4)
    message("p values of individual tests:")
    for(i in 1:ncol(pvals)) message(paste(colnames(pvals)[i],": ", pvals[1,i]))
    message(paste0("adjusted p value of combined tests: ", minp_adj))
    
}
