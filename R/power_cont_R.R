#' Find the power of built-in continuous two sample tests using Rcpp and parallel computing.
#' @param  rxy  function to generate a list with data sets x, y and (optional) vals, weights
#' @param  xparam first argument passed to rxy
#' @param  yparam second argument passed to rxy
#' @param  TS test statistic
#' @param  typeTS which format has TS?
#' @param  TSextra list of items passed TS
#' @param  alpha =0.05, the level of the hypothesis test 
#' @param  B = 1000 number of simulation runs
#' @param  maxProcessor maximum number of cores to use. If maxProcessor=1 no parallel computing is used.
#' @return A numeric vector of power values.
#' @export 
power_cont_R=function(rxy, xparam, yparam, 
                    TS, typeTS, TSextra, alpha=0.05, 
                    B=1000, maxProcessor) {
   dta=rxy(xparam[1], yparam[1])
   if(typeTS==1) TS_data=TS(dta$x, dta$y)
   if(typeTS==2) TS_data=TS(dta$x, dta$y, dta$wx, dta$wy)
   if(typeTS==3) TS_data=TS(dta$x, dta$y)
   if(typeTS==4) TS_data=TS(dta$x, dta$y, TSextra)
   pwr=matrix(0, length(xparam), length(TS_data))
   colnames(pwr)=names(TS_data)
   rownames(pwr)=xparam
   if(maxProcessor>1) {
      tm=timecheck(dta, TS, typeTS, TSextra)
      if(tm*length(xparam)*2*B<20) maxProcessor=1
      message("maxProcessor set to 1 for faster computation")
   }
   if(maxProcessor>1) cl <- parallel::makeCluster(maxProcessor)
   if(maxProcessor==1) { # no parallel processing 
         tmp=power_cont(rxy=rxy, xparam, yparam, TS, 
                             typeTS, TSextra, B=B)
         Data=tmp$Data
         Permuted=tmp$Permuted
    }  
    else { 
        z=parallel::clusterCall(cl, power_cont, 
                              rxy=rxy,  xparam=xparam, yparam=yparam,
                              TS=TS, typeTS=typeTS, TSextra=TSextra, 
                              B=round(B/maxProcessor))
        Permuted=z[[1]][["Permuted"]]
        Data=z[[1]][["Data"]]
        for(i in 2:maxProcessor) {
          Permuted=rbind(Permuted,z[[i]][["Permuted"]])
          Data=rbind(Data,z[[i]][["Data"]])
        }  
    }
    for(i in seq_along(xparam)) {
       tmpD=Data[Data[,1]==xparam[i], -1, drop=FALSE]
       tmpP=Permuted[Permuted[,1]==xparam[i], -1, drop=FALSE]
       crtval=apply(tmpP, 2, quantile, 
                prob=1-alpha, na.rm=TRUE)
       for(j in seq_along(crtval)) 
         pwr[i, j]=sum(tmpD[ ,j]>crtval[j])/nrow(tmpD)
    }   
    round(pwr, 3)
}
