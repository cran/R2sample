#' Find the power of built-in continuous two sample tests using Rcpp and parallel computing.
#' @param  rxy  function to generate a list with data sets x, y and (optional) vals, weights
#' @param  xparam first argument passed to rxy
#' @param  yparam second argument passed to rxy
#' @param  TS test statistic
#' @param  typeTS which format has TS?
#' @param  TSextra list of items passed TS
#' @param  alpha =0.05, the level of the hypothesis test 
#' @param  samplingmethod =independence or MCMC in discrete data case
#' @param  B = 1000 number of simulation runs
#' @param  maxProcessor maximum number of cores to use. If maxProcessor=1 no parallel computing is used.
#' @return A numeric vector of power values.
#' @export 
power_disc_R=function(rxy, xparam, yparam, 
                    TS, typeTS, TSextra, alpha=0.05, 
                    samplingmethod=1,
                    B=1000, maxProcessor) {
   dta=rxy(xparam[1], yparam[1])
   if(typeTS==5) TS_data=TS(dta$x, dta$y, dta$vals, weights(dta))
   if(typeTS==6) TS_data=TS(dta$x, dta$y)
   if(typeTS==7) TS_data=TS(dta$x, dta$y, dta$vals)
   if(typeTS==8) TS_data=TS(dta$x, dta$y, dta$vals, TSextra)
   pwr=matrix(0, length(xparam), length(TS_data))
   colnames(pwr)=names(TS_data)
   rownames(pwr)=xparam
   
   if(maxProcessor>1) cl <- parallel::makeCluster(maxProcessor)
   if(maxProcessor==1) { # no parallel processing 
         tmp=power_disc(rxy=rxy, xparam, yparam, TS, 
                        typeTS, TSextra,samplingmethod, B)
         Data=tmp$Data
         Permuted=tmp$Permuted
    }  
    else { 
        z=parallel::clusterCall(cl, power_disc, 
                   rxy=rxy,  xparam=xparam, yparam=yparam,
                  TS=TS, typeTS=typeTS, TSextra=TSextra,
                  samplingmethod=samplingmethod,
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
