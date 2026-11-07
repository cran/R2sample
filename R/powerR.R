#' Find the power of two sample tests using Rcpp and parallel computing.
#' @param  rxy  function to generate a list with data sets x, y and (optional) vals, weights
#' @param  xparam first argument passed to rxy
#' @param  yparam second argument passed to rxy
#' @param  TS test statistic
#' @param  typeTS which format has TS?
#' @param  TSextra list of items passed TS
#' @param  alpha =0.05, the level of the hypothesis test 
#' @param  B = 1000 number of simulation runs
#' @param  SuppressMessages = FALSE print informative messages?
#' @param  maxProcessor maximum number of cores to use. If maxProcessor=1 no parallel computing is used.
#' @param seed =NULL, optional seed for reproducibility
#' @return A numeric vector of power values.
#' @keywords internal
#' @noRd
powerR=function(rxy, xparam, yparam,TS, typeTS, TSextra, alpha=0.05, 
                    B=1000, SuppressMessages, maxProcessor, seed=NULL) {
   set_seed(seed)
   dta=rxy(xparam[1], yparam[1])
   if(typeTS==5) TSextra$adw=weights(dta)
   TS_data = calcTS(dta, TS, typeTS, TSextra);
   pwr=matrix(0, length(xparam), length(TS_data))
   colnames(pwr)=names(TS_data)
   rownames(pwr)=xparam
   if(maxProcessor > 1) {
     tm = timecheck(dta, TS, typeTS, TSextra)
     maxProcessor = use_parallel(
       tm = tm,
       nparam = length(xparam),
       B = B,
       maxProcessor = maxProcessor
     )
     
     if(maxProcessor == 1) {
       if(!SuppressMessages)
         message("Using one core for faster computation")
     }
     else {
       if(!SuppressMessages)
         message(paste("Using", maxProcessor, "cores.."))
     }
   }
   if(maxProcessor>1) {
     cl <- parallel::makeCluster(maxProcessor)
     on.exit(parallel::stopCluster(cl), add=TRUE)
     if(!is.null(seed))
       parallel::clusterSetRNGStream(cl, iseed=seed)
   }   
   if(maxProcessor==1) { # no parallel processing 
         tmp=powerC(rxy=rxy, xparam, yparam, TS, typeTS, TSextra, B=B)
         Data=tmp$Data
         Simulated=tmp$Simulated
         paramalt=tmp$paramalt
    }  
    else { 
        z=parallel::clusterCall(cl, powerC, 
                              rxy=rxy,  xparam=xparam, yparam=yparam,
                              TS=TS, typeTS=typeTS, TSextra=TSextra, 
                              B=round(B/maxProcessor))
        Simulated=z[[1]][["Simulated"]]
        Data=z[[1]][["Data"]]
        paramalt=z[[1]][["paramalt"]]
        for(i in 2:maxProcessor) {
          Simulated=rbind(Simulated,z[[i]][["Simulated"]])
          Data=rbind(Data,z[[i]][["Data"]])
          paramalt=rbind(paramalt,z[[i]][["paramalt"]])
        }  
    }
    for(i in seq_along(xparam)) {
       Index=c(1:nrow(Data))[paramalt[,1]==xparam[i]&paramalt[,2]==yparam[i]]
       tmpD=Data[Index, , drop=FALSE]
       tmpS=Simulated[Index, , drop=FALSE]
# For each simulated data set find the p value against the simulated null
# distribution, counting ties as in twosample_test, and reject if p<alpha
       for(j in 1:ncol(tmpD)) {
         sims=sort(tmpS[, j])
         obs=tmpD[, j]
         nless=findInterval(obs-1e-9*abs(obs), sims, left.open=TRUE)
         pv=(length(sims)-nless)/length(sims)
         pwr[i, j]=mean(pv<alpha)
       }
    }   
    pwr
}
