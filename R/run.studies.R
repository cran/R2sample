#' Power Comparisons
#' 
#' This function runs the case studies included in the package and compares the
#' power of a new test to those included.
#' 
#' For details consult vignette("R2sample","R2sample")
#' 
#' @param TS routine to calculate test statistics.
#' @param study either the name of the study, or its number. If missing all the studies are run.
#' @param TSextra list passed to TS.
#' @param With.p.value =FALSE does user supplied routine return p values?
#' @param BasicComparison =TRUE if true compares tests on one default value of parameter of the alternative distribution.
#' @param nsample = 500, desired sample size.
#' @param alpha =0.05  type I error
#' @param param_alt (list of) values of parameter under the alternative hypothesis. If missing included values are used.
#' @param maxProcessor number of cores to use for parallel programming
#' @param SuppressMessages = FALSE print informative messages?
#' @param B = 1000
#' @return A (list of ) matrices of power values.
#' @examples
#' #The new test is a simple chisquare test:
#' chitest = function(x, y, TSextra) {
#'    nbins=TSextra$nbins
#'    nx=length(x);ny=length(y);n=nx+ny
#'    xy=c(x,y)
#'    bins=quantile(xy, (0:nbins)/nbins)
#'    Ox=hist(x, bins, plot=FALSE)$counts
#'    Oy=hist(y, bins, plot=FALSE)$counts
#'    tmp=sqrt(sum(Ox)/sum(Oy))
#'    chi = sum((Ox/tmp-Oy*tmp)^2/(Ox+Oy))
#'    pval=1-pchisq(chi, nbins-1)
#'    out=ifelse(TSextra$statistic,chi,pval)
#'    names(out)="ChiSquare"
#'    out
#' }
#' TSextra=list(nbins=5,statistic=FALSE) # Use 5 bins and calculate p values
#' run.studies(chitest,TSextra=TSextra, With.p.value=TRUE, B=100)
#' @export

run.studies <- function(TS, study, TSextra, With.p.value=FALSE, BasicComparison=TRUE, 
          nsample=500, alpha=0.05, param_alt, maxProcessor,
          SuppressMessages = FALSE, B=1000) {

  if(!is.function(TS)) {
      if(missing(TS) || !is.logical(TS)) {
         message("TS should either be a function, or TRUE/FALSE (for continuous or discrete) to run the included tests") 
         return(NULL)
      }   
      Continuous=TS
  }    
  else Continuous=ifelse("vals" %in% names(formals(TS)), FALSE, TRUE)
  I80cont=c(16, 17, 11, 14, 10, 10, 13, 10, 20, 3, 5, 11, 9, 11, 14, 16, 13, 17, 11, 16)
  I80disc=c(16, 12, 12, 14, 10, 10, 12, 11, 21, 5, 6, 12, 10, 11, 17, 16, 13, 18, 12, 13)
  if(Continuous) I80=I80cont
  else I80=I80disc
  if(missing(maxProcessor)) 
    maxProcessor=parallel::detectCores(logical = FALSE)-1
  NewParams=ifelse(alpha==0.05, FALSE, TRUE)
  if(!missing(param_alt)) {
      NewParams=TRUE
      if(!SuppressMessages) message("For new parameter values under the alternative power  values will also be calculated for included tests")
      if(!missing(study) && length(study)>1) {
        if(!is.list(param_alt)) {
           message("param_alt has to be a list with the same length as study")
           return(NULL)
        }
      }
      else if(!is.list(param_alt)) param_alt=list(param_alt)
  }    
  if(NewParams) BasicComparison=FALSE
  list.of.studies=c(
  "uniform.linear",
  "uniform.quadratic",
  "uniform.bump",
  "uniform.sine",
  "beta22.betaaa",
  "beta22.beta2a",
  "normal.shift",
  "normal.stretch",
  "normal.t",
  "normal.outlier1",
  "normal.outlier2",
  "exponential.gamma",
  "exponential.weibull",
  "exponential.bump",
  "gamma.normal",
  "normal.normalmixture",
  "uniform.uniformmixture",
  "uniform.betamixture",
  "chisquare.noncentral",
  "uniform.triangular"
  )
  if(Continuous) list.of.studies=paste0(list.of.studies,".cont")
  else list.of.studies=paste0(list.of.studies,".disc")
  if(missing(TSextra)) TSextra=list(aaa=1)
  if(missing(study)) study=1:20
  else BasicComparison=FALSE
  if(is.numeric(study)) study=list.of.studies[study]
  else {
    if(Continuous) study=paste0(study,".cont")
    else study=paste0(study,".disc")
  }
  out=as.list(seq_along(study))
  names(out)=study
  for(i in seq_along(study)) {
    if(!SuppressMessages) message(paste("Running case study", study[i],"..."))
    pwrold=R2sample::power_studies_results[[list.of.studies[i]]]
    tmp=case.studies(study[i], nsample)
    if(BasicComparison) tmp$param_alt=tmp$param_alt[I80[i]]
    if(NewParams || !is.function(TS)) {
       if(!missing(param_alt)) tmp$param_alt=param_alt[[i]]
       pwrold=R2sample::twosample_power(tmp$f, tmp$param_alt, alpha=alpha, 
                                        maxProcessor=maxProcessor, B=B) 
       if(!is.matrix(pwrold)) pwrold=rbind(pwrold)           
    }  
    if(!is.function(TS)) {out[[i]]=pwrold;next}  
    if(With.p.value) pwr=power_newtest(TS, tmp$f, tmp$param_alt, TSextra, alpha, B[1])
    else  pwr=R2sample::twosample_power(tmp$f, tmp$param_alt, 
                                   TS=TS, TSextra=TSextra, alpha=alpha, 
                                   maxProcessor=maxProcessor, B=B)                                                                    
    if(BasicComparison) out[[i]]=cbind(pwr, pwrold[I80[i], , drop=FALSE])
    else out[[i]]=cbind(pwr, pwrold[, , drop=FALSE])
  }
  if(BasicComparison) {
     A=matrix(0, 20, ncol(out[[1]]))
     for(i in 1:20) A[i, ]= out[[i]][1, ]
     rownames(A) = list.of.studies
     colnames(A) = colnames(out[[1]])
     a1=apply(A, 1, rank)
     message("Average number of studies a method is close to the best:")
     print(sort(apply(a1,1,mean)))
     return(A)
  }
  if(length(study)==1) return(out[[1]])
  out
}
