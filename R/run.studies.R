#' Power Comparisons
#' 
#' This function runs the case studies included in the package and compares the
#' power of a new test to those included.
#' 
#' For details consult vignette("R2sample","R2sample")
#'
#' @param studies =1:40, either the names of the studies, or their numbers. If missing all the studies are run. 
#' @param types =1:5, either the types of studies, or their numbers. If missing all types are run. 
#' @param data_type ="cont", continuous or discrete data
#' @param TS routine to calculate test statistics.
#' @param TSextra list passed to TS.
#' @param With.p.value =FALSE does user supplied routine return p values?
#' @param n = 500, desired sample size.
#' @param alpha =0.05  type I error
#' @param B = 1000
#' @param seed =NULL, optional seed for reproducibility. Case study i in the list of
#'        studies/types that are run uses the seed seed+i, both for the new test and
#'        (if they are recomputed) for the included methods.
#' @param SuppressMessages = TRUE
#' @param RerunIncludedCases =FALSE, recompute the powers of the included methods even
#'        if n=500 and alpha=0.05, instead of using the stored results in power_study
#' @param maxProcessor =1, number of processors used in twosample_power. Results with a
#'        given seed are reproducible only with the same number of processors.
#' @param Null =FALSE, if TRUE both data sets are generated under the null hypothesis
#'        (by rnull), so the result is the type I error. The type "Null" does the same
#'        for case studies 1-20.
#' @param ShowProgress =FALSE, if TRUE a message with the elapsed time and an estimate
#'        of the remaining time is shown after each case study.
#' @return A matrix of power values.
#' @examples
#' TSextra=list(nbins=5, statistic=FALSE) # Use 5 bins and calculate p values
#' run.studies(1:2, 1:2, TS=myTS_cont, TSextra=TSextra, With.p.value=TRUE, B=100)
#' @export
run.studies=function(studies=1:40, types=1:5, data_type="cont",
                     TS, TSextra, With.p.value =FALSE, n=500, 
                     alpha=0.05, B=1000, seed=NULL,
                     SuppressMessages = TRUE, RerunIncludedCases=FALSE,
                     maxProcessor=1, Null=FALSE, ShowProgress=FALSE) {
  # data generator under the null hypothesis for case study tmp
  null_f=function(tmp) {
    if(data_type=="cont") {
      rn=tmp$rnull
      return(function(a) list(x=rn(n), y=rn(n)))
    }
    # discrete: vals are the midpoints of equal-width bins
    v=tmp$vals
    h=diff(v)[1]
    p0=diff(tmp$pnull(c(v[1]-h/2, v+h/2)))
    function(a) {
      x=c(stats::rmultinom(1, n, p0))
      y=c(stats::rmultinom(1, n, p0))
      I=(x+y)>0
      list(x=x[I], y=y[I], vals=v[I])
    }
  }
  
# check if old version is assumed
  if(is.function(studies)) {
    stop("the arguments list of the function run.studies has changed.\n 
          Please check the help file or the vignette",
         call.=FALSE)
  }
  list.of.studies=R2sample::case_studies(ReturnCaseNames = TRUE)
  if(is.numeric(studies[1])) studies=list.of.studies[studies]
  all_types=names(R2sample::power_study[[1]])
  if(is.numeric(types)) types=all_types[types]
  studytype=cbind(rep(studies, length(types)),
                      rep(types, each=length(studies)))
  k=NULL
  for(i in 1:nrow(studytype)) {
    if(studytype[i,1]%in%list.of.studies[1:20] && studytype[i,2]=="mixed")
       k=c(k, i)
    if(studytype[i,1]%in%list.of.studies[21:40]) {
       studytype[i,2]="mixed"
    }  
  }  
  if(!is.null(k)) studytype=studytype[-k, , drop=FALSE]
  studytype=unique(studytype)
# methods in the stored power studies, LR has been removed from the package
  old_names=setdiff(colnames(R2sample::power_study[[data_type]][[1]]), "LR")
  num_old_tests=length(old_names)
  num_studies=nrow(studytype)
  if(missing(TSextra)) TSextra=list(aaa=NULL)
  NewTest=ifelse(missing(TS), FALSE, TRUE)
  num_new_test=0
  out_old_test=matrix(0, nrow(studytype), num_old_tests)
  colnames(out_old_test)=old_names
  rownames(out_old_test)=paste0(studytype[,1], "/", 
                                if(Null) "Null" else studytype[,2])
  out_new_test=NULL
  starttime=Sys.time()
  for(i in 1:nrow(studytype)) {
      ty = studytype[i, 2]
      isNull = Null || ty=="Null"
      # the type Null uses the case study of type mean with null data
      wd = ifelse(ty=="Null", "mean", ty)
      tmp = R2sample::case_studies(
        which=studytype[i,1],
        which_diff = wd,
        n = n,
        data_type = data_type,
        funs_list = R2sample::funs_list[[data_type]][[studytype[i,1]]][[wd]]
      )
      if(isNull) tmp$f = null_f(tmp)
      if(data_type=="cont") tmp$vals=NA
      # stored results: under the null hypothesis only for case studies 1-20
      stored_type = ifelse(isNull, "Null", ty)
      stored = R2sample::power_study[[data_type]][[stored_type]]
      if(NewTest) {
        pwr=R2sample::twosample_power(tmp$f, tmp$apar,
                     TS=TS, TSextra= TSextra,
                    With.p.value=With.p.value, alpha=alpha, 
                    B=B, seed=if(is.null(seed)) NULL else seed+i, 
                    SuppressMessages = SuppressMessages,
                    maxProcessor=maxProcessor)
        num_new_test=ncol(pwr)
        if(i==1) {
          out_new_test=matrix(0, nrow(studytype), num_new_test)
          colnames(out_new_test)=colnames(pwr)
          rownames(out_new_test)=rownames(out_old_test)
        }
        out_new_test[i, ]=pwr 
      }    
      if(n==500&&alpha==0.05&&!RerunIncludedCases&&
         studytype[i, 1]%in%rownames(stored)) {
          out_old_test[i, ]=stored[studytype[i, 1], old_names]
      }
      else {
        pwr_old=R2sample::twosample_power(
                  tmp$f, tmp$apar,                                   
                  alpha=alpha, B=B, seed=if(is.null(seed)) NULL else seed+i,
                  SuppressMessages = SuppressMessages,
                  maxProcessor=maxProcessor)
        out_old_test[i, ]=pwr_old[1, old_names]
      }
      if(ShowProgress) {
        elapsed=as.numeric(difftime(Sys.time(), starttime, units="secs"))
        left=elapsed/i*(num_studies-i)
        fmt=function(s) if(s<120) paste(round(s), "sec") else paste(round(s/60, 1), "min")
        message(i, "/", num_studies, " ", rownames(out_old_test)[i], 
                ": elapsed ", fmt(elapsed), ", remaining about ", fmt(left))
      }
  }  
  out = cbind(out_new_test, out_old_test)
  a1=apply(out, 1, rank)
  message("Average rank of a method:")
  print(round(100*sort(apply(a1,1,mean), decreasing = TRUE)/nrow(a1),1))
  out
}
