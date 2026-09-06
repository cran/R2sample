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
#' @param SuppressMessages = TRUE
#' @return A matrix of power values.
#' @examples
#' TSextra=list(nbins=5, statistic=FALSE) # Use 5 bins and calculate p values
#' run.studies(1:2, 1:2, TS=myTS_cont, TSextra=TSextra, With.p.value=TRUE, B=100)
#' @export
run.studies=function(studies=1:40, types=1:5, data_type="cont",
                     TS, TSextra, With.p.value =FALSE, n=500, 
                     alpha=0.05, B=1000, 
                     SuppressMessages = TRUE) {
  
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
  if(!is.null(k)) studytype=studytype[-k,]
  studytype=unique(studytype)
  num_old_tests=ncol(R2sample::power_study[[data_type]][[1]])
  num_studies=nrow(studytype)
  if(missing(TSextra)) TSextra=list(aaa=NULL)
  NewTest=ifelse(missing(TS), FALSE, TRUE)
  num_new_test=0
  out_old_test=matrix(0, nrow(studytype), num_old_tests)
  if(data_type=="cont") {
    colnames(out_old_test)=colnames(R2sample::power_study[[1]][[1]])
  }  
  else colnames(out_old_test)=colnames(R2sample::power_study[[2]][[1]])
  rownames(out_old_test)=paste0(studytype[,1], "/", studytype[,2])
  out_new_test=NULL
  for(i in 1:nrow(studytype)) {
      tmp = R2sample::case_studies(
        which=studytype[i,1],
        which_diff = studytype[i, 2],
        n = n,
        data_type = data_type,
        funs_list = R2sample::funs_list[[data_type]][[studytype[i,1]]][[studytype[i,2]]]
      )
      if(data_type=="cont") tmp$vals=NA
      if(NewTest) {
        pwr=R2sample::twosample_power(tmp$f, TS=TS, TSextra= TSextra,
                    With.p.value=With.p.value, alpha=alpha, 
                    B=B, SuppressMessages = SuppressMessages)
        num_new_test=length(pwr)
        if(i==1) {
          out_new_test=matrix(0, nrow(studytype), num_new_test)
          colnames(out_new_test)=names(pwr)
          rownames(out_new_test)=names(out_old_test)
        }
        out_new_test[i, ]=pwr 
      }    
      if(n==500&&alpha==0.05) {
          out_old_test[i, ]=
            R2sample::power_study[[data_type]][[studytype[i,2]]][studytype[i, 1], ]
      }
      else {
        out_old_test[i, ]=R2sample::twosample_power(tmp$f, 
                  alpha=alpha, B=B, SuppressMessages = SuppressMessages)
      }
  }  
  out = cbind(out_new_test, out_old_test)
  a1=apply(out, 1, rank)
  message("Average rank of a method:")
  print(round(100*sort(apply(a1,1,mean), decreasing = TRUE)/nrow(a1),1))
  out
}
