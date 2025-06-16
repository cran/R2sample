#' This function checks whether the correct methods have been requested
#' @param  doMethods ="all" Which methods should be included?
#' @param  Continuous is data continuous
#' @param  UseLargeSample should p values be found via large sample theory?
#' @param  WithWeights with weights?
#' @keywords internal
#' @return TRUE or FALSE
#' @export 
test_methods=function(doMethods, Continuous, UseLargeSample, WithWeights) {
    if(doMethods[1]%in%c("all","default")) return(FALSE)
    if(Continuous & !UseLargeSample &!WithWeights) 
       methods=c("KS","Kuiper","CvM","AD","LR","ZA","ZK","ZC",
                 "Wassp1","ES large","ES small","EP large","EP small")
    if(Continuous & !UseLargeSample &WithWeights) 
       methods=c("KS","Kuiper","CvM","AD",
                "ES large","ES small","EP large","EP small")
    if(Continuous & UseLargeSample &!WithWeights) 
       methods=c("KS","Kuiper","CvM","AD")
    if(!Continuous) 
       methods=c("KS","Kuiper","CvM","AD","LR","ZA", "Wassp1","large","small")
    Good=TRUE
    for(i in seq_along(doMethods)) {
      if(!(doMethods[i]%in%methods)) {Good=FALSE;break}
    }
    if(Good) return(FALSE)
    message(paste0(doMethods[i]," is not an included Method!"))
    if(Continuous & !UseLargeSample &!WithWeights) {
         message("For continuous data without weights included methods are")
         message("Method               Code")
         message("Kolmogorov-Smirnov   KS")
         message("Kuiper               Kuiper")
         message("Cramer-vonMises      CvM")
         message("Anderson-Darling     AD")
         message("Lehmann-Rosenblatt   LR")
         message("Zhang's tests        ZA, ZK and ZC")
         message("Wasserstein          Wassp1")
         message("Chi square tests     ES large, ES small")
         message("                     EP large, EP small")
    }
    if(Continuous & !UseLargeSample &WithWeights) {
      message("For continuous data with weights included methods are")
      message("Method               Code")
      message("Kolmogorov-Smirnov   KS")
      message("Kuiper               Kuiper")
      message("Cramer-vonMises      CvM")
      message("Anderson-Darling     AD")
      message("Chi square tests     ES large, ES small")
      message("                     EP large, EP small")
    }
    if(Continuous & UseLargeSample &!WithWeights) {
      message("For continuous data using large sample theory included methods are")
      message("Method               Code")
      message("Kolmogorov-Smirnov   KS")
      message("Kuiper               Kuiper")
      message("Cramer-vonMises      CvM")
      message("Anderson-Darling     AD")
    }
    if(!Continuous) {
      message("For discrete data included methods are")
      message("Method               Code")
      message("Kolmogorov-Smirnov   KS")
      message("Kuiper               Kuiper")
      message("Cramer-vonMises      CvM")
      message("Anderson-Darling     AD")
      message("Lehmann-Rosenblatt   LR")
      message("Zhang's tests        ZA")
      message("Wasserstein          Wassp1")
      message("Chi square tests     large, small")
    }
    TRUE
}
