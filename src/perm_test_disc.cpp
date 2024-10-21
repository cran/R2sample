#include <Rcpp.h>
#include "weights.h"
#include "permute_disc.h"
using namespace Rcpp;

//' run permutation test.
//' 
//' @param x A vector of counts or weights.
//' @param y A vector of counts or weights.
//' @param vals A numeric vector. Indicates discrete data.
//' @param TS routine to calculate test statistics for non-chi-square tests 
//' @param typeTS type of a test statistic
//' @param TSextra additional info passed to TS, if necessary
//' @param samplingmethod =1, 1 for independence sampling, 2 for MCMC
//' @param B Number of simulation runs.
//' @keywords internal
//' @return A list with test statistics and p values
// [[Rcpp::export]]
List perm_test_disc(NumericVector x, 
                   NumericVector y, 
                   NumericVector vals,
                   Function TS, 
                   int typeTS,
                   List TSextra, 
                   int samplingmethod=1, 
                   int B=5000) { 
  
  int i, j;   
/* find out what methods should be run */  
  List dta=List::create(Named("x") = x, Named("y") = y, Named("vals")=vals); 
  NumericVector data_TS, perm_TS,  adw;
  if(typeTS==5) {
    adw=weights(dta);
    data_TS=TS(x, y, vals, adw); 
  }    
  if(typeTS==6)  data_TS=TS(x, y);  
  if(typeTS==7)  data_TS=TS(x, y, vals);  
  if(typeTS==8)  data_TS=TS(x, y, vals, TSextra); 
  int const nummethods=data_TS.size();
  NumericVector pvals(nummethods), stats(nummethods);
  CharacterVector TSmethods=data_TS.names();
  stats.names() = TSmethods; 
  pvals.names() = TSmethods;  

  /*  Find weights and test statistics for data */      

  for(i=0;i<nummethods;++i) {
    stats(i)=data_TS(i);  
    pvals(i)=0.0;  
  }
    
/*  if B=0 return just the test statistics */    
  if(B==0) {
    return List::create(Named("statistics")=stats);
  }  
/* run permutation test  */    
  for(j=0;j<B;++j) {
     GetRNGstate();
     dta=permute_disc(dta, samplingmethod);
     PutRNGstate();
     if(typeTS==5) perm_TS=TS(dta["x"], dta["y"], vals, adw);
     if(typeTS==6) perm_TS=TS(dta["x"], dta["y"]);     
     if(typeTS==7) perm_TS=TS(dta["x"], dta["y"], vals);
     if(typeTS==8) perm_TS=TS(dta["x"], dta["y"], vals, TSextra);
   /*  p value */
     for(i=0;i<nummethods;++i) {
         if(data_TS(i)<perm_TS(i)) pvals(i)=pvals(i)+1.0;
     }
  }

/* find p values and return list with test statistics and p values */  
  for(i=0;i<nummethods;++i) pvals(i)=pvals(i)/B;
  return  List::create(Named("statistics")=stats, Named("p.values")=pvals);
  
}

