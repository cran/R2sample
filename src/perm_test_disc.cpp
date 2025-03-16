#include <Rcpp.h>
#include "permute_disc.h"
#include "teststatistics.h"
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
List perm_test_disc(IntegerVector x, 
                   IntegerVector y, 
                   NumericVector vals,
                   Function TS, 
                   int typeTS,
                   List TSextra, 
                   int samplingmethod=1, 
                   int B=5000) { 
  
  int i, j;   
/* find out what methods should be run */  
  NumericVector TS_data, TS_perm;
  List dta=List::create(Named("x") = x, 
                        Named("y") = y, 
                        Named("vals")=vals);
  TS_data = calcTS(dta, TS, typeTS, TSextra);
  int const nummethods=TS_data.size();
  NumericVector pvals(nummethods), stats(nummethods);
  CharacterVector TSmethods=TS_data.names();
  stats.names() = TSmethods; 
  pvals.names() = TSmethods;  

  /*  Find weights and test statistics for data */      

  for(i=0;i<nummethods;++i) {
    stats(i)=TS_data(i);  
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
     TS_perm = calcTS(dta, TS, typeTS, TSextra);
   /*  p value */
     for(i=0;i<nummethods;++i) {
         if(TS_data(i)<TS_perm(i)) pvals(i)=pvals(i)+1.0;
     }
  }

/* find p values and return list with test statistics and p values */  
  for(i=0;i<nummethods;++i) pvals(i)=pvals(i)/B;
  return  List::create(Named("statistics")=stats, Named("p.values")=pvals);
  
}

