#include <Rcpp.h>
#include "chi_test_disc_cpp.h"
#include "weights_cpp.h"
#include "permute_disc_cpp.h"
#include "TS_disc_cpp.h"
#include "chi_test_cont_cpp.h"
#include "permute_cont_cpp.h"
#include "TS_cont_cpp.h"
using namespace Rcpp;

//' run permutation test.
//' 
//' @param x A numeric vector.
//' @param y A numeric vector.
//' @param vals A numeric vector. Indicates discrete data.
//' @param nbins Two bin numbers for chi square  test.
//' @param B Number of simulation runs.
//' @param doMethod A character vector of methods to include
//' @return A list with test statistics and p values
// [[Rcpp::export]]
List perm_test_cpp(NumericVector x, 
                            NumericVector y, 
                            NumericVector vals=0.0, 
                            IntegerVector nbins= IntegerVector::create(100,10), 
                            int B=5000, 
                CharacterVector doMethod= CharacterVector::create("chi large", "chi small", "t test", "KS", "Kuiper", "CvM", "AD", "LR", "ZA", "ZK", "ZC", "Wassp1")) { 
                
  int const nummethods=10;
  int k=vals.size(), i, j;
  NumericVector tmp1(2), pvals(nummethods+2), chi_TS(2), perm_TS(nummethods), 
         data_TS(nummethods), stats(nummethods+2), adw(k);
  List dta=List::create(Named("x") = x, Named("y") = y, Named("vals")=vals);
  List tmp;
  
  stats.names() = doMethod.names(); 
  pvals.names() = doMethod.names();
   
/* if they are in list of methods, run chi square test routines*/   
  LogicalVector H=in(CharacterVector::create("chi large", "chi small"), doMethod);
  if( (H[0]==TRUE) || (H[1]==TRUE) ) {   
    if(k==1) tmp=chi_test_cont_cpp(dta, nbins);
    else tmp=chi_test_disc_cpp(dta, nbins);
    chi_TS=tmp[0];
    stats(0)=chi_TS(0);
    stats(1)=chi_TS(1);
    tmp1=tmp[1];
    pvals(0)=tmp1[0];
    pvals(1)=tmp1[1];
  }
/*  Find weights and test statistics for data */      
  if(k==1) data_TS=TS_cont_cpp(dta, doMethod);
  else {
     adw=weights_cpp(dta);
     data_TS=TS_disc_cpp(dta, adw, doMethod);  
  }
  
  for(i=2;i<nummethods+2;++i) {
    stats(i)=data_TS(i-2);  
    pvals(i)=0.0;  
  }

/*  if B=0 return just the test statistics */    
  if(B==0) return List::create(Named("statistics")=stats);
  
/* run permutation test  */    
  for(j=0;j<B;++j) {
     if(k==1) perm_TS=TS_cont_cpp(permute_cont_cpp(dta), doMethod);
     else perm_TS=TS_disc_cpp(permute_disc_cpp(dta), adw, doMethod);
     for(i=0;i<nummethods;++i) {
         if(data_TS(i)<perm_TS(i)) pvals(i+2)=pvals(i+2)+1.0;
     }
  }

/* find p values and return list with test statistics and p values */  
  for(i=2;i<nummethods+2;++i) pvals(i)=pvals(i)/B;
  return  List::create(Named("statistics")=stats, Named("p.values")=pvals);
  
}

