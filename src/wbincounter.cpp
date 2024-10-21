#include <Rcpp.h>
using namespace Rcpp;
#include "Cpporder.h"

//' Find counts and/or sum of weights in bins. Useful for power calculations. Replaces hist command from R.
//' 
//' @param x numeric vector
//' @param bins numeric vector
//' @param w numeric vector of weights 
//' @keywords internal
//' @return sum of weights in bins
// [[Rcpp::export]]
Rcpp::NumericMatrix wbincounter(std::vector<double>& x, 
                                std::vector<double>& bins, 
                                std::vector<double>& w) {
  int n=x.size(), m=bins.size(), i, j;
  NumericMatrix wc(m-1, 2);
  
  w = Cpporder(w, x);
  std::sort(x.begin(), x.end());
  i=0;
  j=0;
  while ( (j<m-1) && (i<n) ) {
    if(x[i]<=bins[j+1]) {
       wc(j, 0) = wc(j, 0) + w[i];
       wc(j, 1) = wc(j, 1) + w[i]*w[i];
       ++i;
    }  
    else ++j;
    
  }
  return wc;   
}
