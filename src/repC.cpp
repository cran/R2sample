#include <Rcpp.h>
using namespace Rcpp;

//' cpp version of R routine rep
//' 
//' @param x numeric vector
//' @param times integer vector
//' @keywords internal
//' @return A numeric vector
// [[Rcpp::export]]
Rcpp::NumericVector repC(Rcpp::NumericVector x, Rcpp::IntegerVector times) {
  int n=x.size(), ntimes, i, j, k;
  
  ntimes=times(0);
  for(i=1;i<n;++i) ntimes=ntimes+times(i);
  NumericVector xx(ntimes);

  k=0;
  for(i=0;i<n;++i) {
    if(times(i)>0) {
      for(j=0;j<times(i);++j) {
         xx(k)=x(i);
         ++k;
      }   
    }
  }
  return xx;   
}
