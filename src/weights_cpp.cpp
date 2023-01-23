#include <Rcpp.h>
using namespace Rcpp;

//' find weights for several statistics for discrete data
//' 
//' @param dta A list with vectors x, y and vals.
//' @keywords internal
//' @return A vector of weights
// [[Rcpp::export]]

NumericVector weights_cpp(List dta)  {

  IntegerVector x = as<IntegerVector>(dta["x"]);
  IntegerVector y = as<IntegerVector>(dta["y"]);
  int k=x.size(), n, i, j, l, a, b;
  IntegerVector xy(k);
  NumericVector out(k);

  n=0;
  for(i=0;i<k;++i) {
     xy[i]=x[i]+y[i];
     n=n+xy[i];
  }  
  xy[k-1]=xy[k-1]-1;
  a=0;
  b=0;
  l=0;
  for(i=0;i<k;++i) {
      if(i==0) {a=1;b=xy[0];}
      else {
         a=a+xy[i-1];
         b=b+xy[i];
      }  
      for(j=a;j<=b;++j) {
        l=l+1;
        out(i)=out(i)+1.0/l/(n-l);
      }
  }
  return out;
}
