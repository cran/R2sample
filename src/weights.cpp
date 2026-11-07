#include <Rcpp.h>
using namespace Rcpp;

//' find weights for several statistics for discrete data
//'
//' For class i, occupying positions a,..,b of the combined sorted sample,
//' the AD weight is sum (from a to b) 1/(l(N-l)), with l=N excluded.
//' Uses sum 1/(l(N-l)) = (sum 1/l + sum 1/(N-l))/N and harmonic numbers
//' H_m = digamma(m+1)+Euler's constant, so the cost is O(k), not O(N).
//'
//' @param dta A list with vectors x, y and vals.
//' @keywords internal
//' @return A vector of weights
//' @noRd
// [[Rcpp::export]]

NumericVector weights(List dta)  {

  NumericVector x = as<NumericVector>(dta["x"]);
  NumericVector y = as<NumericVector>(dta["y"]);
  int k=x.size(), i;
  double n=0.0, a, b;
  NumericVector out(k);

  for(i=0;i<k;++i) n=n+floor(x[i])+floor(y[i]);
  b=0.0;
  for(i=0;i<k;++i) {
      a=b+1.0;
      b=b+floor(x[i])+floor(y[i]);
      double bb = (i==k-1) ? b-1.0 : b;   /* exclude l=N */
      if(bb<a) continue;
      /* sum_{l=a}^{bb} 1/l = H_bb - H_{a-1};  sum 1/(N-l) = H_{N-a} - H_{N-bb-1} */
      out(i)=( R::digamma(bb+1.0)-R::digamma(a) +
               R::digamma(n-a+1.0)-R::digamma(n-bb) )/n;
  }
  return out;
}
