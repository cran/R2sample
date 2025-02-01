#include <Rcpp.h>
using namespace Rcpp;

//' a local function needed for the vignette
//' 
//' @param x An integer vector.
//' @param y An integer vector.
//' @param vals A numeric vector with the values of the discrete rv.
//' @return A vector with test statistics
//' @export
// [[Rcpp::export]]
NumericVector myTS2(IntegerVector x, 
                    IntegerVector y,
                    NumericVector vals) {
    
    Rcpp::CharacterVector methods=CharacterVector::create("std t test");    
  int const nummethods=methods.size();
  int k=x.size(), n=sum(x), i;
  double meanx=0.0, meany=0.0, sdx=0.0, sdy=0.0;
  NumericVector TS(nummethods);
  TS.names() =  methods;
  for(i=0;i<k;++i) {
    meanx = meanx + x[i]*vals[i]/n;
    meany = meany + y[i]*vals[i]/n;
  }
  for(i=0;i<k;++i) {
    sdx = sdx + x[i]*(vals[i] - meanx)*(vals[i] - meanx);
    sdy = sdy + y[i]*(vals[i] - meany)*(vals[i] - meany);
  }  
  sdx = sqrt(sdx/(n-1));
  sdy = sqrt(sdy/(n-1));
  TS(0) = std::abs(meanx/sdx - meany/sdy);
  return TS;
}
