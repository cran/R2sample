#include <Rcpp.h>
using namespace Rcpp;

//' permute continuous data
//' 
//' @param dta A list of numeric vectors.
//' @return A list of permuted x and y vectors
// [[Rcpp::export]]
List permute_cont_cpp(List dta) {
  
  NumericVector x = as<NumericVector>(dta["x"]);
  NumericVector y = as<NumericVector>(dta["y"]);
  
  int nx=x.size(), ny=y.size(), n;
  NumericVector z(nx+ny), x1(nx), y1(ny);
  n=nx+ny;
  for(int i=0;i<nx;++i) z[i]=x[i];
  for(int i=0;i<ny;++i) z[i+nx]=y[i];
  z = sample(z, n);
  for(int i=0;i<nx;++i) x1[i]=z[i];
  for(int i=0;i<ny;++i) y1[i]=z[i+nx];
  
  return List::create(Named("x") = x1, 
                      Named("y") = y1);
}
