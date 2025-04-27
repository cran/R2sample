#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
void f() {}

// [[Rcpp::export]]
void test(Function g) {
  for(int j=0;j<5;++j) {
    Rcout<<R::runif(0,1)<<" ";

  }
  Rcout<<"\n";
  for(int j=0;j<5;++j) {
    g();
    GetRNGstate();
    Rcout<<R::runif(0,1)<<" ";
    PutRNGstate();
    
  }

}

