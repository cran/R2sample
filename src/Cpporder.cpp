#include <Rcpp.h>
using namespace Rcpp;

//' sort vector y by values in vector x
//' 
//' @param y numeric vector
//' @param x numeric vector
//' @keywords internal
//' @return numeric vector
// [[Rcpp::export]]  
std::vector<double> Cpporder(std::vector<double>& y, std::vector<double>& x) {
    int n=x.size();
    std::vector<int> Index(n);
    std::iota(Index.begin(), Index.end(), 0);
    std::sort(Index.begin(), Index.end(),
              [&](int A, int B) -> bool {
                return x[A] < x[B];
              });
    std::vector<double> out(n);
    for(int i=0;i<n;++i) out[i]=y[Index[i]];
    return out;
  }  
