#include <Rcpp.h>
#include <cmath>

using namespace Rcpp; 

//' permute discrete data
//' 
//' @param dta A list of numeric vectors.
//' @return A list of permuted vectors
// [[Rcpp::export]]
List permute_disc_cpp(List dta) {

  NumericVector x = as<NumericVector>(dta["x"]);
  NumericVector y = as<NumericVector>(dta["y"]);
  NumericVector vals = as<NumericVector>(dta["vals"]); 
  int nx, ny, n, tmp, xmax, xsum, i, j, k=vals.size();
  
  IntegerVector xy(k),  xnew(k), ynew(k);
  double p;
  nx=0;
  ny=0;
  xmax=0.0;
  xsum=0;  
  for(i=0;i<k;++i) {
    nx+=x[i];
    ny+=y[i];
  }     
  n=nx+ny;
  p=double(nx)/n;
  j=0;
  for(i=0;i<k;++i) {
     xy[i]=x[i]+y[i];
     xsum=xsum+x[i];
     if(x[i]>xmax) {
        xmax=x[i];
        j=i;
     }  
  }
  do {
    for(i=0;i<k;++i) {
      xnew[i] = rbinom(1, xy[i], p)[0];
      ynew[i] = xy[i]-xnew[i];
    }
    tmp=0;
    xnew[j]=0;
    for(i=0;i<k;++i) tmp=tmp+xnew[i];
    xnew[j]=xsum-tmp;
    ynew[j]=xy[j]-xnew[j];
  } while ((xnew[j]<0)|(ynew[j]<0));   
  return List::create(Named("vals") = vals,
                      Named("x") = xnew, 
                      Named("y") = ynew);                    
}
