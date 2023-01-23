#include <Rcpp.h>
using namespace Rcpp;

//' run chi square test for discrete data.
//' 
//' @param dta A list of numeric vectors.
//' @param nbins Integer vector of length 2 with number of bins.
//' @keywords internal
//' @return A list with test statistics, p values and degrees of freedom
// [[Rcpp::export]]
Rcpp::List  chi_test_disc_cpp(Rcpp::List dta, Rcpp::IntegerVector nbins=Rcpp::IntegerVector::create(100, 10)) {
  IntegerVector x = as<IntegerVector>(dta["x"]);
  IntegerVector y = as<IntegerVector>(dta["y"]);
  NumericVector chi(2), pvals(2);
  int i, l, nx, ny, n;
  double pE;

/*  find sample sizes */  
  nx=0;
  ny=0;
  for(i=0;i<x.size();++i) {
    nx+=x[i];
    ny+=y[i];
  }
  n=nx+ny;  
  
/* Make sure no bins have combined counts less than 5 */  
  while (min(x+y)<5) {
    for(i=0;i<x.size();++i) {
      if(x[i]+y[i]<5) {
         if(i==x.size()-1) i=i-1;
         x[i]=x[i]+x[i+1];
         x.erase(i+1);
         y[i]=y[i]+y[i+1];
         y.erase(i+1);
         break;
      }
    }
  }  
/* chi large uses original binning/classes */  
  nbins[0]=x.size(); 
  
/* combine bins if needed and do chi-square  tests */
  for(l=0;l<2;++l) {
    while(x.size()>nbins[l]) {
       i=which_min(x+y);
       if(i==x.size()-1) i=x.size()-2;
       x[i]=x[i]+x[i+1];
       x.erase(i+1);
       y[i]=y[i]+y[i+1];
       y.erase(i+1);     
    }

    chi(l)=0.0;
    for(i=0;i<x.size();++i) {
      pE=(x(i)+y(i))/double(n);
      chi(l)=chi(l)+(x(i)-nx*pE)*(x(i)-nx*pE)/nx/pE;
      chi(l)=chi(l)+(y(i)-ny*pE)*(y(i)-ny*pE)/ny/pE; 
    }
    pvals(l)=1-R::pchisq(chi(l), x.size(), 1, 0);
    nbins[l]=x.size();
  }

  return List::create(Named("statistics")=chi, Named("p.values")=pvals, Named("df")=nbins);
}
