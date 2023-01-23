#include <Rcpp.h>
using namespace Rcpp;

//' run chi square test for continuous data.
//' 
//' @param dta A list of numeric vectors.
//' @param nbins A vector of length 2 of bin lengths.
//' @keywords internal
//' @return A list with test statistics and p values
// [[Rcpp::export]]
Rcpp::List chi_test_cont_cpp(Rcpp::List dta, Rcpp::IntegerVector nbins= IntegerVector::create(100,10)) {
  NumericVector x = as<NumericVector>(dta["x"]);
  NumericVector y = as<NumericVector>(dta["y"]);  
  int nx=x.size(),ny=y.size(),n=nx+ny, i, j, k;
  NumericVector xy(n), bins(nbins(0)+1), chi(2), pvals(2);
  IntegerVector df(2), xc(nbins(0)), yc(nbins(0));
  double pE;


/* for small data sets choose nbins(0) so that the combined 
     data set has expected counts at least five in each bin */
                 
  if(x.size()+y.size()<500) {
     k = (x.size()+y.size())/5.0;
     if(nbins(0)>k) nbins(0) = k;                 
  }
  if(nbins(0)<nbins(1)) nbins(1)=nbins(0);
/* create combined data set xy and sort vectors */  
  for(i=0;i<nx;++i) xy[i]=x[i];
  for(i=0;i<ny;++i) xy[i+nx]=y[i];
  std::sort(xy.begin(), xy.end());
  std::sort(x.begin(), x.end());
  std::sort(y.begin(), y.end());
  
/* for each number of bins, find equal probability bins */  
  for(k=0;k<2;++k) {
    bins(0)=xy(0)-0.000001;
    for(i=1;i<nbins(k);++i) bins(i)=xy(i*n/nbins(k));
    bins(nbins(k))=xy(n-1)+0.000001;
/* find x counts  */      
    i=0;
    j=0;
    while (j<nbins(k))  {
      xc(j)=0;
      while (x(i)<bins(j+1)) {
        xc(j)=xc(j)+1;
        ++i;
        if(i==nx) break;
      }
      ++j;
      if(i==nx) break;
    }
/* find y counts  */      
    i=0;
    j=0;
    while (j<nbins(k))  {
      yc(j)=0;
      while (y(i)<bins(j+1)) {
        yc(j)=yc(j)+1;
        ++i;
        if(i==ny) break;
      }   
      ++j;
      if(i==ny) break;
    }
/*  find test statistics */      
    chi(k)=0.0;
    df(k)=0;
    for(i=0;i<nbins(k);++i) {
       if(xc(i)+yc(i)>0) {
         df(k)=df(k)+1;
         pE=(xc(i)+yc(i))/double(n);     
         chi(k)=chi(k)+(xc(i)-nx*pE)*(xc(i)-nx*pE)/nx/pE;
         chi(k)=chi(k)+(yc(i)-ny*pE)*(yc(i)-ny*pE)/ny/pE; 
       }  
    }

/* find p values using chi-square approximation  */      
    pvals(k)=1-R::pchisq(chi(k), df(k), 1, 0);
  }  
  return List::create(Named("statistics")=chi, Named("p.value")=pvals,  Named("df")=df);
}
