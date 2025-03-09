#include <Rcpp.h>
#include "weights.h"

using namespace Rcpp;

//' Find the power of various continuous tests via permutation.
//' 
//' @param rxy a function that generates x and y data.
//' @param TS routine to calculate test statistics for non-chi-square tests
//' @param xparam  arguments for r1.
//' @param yparam  arguments for r2.
//' @param typeTS indicator for type of test statistics
//' @param TSextra additional info passed to TS, if necessary
//' @param B =1000 number of simulation runs
//' @keywords internal
//' @return A list values of test statistics
// [[Rcpp::export]]
List power_cont(Function rxy, 
                NumericVector xparam,                       
                NumericVector yparam,
                Function TS,
                int typeTS,
                List TSextra,
                int B=1000) { 
/* Find out how many tests are to be done, sample sizes etc. */  
  List dta=rxy(xparam(0),yparam(0));
  NumericVector x = as<NumericVector>(dta["x"]);
  NumericVector y = as<NumericVector>(dta["y"]);
  NumericVector wx(x.length()), wy(y.length());
  if(dta.length()==4) {
    wx = as<NumericVector>(dta["wx"]);
    wy = as<NumericVector>(dta["wy"]);
  }  
  int nx=x.size(), ny=y.size(), n=nx+ny, nl=xparam.size();
  NumericVector TS_data;
  if(typeTS==1) TS_data=TS(x, y);
  if(typeTS==2) TS_data=TS(x, y,  wx, wy);
  if(typeTS==3) TS_data=TS(x, y);
  if(typeTS==4) TS_data=TS(x, y, TSextra);
  int const nummethods=TS_data.size();
  int i, j, m, l, cn=-1;
  NumericVector z(n), wz(n);
  NumericMatrix realdta(B*nl, 1+nummethods), permdta(B*nl, 1+nummethods);
  NumericVector TS_perm(nummethods);
  IntegerVector Index(n);

/*l loop over values in xparam */  
  for(l=0;l<nl;++l) {
/*m loop over simulation runs */  
    for(m=0;m<B;++m) {
      ++cn;
/*  create new data set  */            
      dta=rxy(xparam(l),yparam(l));
/*  find test statistics for data  */                
      x = as<NumericVector>(dta["x"]);
      y = as<NumericVector>(dta["y"]);
      for(i=0;i<nx;++i) {
        Index(i)=i;
        z(i)=x(i);
      }  
      for(i=0;i<ny;++i) {
        Index(i+nx)=i+nx;
        z(i+nx)=y(i);
      }
      if(dta.length()==4) {
        wx = as<NumericVector>(dta["wx"]);
        wy = as<NumericVector>(dta["wy"]);
      } 
      if(typeTS==1) TS_data=TS(x, y);
      if(typeTS==2) TS_data=TS(x, y,  wx, wy);
      if(typeTS==3) TS_data=TS(x, y);
      if(typeTS==4) TS_data=TS(x, y, TSextra);
      realdta(cn,0)=xparam(l);
      for(j=0;j<nummethods;++j) realdta(cn,j+1)=TS_data(j);
/*  find test statistics for permuted data   */ 
      Index = Rcpp::sample(Index, n);
      for(i=0;i<nx;++i) {
        x[i]=z[Index[i]];
        if(dta.length()==4) wx[i]=wz[Index[i]]; 
      }  
      for(i=0;i<ny;++i) {
        y[i]=z[Index[i+nx]];
        if(dta.length()==4) wy[i]=wz[Index[i+nx]];
      }
      if(typeTS==1) TS_perm=TS(x, y); 
      if(typeTS==2) TS_perm=TS(x, y, wx, wy); 
      if(typeTS==3) TS_perm=TS(x, y);   
      if(typeTS==4) TS_perm=TS(x, y, TSextra);
      permdta(cn,0)=xparam(l);
      for(j=0;j<nummethods;++j) permdta(cn,j+1)=TS_perm(j);
    }  
  }
  return List::create(Named("Data")=realdta, 
               Named("Permuted")=permdta);
}

