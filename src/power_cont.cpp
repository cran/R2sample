#include <Rcpp.h>
#include "weights.h"
#include "teststatistics.h"

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
  int nx=x.size(), ny=y.size(), n=nx+ny, nl=xparam.size();
  NumericVector TS_data;
  TS_data = calcTS(dta, TS, typeTS, TSextra);
  int const nummethods=TS_data.size();
  int i, j, m, l, cn=-1;
  NumericVector wx(nx), wy(ny), z(n), wz(n);
  NumericMatrix realdta(B*nl, 1+nummethods), permdta(B*nl, 1+nummethods);
  NumericVector TS_perm(nummethods);
  IntegerVector Index(n);

  bool withweights=false;
  if(typeTS==2) withweights=true;
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
      if(withweights) {
        wx = as<NumericVector>(dta["wx"]);
        wy = as<NumericVector>(dta["wy"]);
        
      }
      for(i=0;i<nx;++i) {
        Index(i)=i;
        z(i)=x(i);
        if(typeTS==2) wz(i)=wx(i);
      }  
      for(i=0;i<ny;++i) {
        Index(i+nx)=i+nx;
        z(i+nx)=y(i);
        if(withweights) wz(i+nx)=wy(i);
      }
      TS_data = calcTS(dta, TS, typeTS, TSextra);
      realdta(cn,0)=xparam(l);
      for(j=0;j<nummethods;++j) realdta(cn,j+1)=TS_data(j);
/*  find test statistics for permuted data   */ 
      Index = Rcpp::sample(Index, n);
      for(i=0;i<nx;++i) {
        x[i]=z[Index[i]];
        if(withweights) wx[i]=wz[Index[i]]; 
      }  
      for(i=0;i<ny;++i) {
        y[i]=z[Index[i+nx]];
        if(withweights) wy[i]=wz[Index[i+nx]];
      }
      dta["x"]=x;
      dta["y"]=y;
      if(withweights) {
        dta["wx"]=wx;
        dta["wy"]=wy;
      }
      TS_perm = calcTS(dta, TS, typeTS, TSextra);
      permdta(cn,0)=xparam(l);
      for(j=0;j<nummethods;++j) permdta(cn,j+1)=TS_perm(j);
    }  
  }
  return List::create(Named("Data")=realdta, 
               Named("Permuted")=permdta);
}

