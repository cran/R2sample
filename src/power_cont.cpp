#include <Rcpp.h>
#include "weights.h"

using namespace Rcpp;

//' Find the power of various discrete tests via permutation.
//' 
//' @param rxy a function that generates x and y data.
//' @param TS routine to calculate test statistics for non-chi-square tests
//' @param typeTS indicator for type of test statistics
//' @param TSextra additional info passed to TS, if necessary
//' @param alpha A numeric constant
//' @param B =c(1000, 1000) number of simulation runs for power and for p values
//' @param xparam  arguments for r1.
//' @param yparam  arguments for r2.
//' @keywords internal
//' @return A numeric matrix of powers
// [[Rcpp::export]]
NumericMatrix power_cont(Function rxy, 
                        Function TS,
                        int typeTS,
                        List TSextra,
                        double alpha=0.05, 
                        NumericVector B=NumericVector::create(1000, 1000), 
                        NumericVector xparam=0.0,                       
                        NumericVector yparam=0.0) { 
/* Find out how many tests are to be done, sample sizes etc. */  
  List dta=rxy(xparam(0),yparam(0));
  NumericVector x = as<NumericVector>(dta["x"]);
  NumericVector y = as<NumericVector>(dta["y"]);
  NumericVector wx(x.length()), wy(y.length());
  if(dta.length()==4) {
    wx = as<NumericVector>(dta["wx"]);
    wy = as<NumericVector>(dta["wy"]);
  }  
  int nx=x.size(), ny=y.size(), n=nx+ny;
  NumericVector TS_data;
  if(typeTS==1) TS_data=TS(x, y);
  if(typeTS==2) TS_data=TS(x, y,  wx, wy);
  if(typeTS==3) TS_data=TS(x, y);
  if(typeTS==4) TS_data=TS(x, y, TSextra);
  int const nummethods=TS_data.size();
  int rp=yparam.size(), i, j, l, k, m;
  NumericVector z(n), wz(n);
  NumericMatrix pwr(rp, nummethods);
  NumericVector TS_perm(nummethods);
  IntegerVector counter(nummethods), Index(nx+ny);
  
  
/*  l loop over values in xparam, yparam */    
  for(l=0;l<rp;++l) {
    for(j=0;j<nummethods;++j) pwr(l,j)=0.0;  
/*  m loop over simulation runs */  

    for(m=0;m<B[0];++m) {
/*    create new data set  */            
      dta=rxy(xparam(l),yparam(l));
/*       find test statistics for data  */                
      x = as<NumericVector>(dta["x"]);
      y = as<NumericVector>(dta["y"]);
      if(dta.length()==4) {
        wx = as<NumericVector>(dta["wx"]);
        wy = as<NumericVector>(dta["wy"]);
      } 
      for(i=0;i<nx;++i) {
        Index(i)=i;
        z[i]=x(i);
        wz[i]=wx[i];
      }  
      for(i=0;i<ny;++i) {
        Index(i+nx)=i+nx;
        z[i+nx]=y(i);
        wz[i+nx]=wy[i];
      }
      if(typeTS==1) TS_data=TS(x, y);
      if(typeTS==2) TS_data=TS(x, y,  wx, wy);
      if(typeTS==3) TS_data=TS(x, y);
      if(typeTS==4) TS_data=TS(x, y, TSextra);
      for(j=0;j<nummethods;++j) counter[j]=0;
      for(k=0;k<B[1];++k) {     
/*    k loop to find test statistics for permuted data   */ 
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
      /* and compare with test statistics from data */       
          for(j=0;j<nummethods;++j) {
             if(TS_perm(j)>TS_data(j)) counter[j]=counter[j]+1;
          }
       } /* end of k loop of p value simulation*/ 
/*    p values < alpha?  */                 
      for(j=0;j<nummethods;++j) 
        if(double(counter[j])/B[1]<alpha) 
            pwr(l, j)=pwr(l, j)+1.0;
    }  /* end of loop of power simulation*/ 

    for(j=0;j<nummethods;++j) pwr(l, j)=pwr(l, j)/B[0];
  }  /* end of l loop over cases */
  return  pwr;
  
}

