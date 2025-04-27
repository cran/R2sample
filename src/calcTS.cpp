#include <Rcpp.h>
using namespace Rcpp;
//' This function calculates the test statistics for continuous data
//' @param  dta data set
//' @param  TS routine
//' @param  typeTS format of TS
//' @param  TSextra list passed to TS function
//' @keywords internal
//' @return A vector of numbers
// [[Rcpp::export]]
NumericVector calcTS(
     Rcpp::List dta, 
     Rcpp::Function TS,
     int typeTS,
     Rcpp::List TSextra) {
   NumericVector TS_data;
   if(typeTS<5) {
     NumericVector xx=dta["x"], x(xx.size()), yy=dta["y"], y(yy.size());
     int nx=xx.size(),ny=yy.size();
     for(int i=0;i<nx;++i) x[i]=xx[i];
     for(int i=0;i<ny;++i) y[i]=yy[i];
     if(typeTS==1) TS_data=TS(x, y); 
     if(typeTS==2) {
       NumericVector wxx=dta["x"], wx(wxx.size()), wyy=dta["y"], wy(wyy.size());
       for(int i=0;i<nx;++i) wx[i]=wxx[i];
       for(int i=0;i<ny;++i) wy[i]=wyy[i];       
       TS_data=TS(x, y, wx, wy); 
     }   
     if(typeTS==3) TS_data=TS(x, y);   
     if(typeTS==4) TS_data=TS(x, y, TSextra);
   }
   else {
     NumericVector vals=dta["vals"];
     IntegerVector xx=dta["x"], x(xx.size()), yy=dta["y"], y(yy.size());
     for(int i=0;i<vals.size();++i) {
       x[i]=xx[i];
       y[i]=yy[i];
     }   
     if(typeTS==5) {
       NumericVector adw=TSextra["adw"];
       TS_data=TS(x, y, vals, adw);  
     }   
     if(typeTS==6) TS_data=TS(x, y, vals);  
     if(typeTS==7) TS_data=TS(x, y, vals, TSextra); 
   }
   return  TS_data;
 }
 
 
