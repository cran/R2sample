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
   if(typeTS==1) TS_data=TS(dta["x"], dta["y"]); 
   if(typeTS==2) TS_data=TS(dta["x"], dta["y"], dta["wx"], dta["wy"]); 
   if(typeTS==3) TS_data=TS(dta["x"], dta["y"]);   
   if(typeTS==4) TS_data=TS(dta["x"], dta["y"], TSextra);
   if(typeTS==5) TS_data=TS(dta["x"], dta["y"], dta["vals"], TSextra["adw"]);  
   if(typeTS==6) TS_data=TS(dta["x"], dta["y"]);  
   if(typeTS==7) TS_data=TS(dta["x"], dta["y"], dta["vals"]);  
   if(typeTS==8) TS_data=TS(dta["x"], dta["y"], dta["vals"], TSextra); 
   return  TS_data;
 }
 
 
#include <Rcpp.h>
using namespace Rcpp;
//' This function calculates the test statistics for continuous data
//' @param  typeTS format of TS
//' @param  x continuous data set
//' @param  y continuous data set
//' @param  TS routine
//' @param  wx weights for x
//' @param  wy weights for y
//' @param  TSextra list passed to TS function
//' @keywords internal
//' @return A vector of numbers
// [[Rcpp::export]]
NumericVector ts_C(
       int typeTS,
       Rcpp::NumericVector x, 
       Rcpp::NumericVector y,
       Rcpp::Function TS, 
       Rcpp::List TSextra,
       Rcpp::NumericVector wx, 
       Rcpp::NumericVector wy) {
  NumericVector TS_data;
  if(typeTS==1) TS_data=TS(x, y); 
  if(typeTS==2) TS_data=TS(x, y, wx, wy); 
  if(typeTS==3) TS_data=TS(x, y);   
  if(typeTS==4) TS_data=TS(x, y, TSextra);
  return  TS_data;
}

//' This function calculates the test statistics for discrete data
//' @param  typeTS format of TS
//' @param  x discrete data set (counts)
//' @param  y discrete data set (counts)
//' @param  vals values of discrete RV
//' @param  TS routine 
//' @param  TSextra list passed to TS function
//' @param  adw vector of weights for Anderson-Darling test
//' @keywords internal
//' @return A vector of numbers
// [[Rcpp::export]]
NumericVector ts_D(
    int typeTS,
    Rcpp::IntegerVector x, 
    Rcpp::IntegerVector y,
    Rcpp::NumericVector vals,
    Rcpp::Function TS, 
    Rcpp::List TSextra,
    Rcpp::NumericVector adw) {
  NumericVector TS_data;
  if(typeTS==5)  TS_data=TS(x, y, vals, adw);  
  if(typeTS==6)  TS_data=TS(x, y);  
  if(typeTS==7)  TS_data=TS(x, y, vals);  
  if(typeTS==8)  TS_data=TS(x, y, vals, TSextra); 
  return  TS_data;
 }
 
 
