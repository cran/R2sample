#ifndef CALCTS_H
#define CALCTS_H

#include <Rcpp.h>
Rcpp::NumericVector calcTS(
       Rcpp::List dta, 
       Rcpp::Function TS, 
       int typeTS,
       Rcpp::List TSextra
        );

#endif

#ifndef TS_C_H
#define TS_C_H

#include <Rcpp.h>
Rcpp::NumericVector ts_C(
       int typeTS,
       Rcpp::NumericVector x, 
       Rcpp::NumericVector y,
       Rcpp::Function TS, 
       Rcpp::List TSextra,
       Rcpp::NumericVector wx, 
       Rcpp::NumericVector wy
        );

#endif

#ifndef TS_D_H
#define TS_D_H

#include <Rcpp.h>
Rcpp::NumericVector ts_D(
    int typeTS,
    Rcpp::IntegerVector x, 
    Rcpp::IntegerVector y,
    Rcpp::NumericVector vals,
    Rcpp::Function TS, 
    Rcpp::List TSextra,
    Rcpp::NumericVector adw
        );

#endif
