#ifndef TS_DISC_H
#define TS_DISC_H

#include <Rcpp.h>
Rcpp::NumericVector TS_disc(
           Rcpp::IntegerVector x,
           Rcpp::IntegerVector y, 
           Rcpp::NumericVector vals,
           Rcpp::NumericVector ADweights=Rcpp::NumericVector::create(2));

#endif
