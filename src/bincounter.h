#ifndef BINCOUNTER_H
#define BINCOUNTER_H

#include <Rcpp.h>
Rcpp::IntegerVector bincounter(
      Rcpp::NumericVector x, 
      Rcpp::NumericVector bins);

#endif
