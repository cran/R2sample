#ifndef TSW_CONT_H
#define TSW_CONT_H

#include <Rcpp.h>
Rcpp::NumericVector TSw_cont(
   std::vector<double>& x, 
   std::vector<double>& y,
   std::vector<double>& wx, 
   std::vector<double>& wy);

#endif
