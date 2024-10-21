#ifndef WBINCOUNTER_H
#define WBINCOUNTER_H

#include <Rcpp.h>
Rcpp::NumericMatrix wbincounter(
      std::vector<double>& x, 
      std::vector<double>& bins,
      std::vector<double>& w);

#endif
