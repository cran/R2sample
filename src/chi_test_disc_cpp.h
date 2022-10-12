#ifndef CHI_TEST_DISC_CPP_H
#define CHI_TEST_DISC_CPP_H

#include <Rcpp.h>
Rcpp::List chi_test_disc_cpp(Rcpp::List dta, Rcpp::IntegerVector nbins=Rcpp::IntegerVector::create(100,10));

#endif
