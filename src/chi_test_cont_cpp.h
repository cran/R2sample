#ifndef CHI_TEST_CONT_CPP_H
#define CHI_TEST_CONT_CPP_H

#include <Rcpp.h>
Rcpp::List chi_test_cont_cpp(Rcpp::List dta, Rcpp::IntegerVector nbins=Rcpp::IntegerVector::create(100,10));

#endif
