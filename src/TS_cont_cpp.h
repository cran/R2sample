#ifndef TS_CONT_CPP_H
#define TS_CONT_CPP_H

#include <Rcpp.h>
Rcpp::NumericVector TS_cont_cpp(Rcpp::List dta,
          Rcpp::CharacterVector doMethod= 
          Rcpp::CharacterVector::create("chi large", "chi small", 
          "t test", "KS", "Kuiper", "CvM", "AD", "LR", "ZA", "ZK", "ZC", "Wassp1"));

#endif
