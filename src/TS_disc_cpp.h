#ifndef TS_DISC_CPP_H
#define TS_DISC_CPP_H

#include <Rcpp.h>
Rcpp::NumericVector TS_disc_cpp(Rcpp::List dta, 
           Rcpp::NumericVector ADweights, 
           Rcpp::CharacterVector doMethod= Rcpp::CharacterVector::create("chi large", 
           "chi small", "t test", "KS", "Kuiper", "CvM", "AD", "LR", 
           "ZA", "ZK", "ZC", "Wassp1"));

#endif
