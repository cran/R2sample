#ifndef GEN_SIM_DATA_H
#define GEN_SIM_DATA_H

#include <Rcpp.h>
Rcpp::List gen_sim_data(Rcpp::List dta, Rcpp::List TSextra);

#endif

#ifndef GEN_CONT_NOWEIGHTS_H
#define GEN_CONT_NOWEIGHTS_H

#include <Rcpp.h>
Rcpp::List gen_cont_noweights(
              Rcpp::NumericVector x,
              Rcpp::NumericVector y, 
              Rcpp::List TSextra);
#endif

#ifndef GEN_CONT_WEIGHTS_H
#define GEN_CONT_WEIGHTS_H

#include <Rcpp.h>
Rcpp::List gen_cont_weights(
              Rcpp::NumericVector x,
              Rcpp::NumericVector y, 
              Rcpp::NumericVector wx,
              Rcpp::NumericVector wy, 
              Rcpp::List TSextra);

#endif

#ifndef GEN_DISC_H
#define GEN_DISC_H

#include <Rcpp.h>
Rcpp::List gen_disc(
              Rcpp::NumericVector x,
              Rcpp::NumericVector y, 
              Rcpp::NumericVector vals,
              Rcpp::List TSextra);

#endif
