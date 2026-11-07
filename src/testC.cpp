#include <Rcpp.h>
#include "calcTS.h"
#include "gen_sim_data.h"
using namespace Rcpp;

//' run test using either simulation or permutation.
//' 
//' @param dta a list with the data
//' @param TS routine to calculate test statistics for non-chi-square tests 
//' @param typeTS type of a test statistic
//' @param TSextra additional info passed to TS, if necessary
//' @param B =5000, number of simulation runs.
//' @keywords internal
//' @return A list with test statistics and p values
//' @noRd
// [[Rcpp::export(rng=false)]]
List testC(List dta, 
          Function TS,
          int typeTS,
          List TSextra,
          int B=5000) { 
               
  NumericVector TS_data;
  TS_data = calcTS(dta, TS, typeTS, TSextra);
  int const nummethods=TS_data.size();
  int i, j;
  NumericVector pvals(nummethods), TS_sim(nummethods), stats(nummethods);

  CharacterVector TSmethods=TS_data.names();
  stats.names() = TSmethods; 
  pvals.names() = TSmethods;
  for(i=0;i<nummethods;++i) {
    stats(i)=TS_data(i);  
    pvals(i)=0.0;  
  }
/*  if B=0 return just the test statistics */    
  if(B==0) return List::create(Named("statistics")=stats);
/* find test statistic for simulated or simuted data  */ 
  for(j=0;j<B;++j) {
    GetRNGstate();
    List simdta=gen_sim_data(dta, TSextra);
    PutRNGstate();
    TS_sim= calcTS(simdta, TS, typeTS, TSextra);
/* count simulated values at least as large as the observed one. Ties are
   common (KS, Kuiper, discrete data) and have to be counted, otherwise the
   test is anti-conservative. A small relative tolerance makes sure values
   that are equal up to rounding errors are treated as ties. */
    for(i=0;i<nummethods;++i) {
       double tol=1e-9*std::abs(TS_data(i));
       if(TS_sim(i)>=TS_data(i)-tol) pvals(i)=pvals(i)+1.0;
    }
  }

/* find p values and return list with test statistics and p values */  
  for(i=0;i<nummethods;++i) pvals(i)=pvals(i)/B;
  return  List::create(Named("statistics")=stats, Named("p.values")=pvals);
}

