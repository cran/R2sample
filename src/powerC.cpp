#include <Rcpp.h>
#include "weights.h"
#include "calcTS.h"
#include "gen_sim_data.h"
using namespace Rcpp;

//' Find the power of various continuous tests via simutation or permutation.
//' 
//' @param rxy a function that generates x and y data.
//' @param TS routine to calculate test statistics for non-chi-square tests
//' @param xparam  first argument for rxy.
//' @param yparam  second argument for rxy.
//' @param typeTS indicator for type of test statistics
//' @param TSextra additional info passed to TS, if necessary
//' @param B =1000 number of simulation runs
//' @keywords internal
//' @return A list values of test statistics
//' @noRd
// [[Rcpp::export(rng=false)]]
List powerC(Function rxy, 
            NumericVector xparam,                       
            NumericVector yparam,
            Function TS,
            int typeTS,
            List TSextra,
            int B=1000) { 
/* Find out how many tests are to be done, sample sizes etc. */  
  int j,l, m, nl=xparam.size();
  NumericVector TS_data;
  List dta=rxy(xparam(0),yparam(0));
  if(typeTS==5) {
    NumericVector adw=weights(dta);
    TSextra["adw"]=adw;
  }
  TS_data = calcTS(dta, TS, typeTS, TSextra);
  int const nummethods=TS_data.size();
  int cn=-1;
  NumericMatrix realdta(B*nl, nummethods), simdta(B*nl, nummethods), 
                paramalt(B*nl,2);
  NumericVector TS_sim(nummethods);
/*l loop over values in xparam */  
  for(l=0;l<nl;++l) {
/*m loop over simulation runs */  
    for(m=0;m<B;++m) {
      ++cn;
/*  create new data set  */            
      dta=rxy(xparam(l),yparam(l)); 
      if(typeTS==5) {
        NumericVector adw=weights(dta);
        TSextra["adw"]=adw;
      }
      TS_data = calcTS(dta, TS, typeTS, TSextra);
      paramalt(cn,0)=xparam(l);
      paramalt(cn,1)=yparam(l);
      for(j=0;j<nummethods;++j) realdta(cn,j)=TS_data(j);
      GetRNGstate();
      List tmp=gen_sim_data(dta, TSextra);
      PutRNGstate();
      TS_sim = calcTS(tmp, TS, typeTS, TSextra);
      for(j=0;j<nummethods;++j) simdta(cn,j)=TS_sim(j);
    }
  }
  return List::create(Named("Data")=realdta, 
               Named("Simulated")=simdta,
               Named("paramalt")=paramalt);
}

