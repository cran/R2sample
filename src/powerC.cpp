#include <Rcpp.h>
#include "weights.h"
#include "calcTS.h"
#include "gen_sim_data.h"
using namespace Rcpp;

//' Find the power of various continuous tests via simutation or permutation.
//' 
//' @param rxy a function that generates x and y data.
//' @param TS routine to calculate test statistics for non-chi-square tests
//' @param xparam  arguments for r1.
//' @param yparam  arguments for r2.
//' @param typeTS indicator for type of test statistics
//' @param TSextra additional info passed to TS, if necessary
//' @param B =1000 number of simulation runs
//' @keywords internal
//' @return A list values of test statistics
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
  NumericMatrix realdta(B*nl, 1+nummethods), simdta(B*nl, 1+nummethods);
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
      realdta(cn,0)=xparam(l);
      for(j=0;j<nummethods;++j) realdta(cn,j+1)=TS_data(j);
      List tmp=gen_sim_data(dta, TSextra);
      TS_sim = calcTS(tmp, TS, typeTS, TSextra);
      simdta(cn,0)=xparam(l);
      for(j=0;j<nummethods;++j) simdta(cn,j+1)=TS_sim(j);
    }  
  }
  return List::create(Named("Data")=realdta, 
               Named("Simulated")=simdta);
}

