#include <Rcpp.h>
#include "weights.h"
#include "permute_disc.h"
#include "teststatistics.h"
using namespace Rcpp;

//' Find the power of various continuous tests via permutation.
//' 
//' @param rxy a function that generates x and y data.
//' @param TS routine to calculate test statistics for non-chi-square tests
//' @param xparam  arguments for r1.
//' @param yparam  arguments for r2.
//' @param typeTS indicator for type of test statistics
//' @param TSextra additional info passed to TS, if necessary
//' @param samplingmethod =1, 1 for independence sampling, 2 for MCMC
//' @param B =1000 number of simulation runs
//' @keywords internal
//' @return A list values of test statistics
// [[Rcpp::export]]
List power_disc(Function rxy, 
                NumericVector xparam,                       
                NumericVector yparam,
                Function TS,
                int typeTS,
                List TSextra,
                int samplingmethod=1,
                int B=1000) { 
/* Find out how many tests are to be done, sample sizes etc. */  
  List dta=rxy(xparam(0),yparam(0));
  int nl=xparam.size();
  NumericVector TS_data;
  if(typeTS==5) TSextra["adw"]=weights(dta);
  TS_data = calcTS(dta, TS, typeTS, TSextra);
  int const nummethods=TS_data.size();
  int j, m, l, cn=-1;
  NumericMatrix realdta(B*nl, 1+nummethods), permdta(B*nl, 1+nummethods);
  NumericVector TS_perm(nummethods);

/*l loop over values in xparam */  
  for(l=0;l<nl;++l) {
/*m loop over simulation runs */  
    for(m=0;m<B;++m) {
      ++cn;
/*  create new data set  */            
      dta=rxy(xparam(l),yparam(l));
/*  find test statistics for data  */
      if(typeTS==5) TSextra["adw"]=weights(dta);
      TS_data = calcTS(dta, TS, typeTS, TSextra);
      realdta(cn,0)=xparam(l);
      for(j=0;j<nummethods;++j) realdta(cn,j+1)=TS_data(j);

/*  find test statistics for permuted data   */ 
      GetRNGstate();
      dta=permute_disc(dta, samplingmethod);
      PutRNGstate();
      TS_perm = calcTS(dta, TS, typeTS, TSextra);
      permdta(cn,0)=xparam(l);
      for(j=0;j<nummethods;++j) permdta(cn,j+1)=TS_perm(j);
    }  
  }
  return List::create(Named("Data")=realdta, 
               Named("Permuted")=permdta);
}

