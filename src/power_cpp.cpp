#include <Rcpp.h>
#include "chi_test_cont_cpp.h"
#include "weights_cpp.h"
#include "permute_cont_cpp.h"
#include "TS_cont_cpp.h"
#include "chi_test_disc_cpp.h"
#include "permute_disc_cpp.h"
#include "TS_disc_cpp.h"
#include "bincounter_cpp.h"
using namespace Rcpp;


//' Find the power of various tests via permutation.
//' 
//' @param rxy a function that generates x and y data.
//' @param nbins Two bin numbers for chi square  test.
//' @param alpha A numeric constant
//' @param B Number of simulation runs.
//' @param xparam  arguments for r1.
//' @param yparam  arguments for r2.
//' @param doMethod A character vector of methods to include
//' @return A numeric matrix of powers
// [[Rcpp::export]]
NumericMatrix power_cpp(Function rxy, 
                        IntegerVector nbins= IntegerVector::create(100,10), 
                        double alpha=0.05, 
                        int B=1000, 
                        NumericVector xparam=0.0,                       
                        NumericVector yparam=0.0, 
                CharacterVector doMethod= CharacterVector::create("chi large", 
       "chi small", "t test", "KS", "Kuiper", "CvM", "AD", "LR", "ZA", "ZK", "ZC", "Wassp1")) { 
 
  int const nummethods=10;                      
  int rp=yparam.size(), i, j, l, k;
  NumericMatrix pwr(rp, nummethods+2),zeromatrix(rp, nummethods+2);
  NumericVector sim_data(nummethods),sim_perm(nummethods), crit_vals(nummethods), tmp(2);
  IntegerVector counter(nummethods);
  colnames(pwr) = CharacterVector::create("chi large", "chi small", "t test", "KS", "Kuiper", 
                                              "CvM", "AD", "LR", "ZA", "ZK", "ZC", "Wassp1");
  CharacterVector rna(rp);
  NumericVector w(1000);
/*   include chi square tests?  */          
  LogicalVector H=in(CharacterVector::create("chi large", "chi small"), doMethod);
  for(l=0;l<rp;++l) rna(l)=std::to_string(yparam[l]);
  rownames(pwr) =  rna; 
  
/*  loop over values in xparam, yparam */    
  for(l=0;l<rp;++l) {
/*  loop over simulation runs */  
    for(i=0;i<B;++i) {
/*    create new data set  */            
      List dta=rxy(xparam(l),yparam(l));
      
/*    if data is discrete check that vectors are of equal length  
      and either x or y is greater than 0 for each value of vals. 
      If not write message and return matrix of zeros    
      dta.size()==3 means discrete data  */                    
      if(dta.size()==3) {
          NumericVector vals = as<NumericVector>(dta["vals"]);    
          IntegerVector x = as<IntegerVector>(dta["x"]);
          IntegerVector y = as<IntegerVector>(dta["y"]);
          if( (x.size()!=vals.size()) | (y.size()!=vals.size()) ) {
               Rcout<<"Data generated has x, y and vals vectors of unequal lengths. Check your rxy function!\n";
               return zeromatrix;            
          }
          for(int i0=0;i0<vals.size();++i0) {
              if(x[i0]+y[i0]==0) {
                  Rcout<<"Data generated has x=0 and y=0 for some value of vals. Check your rxy function!\n";
                  return zeromatrix;                          
              }
          }
      }
/*       find test statistics for data  */              
/*    dta.size()==2 means continuous data  */              
      if(dta.size()==2) {        
         sim_data=TS_cont_cpp(dta, doMethod);
         if( (H[0]==TRUE) | (H[1]==TRUE) )  tmp=chi_test_cont_cpp(dta, nbins)[1];
      }  
      else {
         w=weights_cpp(dta);
         sim_data=TS_disc_cpp(dta, w, doMethod);
         if( (H[0]==TRUE) | (H[1]==TRUE) )  tmp=chi_test_disc_cpp(dta, nbins)[1];
      }
      for(j=0;j<nummethods;++j) counter[j]=0;
      for(k=0;k<B;++k) {     
/*    find test statistics for permuted data   */                       
         if(dta.size()==2) sim_perm=TS_cont_cpp(permute_cont_cpp(dta), doMethod);
         else sim_perm=TS_disc_cpp(permute_disc_cpp(dta), w, doMethod);  
/* and compare with test statistics from data */       
         for(j=0;j<nummethods;++j) {
            if(sim_perm(j)>sim_data(j)) counter[j]=counter[j]+1;
         }
         
      } /* end of loop of p value simulation*/ 
/*    p values of chi square tests < alpha?  */           
      if(tmp(0)<alpha) pwr(l, 0)=pwr(l, 0)+1.0;
      if(tmp(1)<alpha) pwr(l, 1)=pwr(l, 1)+1.0;
/*    p values of other tests < alpha?  */                 
      for(j=0;j<nummethods;++j) 
        if(double(counter[j])/B<alpha) 
            pwr(l, j+2)=pwr(l, j+2)+1.0;
    }  /* end of loop of power simulation*/ 

    for(j=0;j<nummethods+2;++j) pwr(l, j)=pwr(l, j)/B;
  }  /* end of loop over cases */
  return  pwr;
  
}

