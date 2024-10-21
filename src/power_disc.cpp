#include <Rcpp.h>
#include "permute_disc.h"
#include "weights.h"
using namespace Rcpp;

//' Find the power of various tests via permutation.
//' 
//' @param rxy a function that generates x and y data.
//' @param TS routine to calculate test statistics for non-chi-square tests 
//' @param typeTS indicator for type of test statistics
//' @param TSextra additional info passed to TS, if necessary
//' @param alpha A numeric constant
//' @param samplingmethod =2, 1 for independence sampling, 2 for MCMC in discrete data case
//' @param B Number of simulation runs.
//' @param xparam  arguments for r1.
//' @param yparam  arguments for r2.
//' @keywords internal
//' @return A numeric matrix of powers
// [[Rcpp::export]]
NumericMatrix power_disc(Function rxy, 
                          Function TS,
                          int typeTS,
                          List TSextra,
                          double alpha=0.05, 
                          int samplingmethod=2,
                          NumericVector B=NumericVector::create(1000, 2000), 
                          NumericVector xparam=0.0,                       
                          NumericVector yparam=0.0) { 
   /* Find out how many tests are to be done */  
   List dta=rxy(xparam(0),yparam(0));
   IntegerVector x = as<IntegerVector>(dta["x"]);
   IntegerVector y = as<IntegerVector>(dta["y"]);
   NumericVector vals = as<NumericVector>(dta["vals"]);
   NumericVector data_TS, w;
   if(typeTS==5) {
     w=weights(dta);
     data_TS=TS(x, y, vals, w);
   }  
   if(typeTS==6) data_TS=TS(x, y);  
   if(typeTS==7) data_TS=TS(x, y, vals);
   if(typeTS==8) data_TS=TS(x, y, vals, TSextra);
   int const nummethods=data_TS.size();
   
   int rp=yparam.size(), i, j, l, k;
   NumericMatrix pwr(rp, nummethods),zeromatrix(rp, nummethods);
   NumericVector sim_data(nummethods),sim_perm(nummethods);
   IntegerVector counter(nummethods);
   
   /*  loop over values in xparam, yparam */    
   for(l=0;l<rp;++l) {
     /*  loop over simulation runs */  
     for(i=0;i<B[0];++i) {
       /*    create new data set  */            
       dta=rxy(xparam(l),yparam(l));
       /*    check that vectors are of equal length  
        and either x or y is greater than 0 for each value of vals. 
        If not write message and return matrix of zeros   */ 
       
       x = as<IntegerVector>(dta["x"]);
       y = as<IntegerVector>(dta["y"]);
       vals = as<NumericVector>(dta["vals"]);
       if( (x.size()!=vals.size()) || (y.size()!=vals.size()) ) {
         Rcout<<"Data generated has x, y and vals vectors of unequal lengths. Check your rxy function!\n";
         return zeromatrix;            
       }
       for(int i0=0;i0<vals.size();++i0) {
         if(x[i0]+y[i0]==0) {
           Rcout<<"Data generated has x=0 and y=0 for some value of vals. Check your rxy function!\n";
           return zeromatrix;                          
         }
       }
       /*       find test statistics for data  */              
       if(typeTS==5) {
         w=weights(dta);
         sim_data=TS(x, y, vals, w);
       }
       if(typeTS==6) sim_data=TS(x, y);
       if(typeTS==7) sim_data=TS(x, y, vals);
       if(typeTS==8) sim_data=TS(x, y, vals, TSextra);
       for(j=0;j<nummethods;++j) counter[j]=0;
       dta = permute_disc(dta, 1);
       for(k=0;k<B[1];++k) {     
         /*       find test statistics for permuted data   */ 
         GetRNGstate();
         dta = permute_disc(dta, samplingmethod);
         PutRNGstate();
         if(typeTS==5) sim_perm=TS(dta["x"], dta["y"], dta["vals"], w); 
         if(typeTS==6) sim_perm=TS(dta["x"], dta["y"]); 
         if(typeTS==7) sim_perm=TS(dta["x"], dta["y"], dta["vals"]); 
         if(typeTS==8) sim_perm=TS(dta["x"], dta["y"], dta["vals"], TSextra); 
         /* and compare with test statistics from data */       
         for(j=0;j<nummethods;++j) {
           if(sim_perm(j)>sim_data(j)) counter[j]=counter[j]+1;
         }
       } /* end of loop of p value simulation*/ 
         
         /*    p values < alpha?  */     
         
         for(j=0;j<nummethods;++j) 
           if(double(counter[j])/B[1]<alpha) 
             pwr(l, j)=pwr(l, j)+1.0;
     }  /* end of loop of power simulation*/ 
         
         for(j=0;j<nummethods;++j) pwr(l, j)=pwr(l, j)/B[0];
   }  /* end of loop over cases */
         return  pwr;
   
 }
