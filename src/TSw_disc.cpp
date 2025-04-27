#include <Rcpp.h>
using namespace Rcpp;

//' Find test statistics for weighted discrete data
//' 
//' @param x integer vector of counts
//' @param y integer vector of counts
//' @param vals A numeric vector with the values of the discrete rv.
//' @param wx integer vector of weights
//' @param wy integer vector of weights
//' @keywords internal
//' @return A vector with test statistics
// [[Rcpp::export]]
Rcpp::NumericVector TSw_disc(
        Rcpp::IntegerVector x, 
        Rcpp::IntegerVector y,
        Rcpp::NumericVector vals,  
        Rcpp::NumericVector wx,
        Rcpp::NumericVector wy) {
    
  Rcpp::CharacterVector methods=CharacterVector::create("KS", "K", "CvM", "AD");    
  int const nummethods=methods.size();
  int k=x.size(), i;
  NumericVector TS(nummethods), edfx(k), edfy(k), edfz(k);
  double nx, ny, n, tmp;
  TS.names() =  methods;

  /*  Find sample sizes, cumulative sums*/
  
  nx=0.0;
  ny=0.0;
  for(i=0;i<k;++i) {
    nx = nx + x[i]*wx[i];
    ny = ny + y[i]*wy[i];
  }
  n=nx+ny;
  
  /* Edfs */
  
  edfx[0] = x[0]*wx[0]/nx;
  edfy[0] = y[0]*wy[0]/ny;
  edfz[0] = (x[0]*wx[0]+y[0]*wy[0])/n;
  for(i=1;i<k;++i) {
    edfx[i] = edfx[i-1] + x[i]*wx[i]/nx;
    edfy[i] = edfy[i-1] + y[i]*wy[i]/ny;
    edfz[i] = edfz[i-1] + (x[i]*wx[i]+y[i]*wy[i])/n;
  }
  /*  Kolmogorov-Smirnov and Kuiper*/
  double mx = 0;
  double Mx = 0;
  for(i=0;i<k-1;++i) {
      tmp = edfx[i] - edfy[i];
      if(tmp<0 && std::abs(tmp)>std::abs(mx)) mx=std::abs(tmp);
      if(tmp>0 && std::abs(tmp)>std::abs(Mx)) Mx=std::abs(tmp);      
  }
  if(std::abs(mx)>std::abs(Mx)) TS(0)=std::abs(mx);
  else TS(0)=std::abs(Mx);
  TS(1)=Mx+mx; 
 
 /* Cramer-von Mises and  Anderson-Darling */
 
  tmp = (edfx[0]-edfy[0])*(edfx[0]-edfy[0]);
  TS(2)=tmp*edfz[0];
  if(edfz[0]<1) TS(3) = tmp/(1-edfz[0]);
  else TS(3) = 0.0;    
  for(i=1;i<k;++i) {
      tmp = (edfx[i]-edfy[i])*(edfx[i]-edfy[i]);
      TS(2) = TS(2) + tmp*(edfz[i]-edfz[i-1]);
      if( (edfz[i]>0)&&(edfz[i]<1) ) {
        TS(3) = TS(3) + tmp/edfz[i]/(1-edfz[i])*(edfz[i]-edfz[i-1]); 
      }     
  }
  
  return TS;
}
