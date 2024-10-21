#include <Rcpp.h>
using namespace Rcpp;

//' find test statistics for discrete data
//' 
//' @param x  vector of weights data set 1
//' @param y  vector of weights data set 2
//' @keywords internal
//' @return A vector of test statistics
// [[Rcpp::export]]

NumericVector TSw_disc(std::vector<double>& x, 
                       std::vector<double>& y)  {
  
  CharacterVector methods = CharacterVector::create("KS", "Kuiper", "CvM", "AD");
  int const nummethods=4;

  int k=x.size(), i;
  double nx, ny, n;
  NumericVector TS(nummethods);
  TS.names() = methods;
  std::vector<double> xy(k), edfx(k), edfy(k), edfxy(k);

  double tmp;

  /*  sample sizes*/  
  
  nx=0.0;
  ny=0.0;
  for(i=0;i<k;++i) {
    nx+=x[i];
    ny+=y[i];
  } 
  if(nx<1) nx=1;
  if(ny<1) ny=1;
  n=nx+ny;

  /*  combined data, then sorted, cumulative sums of x, y and xy */
  
  edfx[0]=x[0]/nx;
  edfy[0]=y[0]/ny;
  xy[0]=x[0]+y[0];
  edfxy[0]=xy[0]/n;
  for(i=1;i<k;++i) {
    edfx[i]=edfx[i-1]+x[i]/nx;
    edfy[i]=edfy[i-1]+y[i]/ny;
    xy[i]=x[i]+y[i];
    edfxy[i]=edfxy[i-1]+xy[i]/n;    
  }  

  /*  Kolomogorov-Smirnov and Kuiper */
  
     double mx=edfx[0]-edfy[0];
     double Mx=mx;
     for(i=1;i<k;++i) {
       tmp=edfx[i]-edfy[i];
       if(mx>tmp) mx=tmp;
       if(Mx<tmp) Mx=tmp;
     }
     if(-mx>Mx) TS[0]=-mx;
     else TS[0]=Mx;
     TS[1]=Mx-mx;
     
  /* Cramer-vonMises and Anderson-Darling test*/
   
     tmp=edfx[0]-edfy[0];     
     TS[2]=tmp*tmp*xy[0];
     TS[3]=tmp*tmp/(1-edfxy[0]); 
     for(i=1;i<k;++i) {
          tmp=edfx[i]-edfy[i];     
          TS[2]=TS[2]+xy[i]*tmp*tmp;
          if(edfxy[i]<1) 
            TS[3]=TS[3]+tmp*tmp/edfxy[i]/(1-edfxy[i])*(edfxy[i]-edfxy[i-1]);
     }
     TS[2]=TS[2]*nx*ny/n/n;    
     TS[3]=TS[3]*nx*ny; 
 
  return TS;
}
