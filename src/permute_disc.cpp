#include <Rcpp.h>
#include <cmath>
#include <random>

using namespace Rcpp; 

//' permute discrete data
//' 
//' @param dta A list of numeric vectors.
//' @param samplingmethod  =2 , 1 for independence sampling by expansion or 2 for MCMC
//' @keywords internal
//' @return A list of permuted vectors
// [[Rcpp::export]]
Rcpp::List permute_disc(const Rcpp::List dta, int samplingmethod=1) {

  NumericVector dtax = as<NumericVector>(dta["x"]); 
  NumericVector dtay = as<NumericVector>(dta["y"]);

  int nx, ny, n, i, j,  d, k=dtax.size();
  NumericVector xy(k), xout(k), yout(k);
  IntegerVector Index(k);
  std::vector<int> x(k), y(k);
  std::vector<double> rx(k), ry(k);
  nx=0;
  ny=0;
  for(i=0;i<k;++i) {
    Index[i]=i;
    x[i] = floor(dtax[i]);
    y[i] = floor(dtay[i]);
    xy[i] = x[i]+y[i];
    nx=nx+x[i];
    ny=ny+y[i];
    rx[i] = dtax[i]-x[i];
    ry[i] = dtay[i]-y[i];
  }
  n=nx+ny;
  if(samplingmethod==1) {
    IntegerVector z(n), zx(nx), zy(ny);
    int l=0;
    for(i=0;i<k;++i) {
      for(j=0;j<x[i];++j) {z[l]=i;++l;}
    }  
    l=0;
    for(i=0;i<k;++i) {
      for(j=0;j<y[i];++j) {z[nx+l]=i;++l;}
    }
    z = Rcpp::sample(z, n);
    for(i=0;i<nx;++i) zx[i]=z[i];
    for(i=0;i<ny;++i) zy[i]=z[nx+i];
    std::sort(zx.begin(), zx.end());
    std::sort(zy.begin(), zy.end());
    j=0;
    for(i=0;i<k;++i) {
      x[i]=0;
      while (j<nx && zx[j]==i) {
        x[i]=x[i]+1;
        ++j;
      }
    } 
    j=0;
    for(i=0;i<k;++i) {
      y[i]=0;
      while (j<ny && zy[j]==i) {
        y[i]=y[i]+1;
        ++j;
      } 
    } 
  }
  if(samplingmethod==2) {
      IntegerVector ij=Rcpp::sample(Index, 2, false);
      i=ij[0];
      j=ij[1];    
      d=1;
      double tmp=(xy[i]-x[i])/(xy[j]+1.0-x[j])*double(x[j])/(x[i]+1.0);
      if(R::runif(0, 1)<0.5 && x[j]>5 && y[i]>5) {
        d=2;
        tmp=tmp*(xy[i]-x[i]-1.0)/(xy[j]+2.0-x[j])*double(x[j]+1)/(x[i]+2.0);
      }
      if(R::runif(0, 1)<0.25 && x[j]>10 && y[i]>10) {
        d=3;
        tmp=tmp*(xy[i]-x[i]-2.0)/(xy[j]+3.0-x[j])*double(x[j]+2.0)/(x[i]+3.0);
      }      
      if(R::runif(0, 1)>tmp) return dta;
      else {
        x[i]=x[i]+d;
        x[j]=x[j]-d;
        y[i]=y[i]-d;
        y[j]=y[j]+d;
      }
  }
  for(i=0;i<k;++i) {
    xout[i]=x[i];
    yout[i]=y[i];
    if((R::runif(0, 1)<0.5)) xout[i]=xout[i]+rx[i];
    else yout[i]=yout[i]+rx[i];
    if((R::runif(0, 1)<0.5)) xout[i]=xout[i]+ry[i];
    else yout[i]=yout[i]+ry[i];
  }
  
  return List::create(Named("x") = xout, 
                      Named("y") = yout,
                      Named("vals") = dta["vals"]);                    
}
