#include <Rcpp.h>
#include "EppsS.h"
using namespace Rcpp;

//' find test statistics for discrete data
//'
//' @param x integer vector of data set 1
//' @param y integer vector of data set 2
//' @param vals numeric vector of values of discrete data set, sorted increasingly
//' @param ADweights A vector of weights for AD method, see weights()
//' @keywords internal
//' @return A vector of test statistics
//' @noRd
// [[Rcpp::export]]
NumericVector TS_disc(IntegerVector x,
                      IntegerVector y,
                      NumericVector vals,
                      NumericVector ADweights=Rcpp::NumericVector::create(2))  {

  CharacterVector methods = CharacterVector::create("KS", "Kuiper", "CvM", "AD", "ZA", "Wassp1", "EppsS");
  int const nummethods=methods.size();

  int k=vals.size(), nx, ny, n, i;
  NumericVector TS(nummethods), Fx(k), Fy(k);
  IntegerVector xy(k), cx(k), cy(k), cxy(k);
  double tmp;

  TS.names() = methods;
  if((int) ADweights.size()!=k) stop("ADweights must have the same length as vals");

  /*  sample sizes*/
  nx=0;
  ny=0;
  for(i=0;i<k;++i) {
    nx+=x[i];
    ny+=y[i];
  }
  n=nx+ny;

  /*  cumulative sums of x, y and xy, edfs */

  for(i=0;i<k;++i) {
    xy[i]=x[i]+y[i];
    cx[i]=x[i];
    cy[i]=y[i];
    cxy[i]=xy[i];
    if(i>0) {
      cx[i]+=cx[i-1];
      cy[i]+=cy[i-1];
      cxy[i]+=cxy[i-1];
    }
    Fx[i]=double(cx[i])/nx;
    Fy[i]=double(cy[i])/ny;
  }

  /*  Kolmogorov-Smirnov and Kuiper */

     double mx=Fx[0]-Fy[0];
     double Mx=Fx[0]-Fy[0];
     for(i=1;i<k;++i) {
       tmp=Fx[i]-Fy[i];
       if(mx>tmp) mx=tmp;
       if(Mx<tmp) Mx=tmp;
     }
     if(-mx>Mx) TS(0)=-mx;
     else TS(0)=Mx;
     TS(1)=Mx-mx;

  /* Cramer-vonMises and Anderson-Darling test*/

     TS(2)=0.0;
     TS(3)=0.0;
     for(i=0;i<k;++i) {
          tmp=Fx(i)-Fy(i);
          TS(2)=TS(2)+xy(i)*tmp*tmp;
          TS(3)=TS(3)+ADweights(i)*tmp*tmp;
     }
     TS(2)=TS(2)*nx*ny/n/n;
     TS(3)=TS(3)*nx*ny;

  /* Zhang's ZA test. Within a class the "midpoint" edfs
     px = (cx[i-1]+x[i]/2)/nx and py = (cy[i-1]+y[i]/2)/ny are used at each
     of the tied positions j=cxy[i-1]+1,..,cxy[i] of the combined sample.
     Without ties this is exactly Zhang's (2006) statistic, as in TS_cont.
     The sign is changed so that large values lead to rejection. */

     TS(4)=0.0;
     for(i=0;i<k;++i) {
       double px=(double(cx[i])-x[i]/2.0)/nx;
       double py=(double(cy[i])-y[i]/2.0)/ny;
       double ex=px*log(px+1e-10)+(1-px)*log(1-px+1e-10);
       double ey=py*log(py+1e-10)+(1-py)*log(1-py+1e-10);
       /* w = sum_{j=a}^{b} 1/((j-0.5)(n-j+0.5))
            = ( sum 1/(j-0.5) + sum 1/(n-j+0.5) )/n, via digamma */
       int a=cxy[i]-xy[i]+1, b=cxy[i];
       double w=( R::digamma(b+0.5)-R::digamma(a-0.5) +
                  R::digamma(n-a+1.5)-R::digamma(n-b+0.5) )/n;
       TS(4)=TS(4)+(nx*ex+ny*ey)*w;
     }

 /*  Wasserstein */

     TS(5)=0.0;
     for(i=0;i<k-1;++i) {
       TS(5)=TS(5)+std::abs(Fx(i)-Fy(i))*(vals(i+1)-vals(i));
     }

 /* Epps-Singleton test (Epps and Singleton 1986) */

     {
       std::vector<double> v(k), wx(k), wy(k);
       for(i=0;i<k;++i) {v[i]=vals[i]; wx[i]=x[i]; wy[i]=y[i];}
       TS(6)=eppss_core(v, wx, wy, NULL);
     }

  return TS;
}
