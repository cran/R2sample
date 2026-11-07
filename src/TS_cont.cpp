#include <Rcpp.h>
#include "EppsS.h"
using namespace Rcpp;

//' find test statistics for continuous data
//'
//' @param x_input first continuous data set
//' @param y_input second continuous data set
//' @keywords internal
//' @return A vector of test statistics
//' @noRd
// [[Rcpp::export]]
NumericVector TS_cont(const NumericVector& x_input,
                      const NumericVector& y_input) {

  // Rcpp vectors can share storage with the R objects supplied by the caller.
  // This routine sorts x and y, so work on private copies to avoid modifying
  // the user's input vectors.
  NumericVector x = clone(x_input);
  NumericVector y = clone(y_input);

  CharacterVector methods = CharacterVector::create("KS", "Kuiper", "CvM", "AD", "ZA", "ZK", "ZC", "Wassp1", "EppsS");
  int const nummethods=methods.size();

  int nx=x.size(), ny=y.size(), n=nx+ny, i;
  NumericVector TS(nummethods), Fx(n), Fy(n), px(n), py(n), sxy(n);
  double tmp;

  TS.names() =  methods;

  /*  sort data, combined sorted sample */

  std::sort(x.begin(), x.end());
  std::sort(y.begin(), y.end());
  for(i=0;i<nx;++i) sxy[i]=x[i];
  for(i=0;i<ny;++i) sxy[i+nx]=y[i];
  std::sort(sxy.begin(), sxy.end());

  /* empirical distribution functions of x and y evaluated on the
     combined sorted data set: Fx[i] = #{x <= sxy[i]}/nx.
     px, py are the "midpoint" edfs used by Zhang's tests:
     px[i] = ( #{x < sxy[i]} + #{x == sxy[i]}/2 )/nx.
     Without ties this is Fx[i]-0.5/nx if sxy[i] is an x and Fx[i] otherwise,
     as in Zhang (2006). */

  for(i=0;i<n;++i) {
    double lx = std::lower_bound(x.begin(), x.end(), sxy[i]) - x.begin();
    double ux = std::upper_bound(x.begin(), x.end(), sxy[i]) - x.begin();
    double ly = std::lower_bound(y.begin(), y.end(), sxy[i]) - y.begin();
    double uy = std::upper_bound(y.begin(), y.end(), sxy[i]) - y.begin();
    Fx[i] = ux/nx;
    Fy[i] = uy/ny;
    px[i] = (lx+ux)/2.0/nx;
    py[i] = (ly+uy)/2.0/ny;
  }

  /*  Kolmogorov-Smirnov and Kuiper Tests*/

    double mx=Fx[0]-Fy[0], Mx=Fx[0]-Fy[0];
    for(i=1;i<n;++i) {
      tmp = Fx[i]-Fy[i];
      if(mx>tmp) mx=tmp;
      if(Mx<tmp) Mx=tmp;
    }
    if(-mx>Mx) TS(0)=-mx;
    else TS(0)=Mx;
    TS(1)=Mx-mx;

  /* Cramer-vonMises test (Anderson 1962) */

    TS(2)=sum((Fx-Fy)*(Fx-Fy))*nx*ny/n/n;

  /*  Anderson-Darling test (Pettitt 1976):
      nm sum_{j=1}^{N-1} (F_j-G_j)^2/(j(N-j)) */

    TS(3)=0.0;
    for(i=0;i<n-1;++i) {
      tmp=Fx[i]-Fy[i];
      TS(3)=TS(3)+tmp*tmp/(double(i)+1.0)/(double(n)-double(i)-1.0);
    }
    TS(3)=TS(3)*nx*ny;

  /*   Zhang's tests (Zhang 2006). ZA and ZC are returned with the sign
       changed, so that for all tests large values lead to rejection. */

    TS(4)=0.0;
    TS(5)=0.0;
    for(i=0;i<n;++i) {
      double ex=px[i]*log(px[i]+1e-10)+(1-px[i])*log(1-px[i]+1e-10);
      double ey=py[i]*log(py[i]+1e-10)+(1-py[i])*log(1-py[i]+1e-10);
      double p=(i+0.5)/n;
      TS(4)=TS(4)+(nx*ex+ny*ey)/((i+0.5)*(n-i-0.5));
      tmp=nx*ex+ny*ey-n*(p*log(p)+(1-p)*log(1-p));
      if(tmp>TS(5))  TS(5)=tmp;
    }

  /*  ZC uses the (mid)ranks of x and y in the combined sample */

    tmp=0.0;
    for(i=0;i<nx;++i) {
      double l = std::lower_bound(sxy.begin(), sxy.end(), x[i]) - sxy.begin();
      double u = std::upper_bound(sxy.begin(), sxy.end(), x[i]) - sxy.begin();
      double R = (l+u+1.0)/2.0;
      tmp=tmp+log(nx/(i+0.5)-1)*log(n/(R-0.5)-1);
    }
    for(i=0;i<ny;++i) {
      double l = std::lower_bound(sxy.begin(), sxy.end(), y[i]) - sxy.begin();
      double u = std::upper_bound(sxy.begin(), sxy.end(), y[i]) - sxy.begin();
      double R = (l+u+1.0)/2.0;
      tmp=tmp+log(ny/(i+0.5)-1)*log(n/(R-0.5)-1);
    }
    TS(6)=-tmp/n;

 /* Wasserstein p=1:  W1 = integral |Fx(t)-Fy(t)| dt
    = sum_i |Fx(z_i)-Fy(z_i)| (z_{i+1}-z_i) over the combined sorted sample */

    TS(7)=0.0;
    for(i=0;i<n-1;++i)
      TS(7)=TS(7)+std::abs(Fx[i]-Fy[i])*(sxy[i+1]-sxy[i]);

 /* Epps-Singleton test (Epps and Singleton 1986): merge the sorted samples
    into sorted values with counts for x and y */

    {
      std::vector<double> v(n), cx(n), cy(n);
      int a=0, b=0;
      for(i=0;i<n;++i) {
        if(b>=ny || (a<nx && x[a]<=y[b])) {v[i]=x[a]; cx[i]=1.0; cy[i]=0.0; ++a;}
        else {v[i]=y[b]; cx[i]=0.0; cy[i]=1.0; ++b;}
      }
      TS(8)=eppss_core(v, cx, cy, NULL);
    }

  return TS;
}
