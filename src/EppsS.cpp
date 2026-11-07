#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include "EppsS.h"
using namespace Rcpp;

/* Epps-Singleton two-sample test (Epps and Singleton 1986), as in
   scipy.stats.epps_singleton_2samp: the empirical characteristic functions of
   the two samples are compared at the points t/sigma, t=(0.4, 0.8), where sigma
   is the semi-interquartile range of the pooled sample.

   The data are given as sorted values v with counts cx and cy of the two
   samples, so the same code works for continuous data (each value has count
   1 in one of the samples) and for discrete data. */

/* type 7 quantile (the R default) of the pooled sample */
static double pooled_quantile(const std::vector<double>& v,
                              const std::vector<double>& c, double N, double p) {
  double h=(N-1.0)*p;
  double lo=std::floor(h);
  /* value of the order statistic with 0-based index j */
  auto value_at=[&](double j) {
    double cum=0.0;
    for(std::size_t i=0;i<v.size();++i) {
      cum+=c[i];
      if(j<cum) return v[i];
    }
    return v.back();
  };
  double a=value_at(lo), b=value_at(std::min(lo+1.0, N-1.0));
  return a+(h-lo)*(b-a);
}

/* eigen decomposition of a symmetric k x k matrix by Jacobi rotations */
static void jacobi_eigen(std::vector<double> A, int k,
                         std::vector<double>& evals, std::vector<double>& evecs) {
  evecs.assign(k*k, 0.0);
  for(int i=0;i<k;++i) evecs[i*k+i]=1.0;
  for(int sweep=0;sweep<100;++sweep) {
    double off=0.0;
    for(int p=0;p<k;++p) for(int q=p+1;q<k;++q) off+=A[p*k+q]*A[p*k+q];
    if(off<1e-30) break;
    for(int p=0;p<k;++p) {
      for(int q=p+1;q<k;++q) {
        double apq=A[p*k+q];
        if(std::abs(apq)<1e-300) continue;
        double theta=(A[q*k+q]-A[p*k+p])/(2.0*apq);
        double t=(theta>=0 ? 1.0 : -1.0)/(std::abs(theta)+std::sqrt(theta*theta+1.0));
        double cs=1.0/std::sqrt(t*t+1.0), sn=t*cs;
        for(int r=0;r<k;++r) {         /* A = A J */
          double arp=A[r*k+p], arq=A[r*k+q];
          A[r*k+p]=cs*arp-sn*arq;
          A[r*k+q]=sn*arp+cs*arq;
        }
        for(int r=0;r<k;++r) {         /* A = J^T A */
          double apr=A[p*k+r], aqr=A[q*k+r];
          A[p*k+r]=cs*apr-sn*aqr;
          A[q*k+r]=sn*apr+cs*aqr;
        }
        for(int r=0;r<k;++r) {         /* V = V J */
          double vrp=evecs[r*k+p], vrq=evecs[r*k+q];
          evecs[r*k+p]=cs*vrp-sn*vrq;
          evecs[r*k+q]=sn*vrp+cs*vrq;
        }
      }
    }
  }
  evals.assign(k, 0.0);
  for(int i=0;i<k;++i) evals[i]=A[i*k+i];
}

double eppss_core(const std::vector<double>& v,
                  const std::vector<double>& cx,
                  const std::vector<double>& cy,
                  int* rank) {
  const int k=4;
  std::size_t K=v.size();
  double nx=0.0, ny=0.0;
  std::vector<double> c(K);
  for(std::size_t i=0;i<K;++i) {
    nx+=cx[i]; ny+=cy[i]; c[i]=cx[i]+cy[i];
  }
  double N=nx+ny;
  if(rank) *rank=0;
  if(nx<1 || ny<1) return 0.0;
  double sigma=(pooled_quantile(v, c, N, 0.75)-pooled_quantile(v, c, N, 0.25))/2.0;
  if(!(sigma>0.0)) {
    /* interquartile range 0 (e.g. discrete data with few values): use the
       standard deviation of the pooled sample instead */
    double m=0.0, s=0.0;
    for(std::size_t i=0;i<K;++i) m+=c[i]*v[i];
    m/=N;
    for(std::size_t i=0;i<K;++i) s+=c[i]*(v[i]-m)*(v[i]-m);
    sigma=std::sqrt(s/N);
    if(!(sigma>0.0)) return 0.0;
  }
  double ts[2]={0.4/sigma, 0.8/sigma};
  /* g(v) = (cos(t1 v), cos(t2 v), sin(t1 v), sin(t2 v)) */
  std::vector<double> mx(k,0.0), my(k,0.0), Sx(k*k,0.0), Sy(k*k,0.0), g(k);
  for(std::size_t i=0;i<K;++i) {
    /* ts[1]=2*ts[0], so the double angle formulas save two trigonometric calls */
    double co=std::cos(ts[0]*v[i]), si=std::sin(ts[0]*v[i]);
    g[0]=co; g[1]=2.0*co*co-1.0;
    g[2]=si; g[3]=2.0*si*co;
    for(int a=0;a<k;++a) {
      mx[a]+=cx[i]*g[a];
      my[a]+=cy[i]*g[a];
      for(int b=0;b<k;++b) {
        Sx[a*k+b]+=cx[i]*g[a]*g[b];
        Sy[a*k+b]+=cy[i]*g[a]*g[b];
      }
    }
  }
  for(int a=0;a<k;++a) {mx[a]/=nx; my[a]/=ny;}
  /* S = (N/nx) Cov_x + (N/ny) Cov_y, covariances with divisor n */
  std::vector<double> S(k*k);
  for(int a=0;a<k;++a)
    for(int b=0;b<k;++b)
      S[a*k+b]=N/nx*(Sx[a*k+b]/nx-mx[a]*mx[b]) + N/ny*(Sy[a*k+b]/ny-my[a]*my[b]);
  std::vector<double> ev, V;
  jacobi_eigen(S, k, ev, V);
  double emax=0.0;
  for(int a=0;a<k;++a) if(ev[a]>emax) emax=ev[a];
  /* w = N d' S^+ d with the pseudo-inverse of S */
  double w=0.0;
  int r=0;
  for(int a=0;a<k;++a) {
    if(!(ev[a]>emax*1e-12) || emax<=0.0) continue;
    ++r;
    double proj=0.0;
    for(int b=0;b<k;++b) proj+=V[b*k+a]*(mx[b]-my[b]);
    w+=proj*proj/ev[a];
  }
  w*=N;
  /* small sample correction */
  if(std::max(nx, ny)<25.0)
    w/=(1.0+std::pow(N, -0.45)+10.1*(std::pow(nx, -1.7)+std::pow(ny, -1.7)));
  if(rank) *rank=r;
  return w;
}

//' Epps-Singleton test statistic and its degrees of freedom
//'
//' @param x first continuous data set, or the counts of the first data set
//' @param y second continuous data set, or the counts of the second data set
//' @param vals values of discrete data, NULL for continuous data
//' @keywords internal
//' @return vector with the statistic and the degrees of freedom of its
//'   asymptotic chi-square null distribution
//' @noRd
// [[Rcpp::export]]
NumericVector EppsS(NumericVector x, NumericVector y,
                    Rcpp::Nullable<NumericVector> vals=R_NilValue) {
  std::vector<double> v, cx, cy;
  if(vals.isNull()) {
    int nx=x.size(), ny=y.size();
    std::vector< std::pair<double,int> > z(nx+ny);
    for(int i=0;i<nx;++i) z[i]=std::make_pair(x[i], 1);
    for(int i=0;i<ny;++i) z[nx+i]=std::make_pair(y[i], 0);
    std::sort(z.begin(), z.end());
    for(std::size_t i=0;i<z.size();++i) {
      v.push_back(z[i].first);
      cx.push_back(z[i].second);
      cy.push_back(1-z[i].second);
    }
  } else {
    NumericVector vv(vals);
    for(int i=0;i<vv.size();++i) {
      v.push_back(vv[i]); cx.push_back(x[i]); cy.push_back(y[i]);
    }
  }
  int r;
  double w=eppss_core(v, cx, cy, &r);
  return NumericVector::create(Named("statistic")=w, Named("df")=r);
}
