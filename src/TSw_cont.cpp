#include <Rcpp.h>
using namespace Rcpp;

//' find test statistics for continuous data with weights
//'
//' Weighted edfs Fx, Fy and weighted pooled edf H are used:
//' KS and Kuiper are based on D=Fx-Fy,
//' CvM = nm/N sum D^2 dH and AD = nm/N sum D^2 dH/(H(1-H)).
//' With all weights equal these are identical to the statistics in TS_cont.
//'
//' @param x first continuous data set
//' @param y second continuous data set
//' @param wx weights of x
//' @param wy weights of y
//' @keywords internal
//' @return A vector of test statistics
//' @noRd
// [[Rcpp::export]]

NumericVector TSw_cont(std::vector<double> x,
                       std::vector<double> y,
                       std::vector<double> wx,
                       std::vector<double> wy) {

  CharacterVector methods = CharacterVector::create("KS", "Kuiper", "CvM", "AD");
  int const nummethods=methods.size();

  int nx=x.size(), ny=y.size(), n=nx+ny, i, j;
  NumericVector TS(nummethods);
  TS.names() =  methods;

  /* combine data, remember sample membership, sort by value */

  std::vector<double> z(n), w(n);
  std::vector<int> g(n), idx(n);
  double Wx=0.0, Wy=0.0;
  for(i=0;i<nx;++i) {z[i]=x[i]; w[i]=wx[i]; g[i]=1; Wx+=wx[i];}
  for(i=0;i<ny;++i) {z[nx+i]=y[i]; w[nx+i]=wy[i]; g[nx+i]=2; Wy+=wy[i];}
  double W=Wx+Wy;
  for(i=0;i<n;++i) idx[i]=i;
  std::stable_sort(idx.begin(), idx.end(),
                   [&](int A, int B) -> bool {return z[A] < z[B];});

  /* run through groups of tied values. D is evaluated after the whole
     group has been added, as in TS_cont. */

  double Fx=0.0, Fy=0.0, H=0.0, m=0.0, M=0.0, cvm=0.0, ad=0.0;
  bool first=true;
  i=0;
  while(i<n) {
    j=i;
    while(j<n && z[idx[j]]==z[idx[i]]) {
      if(g[idx[j]]==1) Fx+=w[idx[j]]/Wx;
      else Fy+=w[idx[j]]/Wy;
      ++j;
    }
    double D=Fx-Fy;
    if(first) {m=D; M=D; first=false;}
    if(D<m) m=D;
    if(D>M) M=D;
    for(int l=i;l<j;++l) {
      double dH=w[idx[l]]/W;
      H+=dH;
      cvm+=D*D*dH;
      if(l<n-1 && H>0.0 && H<1.0) ad+=D*D*dH/(H*(1.0-H));
    }
    i=j;
  }
  if(-m>M) TS(0)=-m;
  else TS(0)=M;
  TS(1) = M-m;
  TS(2) = cvm*nx*ny/n;
  TS(3) = ad*nx*ny/n;

  return TS;

}
