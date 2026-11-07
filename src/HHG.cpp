#include <Rcpp.h>
#include <algorithm>
#include <vector>
using namespace Rcpp;

/* HHG-type K-sample test (Heller, Heller, Kaufman, Brill and Gorfine,
   JMLR 17, 2016), equipartition version, for two samples.

   The pooled sample is divided into A "atoms": for continuous data A groups of
   consecutive order statistics of (nearly) equal size, for discrete data the A
   distinct values. For each number of cells m=2,..,mmax the statistic S_m is the
   average, over all partitions of the atoms into m contiguous cells, of the
   likelihood ratio statistic of the resulting 2 x m table of counts.
   The S_m are combined with the MinP rule: each S_m is turned into a p value
   using a stored sample of its null distribution, and the test statistic is
   -log(min_m p_m), so large values lead to rejection. */

/* lchoose table: lc[n][k] = log(choose(n,k)), n,k=0..A. It is kept between
   calls because the same A is used for every data set of a test. */
static const std::vector< std::vector<double> >& lchoose_table(int A) {
  static int lastA=-1;
  static std::vector< std::vector<double> > lc;
  if(A!=lastA) {
    lc.assign(A+1, std::vector<double>(A+1, R_NegInf));
    for(int n=0;n<=A;++n)
      for(int k=0;k<=n;++k) lc[n][k]=R::lchoose(n, k);
    lastA=A;
  }
  return lc;
}

static double xlogx_term(double c, double num) {
  /* c*log(num), with 0*log(0)=0 */
  if(c<=0.0) return 0.0;
  return c*std::log(num);
}

/* S_m, m=2,..,mmax, from cumulative counts.
   cx[j], ct[j]: number of x's and of all observations in atoms 1..j, j=0..A */
static NumericVector hhg_stats(const std::vector<double>& cx,
                               const std::vector<double>& ct,
                               int A, int mmax) {
  double nx=cx[A], N=ct[A], ny=N-nx;
  if(mmax<2 || A<2) return NumericVector(0);
  NumericVector S(mmax-1);
  const std::vector< std::vector<double> >& lc = lchoose_table(A);
  for(int a=0;a<A;++a) {
    for(int b=a+1;b<=A;++b) {
      double X=cx[b]-cx[a], T=ct[b]-ct[a], Y=T-X;
      if(T<=0.0) continue;
      /* likelihood ratio contribution of the cell made of atoms a+1,..,b */
      double score = 2.0*( xlogx_term(X, X*N/(nx*T)) + xlogx_term(Y, Y*N/(ny*T)) );
      int r = (a>0) + (b<A);                       /* cuts fixed by the cell */
      int f = (a>0 ? a-1 : 0) + (b<A ? A-1-b : 0); /* free cut positions */
      for(int m=2;m<=mmax;++m) {
        int k=m-1-r;
        if(k<0 || k>f) continue;
        /* fraction of all partitions into m cells that contain this cell */
        S[m-2] += score*std::exp(lc[f][k]-lc[A-1][m-1]);
      }
    }
  }
  return S;
}

/* boundaries of the atoms for continuous data: b_j = round(j*N/A) */
static std::vector<int> atom_bounds(int N, int A) {
  std::vector<int> b(A+1);
  for(int j=0;j<=A;++j) b[j]=(int) std::floor((double)j*N/A+0.5);
  b[0]=0; b[A]=N;
  return b;
}

/* -log(min_m p_m), p_m=(1+#{null_m>=S_m})/(B+1), null columns sorted increasingly.
   If S is itself one of the rows of the null sample (self=true), that row is left
   out of its own comparison: p_m=(1+#{other null_m>=S_m})/B. */
static double minp_stat(const NumericVector& S, const NumericMatrix& nullsorted,
                        bool self=false) {
  int B=nullsorted.nrow(), M=S.size();
  double minp=1.0;
  for(int m=0;m<M && m<nullsorted.ncol();++m) {
    NumericMatrix::ConstColumn col=nullsorted(_, m);
    double tol=1e-9*std::abs(S[m]);
    int nless = std::lower_bound(col.begin(), col.end(), S[m]-tol) - col.begin();
    double p = self ? (B-nless)/double(B) : (1.0+B-nless)/(B+1.0);
    if(p<minp) minp=p;
  }
  return -std::log(minp);
}

//' Sizes of the atoms of the HHG test for continuous data
//' @param N total sample size
//' @param A number of atoms
//' @keywords internal
//' @return integer vector of length A
//' @noRd
// [[Rcpp::export]]
IntegerVector HHG_atoms(int N, int A) {
  std::vector<int> b=atom_bounds(N, A);
  IntegerVector t(A);
  for(int j=0;j<A;++j) t[j]=b[j+1]-b[j];
  return t;
}

//' Null distribution of the HHG statistics S_m
//'
//' The x-counts of the atoms are drawn from the multivariate hypergeometric
//' distribution, which is their distribution under random relabeling of the
//' pooled sample.
//' @param t integer vector, number of observations in each atom
//' @param nx size of the first sample
//' @param mmax largest number of cells
//' @param B number of simulation runs
//' @keywords internal
//' @return a B x (mmax-1) matrix, each column sorted increasingly, with attribute
//'   "minp", the sorted null values of the combined statistic -log(min p)
//' @noRd
// [[Rcpp::export]]
NumericMatrix HHG_null(IntegerVector t, int nx, int mmax, int B) {
  int A=t.size(), N=0;
  for(int j=0;j<A;++j) N+=t[j];
  mmax=std::min(mmax, A);
  if(mmax<2) stop("HHG needs at least two atoms");
  NumericMatrix out(B, mmax-1);
  std::vector<double> cx(A+1), ct(A+1);
  ct[0]=0.0;
  for(int j=0;j<A;++j) ct[j+1]=ct[j]+t[j];
  for(int i=0;i<B;++i) {
    int left=nx, Nleft=N;
    cx[0]=0.0;
    for(int j=0;j<A;++j) {
      int xj = (j==A-1) ? left : (int) R::rhyper(t[j], Nleft-t[j], left);
      cx[j+1]=cx[j]+xj;
      left-=xj;
      Nleft-=t[j];
    }
    NumericVector S=hhg_stats(cx, ct, A, mmax);
    for(int m=0;m<mmax-1;++m) out(i, m)=S[m];
  }
  NumericMatrix raw=clone(out);
  for(int m=0;m<mmax-1;++m) {
    NumericMatrix::Column col=out(_, m);
    std::sort(col.begin(), col.end());
  }
  /* null distribution of the combined (MinP) statistic itself, used for p values
     without permutation (large sample sizes, continuous data) */
  NumericVector minp(B);
  for(int i=0;i<B;++i) minp[i]=minp_stat(NumericVector(raw(i, _)), out, true);
  std::sort(minp.begin(), minp.end());
  out.attr("minp")=minp;
  return out;
}

//' HHG statistics S_m for continuous data
//' @param x first data set
//' @param y second data set
//' @param A number of atoms
//' @param mmax largest number of cells
//' @keywords internal
//' @return numeric vector of S_m, m=2,..,mmax
//' @noRd
// [[Rcpp::export]]
NumericVector HHG_stats_cont(const NumericVector& x, const NumericVector& y,
                             int A, int mmax) {
  int nx=x.size(), ny=y.size(), N=nx+ny;
  std::vector< std::pair<double,int> > z(N);
  for(int i=0;i<nx;++i) z[i]=std::make_pair(x[i], 1);
  for(int i=0;i<ny;++i) z[i+nx]=std::make_pair(y[i], 0);
  std::sort(z.begin(), z.end());
  std::vector<int> b=atom_bounds(N, A);
  std::vector<double> cx(A+1), ct(A+1);
  cx[0]=0.0; ct[0]=0.0;
  for(int j=0;j<A;++j) {
    double c=0.0;
    for(int i=b[j];i<b[j+1];++i) c+=z[i].second;
    cx[j+1]=cx[j]+c;
    ct[j+1]=b[j+1];
  }
  return hhg_stats(cx, ct, A, mmax);
}

//' HHG statistics S_m for discrete data
//' @param x counts of first data set
//' @param y counts of second data set
//' @param mmax largest number of cells
//' @keywords internal
//' @return numeric vector of S_m, m=2,..,mmax
//' @noRd
// [[Rcpp::export]]
NumericVector HHG_stats_disc(const NumericVector& x, const NumericVector& y,
                             int mmax) {
  int A=x.size();
  std::vector<double> cx(A+1), ct(A+1);
  cx[0]=0.0; ct[0]=0.0;
  for(int j=0;j<A;++j) {
    cx[j+1]=cx[j]+x[j];
    ct[j+1]=ct[j]+x[j]+y[j];
  }
  return hhg_stats(cx, ct, A, std::min(mmax, A));
}

//' HHG test statistic for continuous data
//' @param x first data set
//' @param y second data set
//' @param A number of atoms
//' @param mmax largest number of cells
//' @param nullsorted matrix from HHG_null
//' @keywords internal
//' @return the test statistic -log(min p)
//' @noRd
// [[Rcpp::export]]
double HHG_cont(const NumericVector& x, const NumericVector& y,
                int A, int mmax, const NumericMatrix& nullsorted) {
  return minp_stat(HHG_stats_cont(x, y, A, mmax), nullsorted);
}

//' HHG test statistic for discrete data
//' @param x counts of first data set
//' @param y counts of second data set
//' @param mmax largest number of cells
//' @param nullsorted matrix from HHG_null
//' @keywords internal
//' @return the test statistic -log(min p)
//' @noRd
// [[Rcpp::export]]
double HHG_disc(const NumericVector& x, const NumericVector& y,
                int mmax, const NumericMatrix& nullsorted) {
  NumericVector S=HHG_stats_disc(x, y, mmax);
  return minp_stat(S, nullsorted);
}
