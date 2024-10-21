#include <Rcpp.h>
#include <ctime>
using namespace Rcpp;

//' run permutation test.
//' 
//' @param x a numeric vector.
//' @param y a numeric vector.
//' @param TS routine to calculate test statistics for non-chi-square tests 
//' @param typeTS type of a test statistic
//' @param TSextra additional info passed to TS, if necessary
//' @param wx a numeric vector of weights of x.
//' @param wy a numeric vector of weights of y.
//' @param B =5000, number of simulation runs.
//' @keywords internal
//' @return A list with test statistics and p values
// [[Rcpp::export]]
List perm_test_cont(NumericVector x, 
                   NumericVector y, 
                   Function TS,
                   int typeTS,
                   List TSextra,
                   NumericVector wx, 
                   NumericVector wy,
                   int B=5000) { 
               
  NumericVector TS_data;
  NumericVector permx(x.size()),permwx(wx.size()),
                permy(y.size()), permwy(wy.size());
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  int seed=round(100000*R::runif(0, 1));
  set_seed_r(seed);
  if(typeTS==1) TS_data=TS(x, y); 
  if(typeTS==2) TS_data=TS(x, y, wx, wy); 
  if(typeTS==3) TS_data=TS(x, y);   
  if(typeTS==4) TS_data=TS(x, y, TSextra); 
  int const nummethods=TS_data.size();
  int i, j;
  NumericVector pvals(nummethods), TS_perm(nummethods), stats(nummethods);
  List dta=List::create(Named("x") = x, Named("y") = y,
                        Named("wx") = wx, Named("wy") = wy);
  CharacterVector TSmethods=TS_data.names();
  stats.names() = TSmethods; 
  pvals.names() = TSmethods;

/*  Find test statistics for data */      
  
  for(i=0;i<nummethods;++i) {
    stats(i)=TS_data(i);  
    pvals(i)=0.0;  
  }

/*  if B=0 return just the test statistics */    
  if(B==0) return List::create(Named("statistics")=stats);
  
/* run permutation test  */   
  int nx=x.size(), ny=y.size(), n;
  NumericVector z(nx+ny), wz(nx+ny);
  IntegerVector Index(nx+ny);
  n=nx+ny;
  for(int i=0;i<nx;++i) {
    Index(i)=i;
    z[i]=x(i);
    wz[i]=wx[i];
  }  
  for(int i=0;i<ny;++i) {
    Index(i+nx)=i+nx;
    z[i+nx]=y(i);
    wz[i+nx]=wy[i];
  } 
  for(j=0;j<B;++j) {
    Index = Rcpp::sample(Index, n); 
    for(int i=0;i<nx;++i) {
      permx[i]=z[Index[i]];
      permwx[i]=wz[Index[i]]; 
    }  
    for(int i=0;i<ny;++i) {
      permy[i]=z[Index[i+nx]];
      permwy[i]=wz[Index[i+nx]];
    }
    if(typeTS==1) TS_perm=TS(permx, permy); 
    if(typeTS==2) TS_perm=TS(permx, permy, permwx, permwy); 
    if(typeTS==3) TS_perm=TS(permx, permy);   
    if(typeTS==4) TS_perm=TS(permx, permy, TSextra); 
    for(i=0;i<nummethods;++i) {
       if(TS_data(i)<TS_perm(i)) pvals(i)=pvals(i)+1.0;
    }
  }

/* find p values and return list with test statistics and p values */  
  for(i=0;i<nummethods;++i) pvals(i)=pvals(i)/B;
  return  List::create(Named("statistics")=stats, Named("p.values")=pvals);
  
}

