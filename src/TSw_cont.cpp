#include <Rcpp.h>
using namespace Rcpp;
#include "Cpporder.h"

//' find test statistics for continuous data with weights
//' 
//' @param x first continuous data set
//' @param y second continuous data set
//' @param wx weights of x
//' @param wy weights of y
//' @keywords internal
//' @return A vector of test statistics
// [[Rcpp::export]]

NumericVector TSw_cont(std::vector<double>& x, 
                       std::vector<double>& y, 
                       std::vector<double>& wx, 
                       std::vector<double>& wy) {
  
  CharacterVector methods = CharacterVector::create("KS", "Kuiper", "CvM", "AD");  
  int const nummethods=methods.size();
  
  int nx=x.size(),ny=y.size(),n=nx+ny, i;
  std::vector<double> xy(n), cwx(nx), cwy(ny), cw(n), wxy(n), cwxy(n), Ixy(n);
  NumericVector TS(nummethods);
  double tmp1=0.0, tmp2=0.0, m=0.0, M=0.0;
  
  TS.names() =  methods;
  
  /*  create index vector for sorting x, wx by x and y, wy by y*/
  
  wx = Cpporder(wx, x);
  std::sort(x.begin(), x.end());
  tmp1=wx[0];
  cwx[0]=wx[0];
  for(i=1;i<nx;++i) {
    tmp1=tmp1+wx[i];
    cwx[i] = cwx[i-1]+wx[i];
  }  
  for(i=0;i<nx;++i) cwx[i]=cwx[i]/tmp1;
  wy=Cpporder(wy, y);
  std::sort(y.begin(), y.end());
  tmp1=wy[0];
  cwy[0]=wy[0];
  for(i=1;i<ny;++i) {
    tmp1=tmp1+wy[i];
    cwy[i] = cwy[i-1]+wy[i];
  }  
  for(i=0;i<ny;++i) cwy[i]=cwy[i]/tmp1;
  for(i=0;i<nx;++i) cw[i]=cwx[i];
  for(i=0;i<ny;++i) cw[nx+i]=cwy[i];
  
  /*join vectors x,y,  create index vector for sorting */  
  for(i=0;i<nx;++i) {
    xy[i]=x[i];
    wxy[i]=wx[i];
    Ixy[i]=1.0;
  }  
  for(i=0;i<ny;++i) {
    xy[nx+i]=y[i];
    wxy[nx+i]=wy[i];
    Ixy[nx+i]=2.0;
  }
  Ixy = Cpporder(Ixy, xy);
  cw=Cpporder(cw, xy);
  wxy=Cpporder(wxy, xy);
  cwxy[0]=wxy[0];
  for(i=1;i<nx;++i) cwxy[i]=cwxy[i-1]+wxy[i];
  for(i=0;i<ny;++i) cwxy[nx+i]=cwxy[nx+i-1]+wxy[i];
  tmp1=0.0;
  tmp2=0.0;
  for(i=0;i<n;++i) {
    if(Ixy[i]==2.0) tmp1=cw[i];
    else tmp2=cw[i];
    if(tmp1-tmp2>m) m=tmp1-tmp2;    
    if(tmp1<tmp2 && tmp1-tmp2<M) M=tmp1-tmp2;
    double tmp=(tmp1-tmp2)*(tmp1-tmp2);
    TS(2)=TS(2)+tmp;
    TS(3)=TS(3)+tmp/cwxy[i]/cwxy[n-i-1];
  }
  if(std::abs(m)>std::abs(M)) TS(0)=std::abs(m);
  else TS(0)=std::abs(M);
  TS(1) = std::abs(m)+std::abs(M);
  TS(2) = TS(2)*nx*ny/n/n;
  TS(3) = TS(3)*nx*ny;
  
  return TS;
  
}



