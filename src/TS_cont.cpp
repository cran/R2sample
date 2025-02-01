#include <Rcpp.h>
#include "repC.h"
#include "bincounter.h"
using namespace Rcpp;

//' find test statistics for continuous data
//' 
//' @param x first continuous data set
//' @param y second continuous data set
//' @keywords internal
//' @return A vector of test statistics
// [[Rcpp::export]]
NumericVector TS_cont(NumericVector x, NumericVector y) {
  
  CharacterVector methods = CharacterVector::create("KS", "Kuiper", "CvM", "AD", "LR", "ZA", "ZK", "ZC", "Wassp1");  
  int const nummethods=methods.size();
  
  int nx=x.size(),ny=y.size(),n=nx+ny, i, j;
  NumericVector xy(n), r(n), p(n);
  NumericVector TS(nummethods), Fx(n), Fy(n), w(n), px(n), py(n), sxy(n);
  IntegerVector Rx(nx),Ry(ny), D(n+2);
  double tmp;

  TS.names() =  methods;

  /*  sort data */
  
  std::sort(x.begin(), x.end());
  std::sort(y.begin(), y.end()); 


  /*  Data in one vector*/  
  
  for(i=0;i<nx;++i) xy[i]=x[i];
  for(i=0;i<ny;++i) xy[i+nx]=y[i];
  
  
  
  /* order(data)-1 and rank(data)-1 */
   
  for(i=0;i<n;++i) sxy[i]=xy[i]; 
  std::sort(sxy.begin(), sxy.end());
  Function Rorder("order");
  IntegerVector idx=Rorder(xy);
  Function Rrank("rank");
  IntegerVector R=Rrank(xy);
  idx=idx-1;
  R=R-1;
  for(i=0;i<nx;++i) Rx[i]=R[i]+1;
  std::sort(Rx.begin(), Rx.end());
  for(i=0;i<ny;++i) Ry[i]=R[i+nx]+1;
  std::sort(Ry.begin(), Ry.end());    
 

  /* empirical distribution functions of x and y evaluated on 
     combined data set*/
     
     for(int i=0;i<n;++i) {
       if(i==0) {
         Fx[0]=0;
         j=0;
       }  
       else 
         Fx[i]=Fx[i-1];  
       while ( (j<nx) && (x[j]<=sxy[i]) ) {
          Fx[i]=Fx[i]+1.0/nx;
          ++j;
       }
     }
     for(int i=0;i<n;++i) {
       if(i==0) {
         Fy[0]=0;
         j=0;
       }  
       else 
         Fy[i]=Fy[i-1];  
       while ( (j<ny) && (y[j]<=sxy[i]) ) {
          Fy[i]=Fy[i]+1.0/ny;
          ++j;
       }
     }        
  

  /*  Kolmogorov-Smirnov and Kuiper Tests*/  
    tmp = Fx[0]-Fy[0];
    double mx=-tmp, Mx=tmp;
    TS(0)=std::abs(mx);
    for(i=1;i<n;++i) {
      if(std::abs(Fx[i]-Fy[i])>TS(0)) TS(0)=std::abs(Fx[i]-Fy[i]);
      tmp = Fx[i]-Fy[i]; 
      if(mx>tmp) mx=tmp;
      if(Mx<tmp) Mx=tmp;      
    }
    if(-mx>Mx) TS(0)=-mx;
    else TS(0)=Mx;
    TS(1)=Mx-mx;

  /* Cramer-vonMises and Anderson-Darling test*/
   
       TS(2)=sum((Fx-Fy)*(Fx-Fy))*nx*ny/n/n;       

  /*    Anderson-Darling test*/  
  
       TS(3)=0.0;
       for(i=0;i<n-1;++i) {
         tmp=Fx[i]-Fy[i];
         TS(3)=TS(3)+tmp*tmp/(double(i)+1.0)/(double(n)-i);
       } 
       TS(3)=TS(3)*nx*ny;
 

  
  /*  Lehmann-Rosenblatt test*/  
  
      TS(4)=0.0;
      for(i=0;i<nx;++i) {
        tmp=double(Rx[i])-double(n)/double(nx)*(double(i)+1.0);
        TS(4)=TS(4)+ny*tmp*tmp;
      } 
      for(i=0;i<ny;++i) {
        tmp=double(Ry[i])-double(n)/double(ny)*(double(i)+1.0);          
        TS(4)=TS(4)+nx*tmp*tmp;
      }
      TS(4)=TS(4)/double(n)/nx/ny; 

  
  /*   Zhangs tests*/  
  
      D[0]=1;
      for(i=0;i<nx;++i) D[i+1]=Rx[i];
      D[nx+1]=n+1;
      int k1=0, k2=0;
      for(i=0;i<nx+1;++i) {
        if(D[i+1]!=D[i]) {
          for(j=0;j<D[i+1]-D[i];++j) {
           px[k1]=k2;
           k1=k1+1;
         }
       } 
       k2=k2+1;
      }
      for(i=0;i<nx;++i) px[Rx[i]-1]=px[Rx[i]-1]-0.5;
      for(i=0;i<n;++i) px[i]=px[i]/nx;
  
      D[0]=1;
      for(i=0;i<ny;++i) D[i+1]=Ry[i];
      D[ny+1]=n+1;
      k1=0, k2=0;
      for(i=0;i<ny+1;++i) {
        if(D[i+1]!=D[i]) {
          for(j=0;j<D[i+1]-D[i];++j) {
            py[k1]=k2;
            k1=k1+1;
          }
        }  
        k2=k2+1;
      }
      for(i=0;i<ny;++i) py[Ry[i]-1]=py[Ry[i]-1]-0.5;
      for(i=0;i<n;++i) py[i]=py[i]/ny;
      
      TS(5)=0.0;
      for(i=0;i<n;++i) {
        TS(5)=TS(5)+
          (nx*(px(i)*log(px(i)+1e-10)+(1-px(i))*log(1-px(i)+1e-10))+
          ny*(py(i)*log(py(i)+1e-10)+(1-py(i))*log(1-py(i)+1e-10)))/((i+0.5)*(n-i-0.5));
      }

      TS(6)=0;
      for(i=0;i<n;++i) {
            p(i)=(i+0.5)/n;
            tmp=nx*(px(i)*log(px(i)+1e-10)+(1-px(i))*log(1-px(i)+1e-10))+
              ny*(py(i)*log(py(i)+1e-10)+(1-py(i))*log(1-py(i)+1e-10))-
              n*(p(i)*log(p(i))+(1-p(i))*log(1-p(i)));  
            if(tmp>TS(6))  TS(6)=tmp;
      }

      tmp=0;
      for(i=0;i<nx;++i)  tmp=tmp+log(nx/(i+0.5)-1)*log(n/(Rx[i]-0.5)-1);
      for(i=0;i<ny;++i)  tmp=tmp+log(ny/(i+0.5)-1)*log(n/(Ry[i]-0.5)-1);  
      TS(7)=-tmp/n;

 /* Wasserstein p=1*/

    NumericVector cux1(nx+1), cuy1(ny+1),cux2(nx-1), cuy2(ny-1), uu(nx+ny-2);
    NumericVector xx(nx+ny+1), yy(nx+ny+1);
    IntegerVector xrepC(nx), yrepC(ny);
    if(nx==ny) {
      for(i=0;i<nx;++i) {
        TS(8)=TS(8)+std::abs(x(i)-y(i));
      }
      TS(8)=TS(8)/nx;
    }
    else {
      cux1(0)=0.0;
      cux1(1)=1.0/nx;
      cux2(0)=1.0/nx;
      for(i=2;i<nx;++i) {
        cux1(i)=cux1(i-1)+1.0/nx;
        cux2(i-1)=cux1(i); 
      }  
      cux1(nx)=1.0;
      cuy1(0)=0.0;
      cuy1(1)=1.0/ny;
      cuy2(0)=1.0/ny;
      for(i=2;i<ny;++i) {
        cuy1(i)=cuy1(i-1)+1.0/ny;
        cuy2(i-1)=cuy1(i); 
      }  
      cuy1(ny)=1.0;
      xrepC=bincounter(cuy2, cux1)+1;
      yrepC=bincounter(cux2, cuy1)+1;
      xx=repC(x, xrepC);
      yy=repC(y, yrepC);

      for(i=0;i<nx-1;++i) uu(i)=cux2(i);  
      for(i=0;i<ny-1;++i) uu(nx+i-1)=cuy2(i);
      std::sort(uu.begin(), uu.end());
      TS(8) = uu(0)*std::abs(xx(0)-yy(0));
      for(i=1;i<nx+ny-2;++i) {
        TS(8) = TS(8)+(uu(i)-uu(i-1))*std::abs(xx(i)-yy(i));
      }    
      TS(8) = TS(8)+(1.0-uu(nx+ny-3))*std::abs(xx(nx+ny-2)-yy(nx+ny-2));
    
   }

  return TS;
}

