#include <Rcpp.h>
using namespace Rcpp;

//' find test statistics for discrete data
//' 
//' @param dta A list
//' @param ADweights A vector of weights for AD method
//' @param doMethod A character vector of methods to include
//' @return A vector of test statistics
// [[Rcpp::export]]

NumericVector TS_disc_cpp(List dta, 
                NumericVector ADweights,
                CharacterVector doMethod= CharacterVector::create("chi large", 
                "chi small", "t test", "KS", "Kuiper", "CvM", "AD", "LR", "ZA", "ZK", "ZC", "Wassp1"))  {
  
  int const nummethods=10;
  NumericVector vals = as<NumericVector>(dta["vals"]);    
  IntegerVector x = as<IntegerVector>(dta["x"]);
  IntegerVector y = as<IntegerVector>(dta["y"]);
  int k=vals.size(), nx, ny, n, i, j;
  NumericVector TS(nummethods), R(k), p(k+1), z(k+1), p1(k), z1(k);
  NumericVector tmpx(2*k+1), tmpy(2*k+1), tmpxy(2*k+1), w(2*k+1), Fx(k), Fy(k);  
  IntegerVector xy(k), cx(k), cy(k), sxy(k), cxy(k), D(k+1), D1(k+1);

  double tmp, sx, sy, a;
  

  TS.names() = CharacterVector::create("t test", "KS", "Kuiper", "CvM", "AD", "LR", "ZA", "ZK", "ZC", "Wassp1");  
 
  /*  sample sizes*/  
  
  nx=0;
  ny=0;
  for(i=0;i<k;++i) {
    nx+=x[i];
    ny+=y[i];
  }   
  n=nx+ny;

  /*  combined data, then sorted, cumulative sums of x, y and xy */
  
  for(i=0;i<k;++i) {
    if(i==0) cx[i]=x[i];
    else cx[i]=cx[i-1]+x[i];
    if(i==0) cy[i]=y[i];
    else cy[i]=cy[i-1]+y[i];    
    xy[i]=x[i]+y[i];
    sxy[i]=xy[i];
    if(i==0) cxy[i]=xy[i];
    else cxy[i]=cxy[i-1]+xy[i];
  }  
  std::sort(sxy.begin(), sxy.end());

    /*  Ranks */
  
    for(j=0;j<k;++j) {
      R[j]=0.0;
      for(i=1;i<=xy[j];++i) {
        R[j]=R[j]+i;
      }
      R[j]=R[j]/xy[j];
      if(j>0) R[j]=R[j]+cxy[j-1];
    }
    
    
   /*  t test  */  

   if(in(CharacterVector::create("t test"), doMethod)[0]==TRUE) {
     sx=0.0; 
     sy=0.0;
     for(i=0;i<k;++i) {
       sx+=vals[i]*x[i];
       sy+=vals[i]*y[i];
     }  
     TS(0)=std::abs(sx/nx-sy/ny);
   }

  /*  Kolomogorov-Smirnov and Kuiper */
  
 
  LogicalVector H=in(CharacterVector::create("KS", "Kuiper"), doMethod);
  if( (H[0]==TRUE) || (H[1]==TRUE) ) {
     double mx=double(cx[0])/nx-double(cy[0])/ny;
     double Mx=double(cx[0])/nx-double(cy[0])/ny;
     for(i=1;i<k;++i) {
       tmp=double(cx[i])/nx-double(cy[i])/ny;
       if(mx>tmp) mx=tmp;
       if(Mx<tmp) Mx=tmp;
     }
     if(H[0]==TRUE) {
         if(-mx>Mx) TS(1)=-mx;
         else TS(1)=Mx;
     }
     if(H[1]==TRUE) TS(2)=Mx-mx;
   }

  /* Cramer-vonMises and Anderson-Darling test*/
   
   H=in(CharacterVector::create("CvM", "AD"), doMethod);
   if( (H[0]==TRUE) || (H[1]==TRUE) ) {
     Fx(0)=x(0)/double(nx);
     Fy(0)=y(0)/double(ny);
     for(i=1;i<k;++i) {
        Fx(i)=Fx(i-1)+x(i)/double(nx);
        Fy(i)=Fy(i-1)+y(i)/double(ny);         
     } 
     
     if(H[0]==TRUE) {
        TS(3)=0.0;
        for(i=0;i<k;++i) {
          tmp=Fx(i)-Fy(i);     
          TS(3)=TS(3)+xy(i)*tmp*tmp;
        }
        TS(3)=TS(3)*nx*ny/n/n;    
     }
   
     if(H[1]==TRUE) {
        TS(4)=0.0;
        for(i=0;i<k;++i) {
          tmp=Fx(i)-Fy(i);     
          TS(4)=TS(4)+ADweights(i)*tmp*tmp;
        }
        TS(4)=TS(4)*nx*ny; 
     } 
   }

  /*    Lehmann-Rosenblatt test*/   
   
  
  if(in(CharacterVector::create("LR"), doMethod)[0]==TRUE) {
    double xt=double(n)/double(nx),yt=double(n)/double(ny);
    TS(5)=double(n)*n*( (nx+1.0)*ny/nx*(2*nx+1.0)+(ny+1.0)*nx/ny*(2*ny+1.0) )/6.0;
    for(i=0;i<k;++i) {
      TS(5)=TS(5)+
        double(ny)*(x[i]*R[i]*R[i]-xt*(cx[i]*(cx[i]+1.0)-cx[i-1]*(cx[i-1]+1.0))*R[i])+
        double(nx)*(y[i]*R[i]*R[i]-yt*(cy[i]*(cy[i]+1.0)-cy[i-1]*(cy[i-1]+1.0))*R[i]);
    }  
    TS(5)=TS(5)/n/nx/ny;
  }  

    /* Zhang's tests */
  
  H=in(CharacterVector::create("ZA", "ZK", "ZC"), doMethod);
  if( (H[0]==TRUE) || (H[1]==TRUE) || (H[2]==TRUE) ) {    D[0]=ceil(R[0])-1;
    for(i=0;i<k-1;++i) D[i+1]=ceil(R[i+1])-ceil(R[i]);
    D[k]=n+1-ceil(R[k-1]);
    p[0]=0.0;
    for(i=1;i<k+1;++i) p[i]=double(cx[i-1])/nx;
    for(i=0;i<k+1;++i) 
       z[i]=nx*(p[i]*log(p[i]+1e-10)+(1-p[i])*log(1-p[i]+1e-10));
    for(i=0;i<k;++i) {
       p1[i]=(double(cx[i])-0.5)/nx;
       if(p1[i]<0) p1[i]=0;  
    }   
    for(i=0;i<k;++i) 
      z1[i]=nx*(p1[i]*log(p1[i]+1e-10)+(1-p1[i])*log(1-p1[i]+1e-10));
    tmpx[0]=z[0];
    for(i=0;i<k;++i) {
      tmpx[2*i+1]=z1[i];
      tmpx[2*i+2]=z[i+1];
    }
    p[0]=0.0;
    for(i=1;i<k+1;++i) p[i]=double(cy[i-1])/ny;
    for(i=0;i<k+1;++i) 
      z[i]=ny*(p[i]*log(p[i]+1e-10)+(1-p[i])*log(1-p[i]+1e-10));
    for(i=0;i<k;++i) {
       p1[i]=(double(cy[i])-0.5)/ny;
       if(p1[i]<0) p1[i]=0;
    }   
    for(i=0;i<k;++i) 
      z1[i]=ny*(p1[i]*log(p1[i]+1e-10)+(1-p1[i])*log(1-p1[i]+1e-10));
    tmpy[0]=z[0];
    for(i=0;i<k;++i) {
      tmpy[2*i+1]=z1[i];
      tmpy[2*i+2]=z[i+1];
    }
    for(i=0;i<2*k+1;++i) tmpxy[i]=tmpx[i]+tmpy[i];
    D1[0]=D[0];
    for(i=1;i<k+1;++i) D1[i]=D[i]+D1[i-1];

     w[0]=0.0;
     for(j=1;j<=D1[0];++j) w[0]=w[0]+1.0/( (j-0.5)*(n-j+0.5) );
     for(i=1;i<=k;++i) {
        w[2*i-1]=1.0/( (D1[i-1]+0.5)*(n-D1[i-1]-0.5) );
        w[2*i]=0.0;
        for(j=(D1[i-1]+2);j<=D1[i];++j) 
            w[2*i]=w[2*i]+1.0/( (j-0.5)*(n-j+0.5) );
     } 
     
     if(H[0]==TRUE) {    
        TS(6)=0.0;
        for(i=0;i<2*k+1;++i) TS(6)=TS(6)+tmpxy[i]*w[i];
     }
     
     if(H[1]==TRUE) {
      a=(D1[0]-0.5)/n;
      TS(7)=tmpxy[0]-n*(a*log(a)+(1-a)*log(1-a));
      for(i=1;i<=k;++i) {
        a=(D1[i]+0.5)/n;
        tmp=tmpxy[2*i]-n*(a*log(a)+(1-a)*log(1-a));
        if(tmp>TS(7)) TS(7)=tmp;
        a=(D1[i]+1.5)/n;
        tmp=tmpxy[2*i+1]-n*(a*log(a)+(1-a)*log(1-a));
        if(tmp>TS(7)) TS(7)=tmp;
        if(i<k) {
          a=(D1[i+1]-0.5)/n;
          tmp=tmpxy[2*i+1]-n*(a*log(a)+(1-a)*log(1-a));
          if(tmp>TS(7)) TS(7)=tmp;
        }  
      }
     }     
     
     if(H[2]==TRUE) {
      TS(8)=0.0;
      j=0;
      for(i=1;i<=nx;++i) {
        TS(8)=TS(8)+log(nx/(i-0.5)-1)*log( n/(R[j]-0.5)-1 );
        if(i==cx[j]) ++j;
      }
      j=0;
      for(i=1;i<=ny;++i) {
       TS(8)=TS(8)+log(ny/(i-0.5)-1)*log( n/(R[j]-0.5)-1 );
       if(i==cy[j]) ++j;
      }
      TS(8)=-TS(8)/n;
    }
  } 
 
 /*  Wasserstein */

  if(in(CharacterVector::create("Wassp1"), doMethod)[0]==TRUE) {
     TS(9)=0.0;  
     for(i=0;i<k-1;++i) {
       TS(9)=TS(9)+std::abs(cx(i)/double(nx)-cy(i)/double(ny))*(vals(i+1)-vals(i));
     }
  }  
  return TS;
}
