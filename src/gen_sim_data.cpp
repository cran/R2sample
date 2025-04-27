#include <Rcpp.h>
#include "gen_sim_data.h"
using namespace Rcpp; 

//' simulate continuous data without weights
//' @param dta data set
//' @param TSextra extra stuff
//' @keywords internal
//' @return A list of permuted vectors
// [[Rcpp::export]]
Rcpp::List gen_sim_data(List dta, List TSextra) {
   List newdta; 
   if(dta.size()==2) {
     NumericVector x=dta["x"];
     NumericVector y=dta["y"];
     newdta=gen_cont_noweights(x, y, TSextra);
   }
   if(dta.size()==3) {
     IntegerVector x=dta["x"];
     IntegerVector y=dta["y"];
     NumericVector vals=dta["vals"];
     newdta=gen_disc(x, y, vals, TSextra);
   }
   if(dta.size()==4) {
     NumericVector x=dta["x"];
     NumericVector y=dta["y"];
     NumericVector wx=dta["wx"];
     NumericVector wy=dta["wy"];
     newdta=gen_cont_weights(x, y, wx, wy, TSextra);
   }
   return newdta;
} 
                               
//' simulate continuous data without weights
//' @param x first data set
//' @param y second data set
//' @param TSextra extra stuff
//' @keywords internal
//' @return A list of permuted vectors
// [[Rcpp::export]]
Rcpp::List gen_cont_noweights(NumericVector x, 
                     NumericVector y, List TSextra) {
  CharacterVector A= CharacterVector::create("rnull");
  CharacterVector B=TSextra.names();
  if(in(A, B)[0]) { 
    List simdta=List::create(Named("x")=x,Named("y")=y);
    Function rnull=TSextra["rnull"];
    List out=rnull(simdta);
    return out;
  }
  int nx, ny, n, i;
  nx=x.size();
  ny=y.size();
  n=nx+ny;
  NumericVector permx(nx), permy(ny), z(n);
  IntegerVector Index(n);
  List permdta=List::create(Named("x")=x,Named("y")=y);
  for(i=0;i<nx;++i) {
      Index(i)=i;
      z[i]=x(i);
  }  
  for(i=0;i<ny;++i) {
      Index(i+nx)=i+nx;
      z[i+nx]=y(i);
  } 
  Index = Rcpp::sample(Index, n);
  for(int i=0;i<nx;++i) permx[i]=z[Index[i]];
  for(i=0;i<ny;++i) permy[i]=z[Index[i+nx]];
  permdta["x"]=permx;
  permdta["y"]=permy;
  return permdta;
}

//' simulate continuous data with weights
//' @param x first data set
//' @param y second data set
//' @param wx weights of first data set
//' @param wy weights of second data set
//' @param TSextra extra stuff
//' @keywords internal
//' @return A list of permuted vectors
// [[Rcpp::export]]
Rcpp::List gen_cont_weights(NumericVector x, 
                     NumericVector y,
                     NumericVector wx,
                     NumericVector wy,
                     List TSextra) {
  
  CharacterVector A= CharacterVector::create("rnull");
  CharacterVector B=TSextra.names();
  if(in(A, B)[0]) { 
    List simdta=List::create(Named("x")=x,Named("y")=y,
                             Named("wx")=wx,Named("wy")=wy);
    Function rnull=TSextra["rnull"];
    List out=rnull(simdta);
    return out;
  }
  List permdta=List::create(Named("x")=x,Named("y")=y,
                            Named("wx")=wx,Named("wy")=wy);
  int nx, ny, n, i;
  nx=x.size();
  ny=y.size();
  n=nx+ny;
  NumericVector permx(nx),permwx(nx), permy(ny), permwy(ny),  
                z(n),  wz(n);
  IntegerVector Index(n);
  for(i=0;i<nx;++i) {
      Index(i)=i;
      z[i]=x(i);
      wz[i]=wx[i];
  }  
  for(i=0;i<ny;++i) {
      Index(i+nx)=i+nx;
      z[i+nx]=y(i);
      wz[i+nx]=wy[i];
  } 
  Index = Rcpp::sample(Index, n); 
  for(int i=0;i<nx;++i) {
      permx[i]=z[Index[i]];
      permwx[i]=wz[Index[i]]; 
  }  
  for(i=0;i<ny;++i) {
      permy[i]=z[Index[i+nx]];
      permwy[i]=wz[Index[i+nx]];
  }
  permdta["x"]=permx;
  permdta["y"]=permy;
  permdta["wx"]=permwx;
  permdta["wy"]=permwy;
  return permdta;
}

//' simulate new discrete data
//' @param dtax first data set, counts
//' @param dtay second data set, counts
//' @param vals values of discrete random variable
//' @param TSextra extra stuff
//' @keywords internal
//' @return A list of permuted vectors
// [[Rcpp::export]]
Rcpp::List gen_disc(IntegerVector dtax, 
                    IntegerVector dtay,
                    NumericVector vals,
                    List TSextra) {
  
  CharacterVector A= CharacterVector::create("rnull");
  CharacterVector B=TSextra.names();
  if(in(A, B)[0]) { 
    List simdta=List::create(Named("x")=dtax,
                             Named("y")=dtay,
                             Named("vals")=vals);
    Function rnull=TSextra["rnull"];
    List out=rnull(simdta);
    return out;
  }
  List out;
  int nx, ny, n, i, j,  d, k=dtax.size();
  NumericVector xy(k), xout(k), yout(k);
  IntegerVector Index(k);
  int samplingmethod=TSextra["samplingmethod"];
    
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
      List dta=List::create(Named("dtax")=dtax,
                            Named("dtay")=dtay,
                            Named("vals")=vals);
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
  out = List::create(Named("x") = xout, 
                       Named("y") = yout,
                       Named("vals") = vals);                    
  return out;
}

