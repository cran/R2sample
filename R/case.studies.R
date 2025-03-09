#' This function creates the functions needed to run the various case studies.
#' @param which name of the case study.
#' @param nsample =500, sample size.
#' @return a list of functions

case.studies=function(which, nsample=500) {

  if(which=="uniform.linear.cont") {
     return(list(
       f=function(slope) {
        if(slope==0) y=runif(nsample)
        else y=(slope-1+sqrt((1-slope)^2+4*slope*runif(nsample)))/2/slope
        list(x=runif(nsample), y=y)
       },
       param_alt=round(seq(0, 0.5, length=25), 3)
     ))
   } 
    
   if(which=="uniform.linear.disc") {
     return(list(
       f = function(slope) {
         if(slope==0) p=rep(1/50, 50)
         else p=diff(slope * (0:50/50)^2 + (1 - slope) * 0:50/50)
         list(
           x=c(rmultinom(1, nsample, rep(1/50,50))),
           y=c(rmultinom(1, nsample, p)), 
           vals=0:49/50+1/100 # midpoints of intervals
         )
      },   
      param_alt=round(seq(0, 0.5, length=25), 3)
    ))
  }     
  
  if(which=="uniform.quadratic.cont") {
     return(list(
       f = function(a) {
          rquad=function(n,a) {
            if(a==0) return(runif(n))
            y=rep(0,n)
            for(i in 1:n) {
             repeat {
              x=runif(1)
              if(runif(1)<(1+a*(x-0.5)^2)/(1+a/4)) {
                 y[i]=x
                 break
              }
             } 
           }
           y
          }
          list(x=rquad(nsample, 0), y=rquad(nsample, a))
       },
       param_alt=round(seq(0, 6, length=25), 2)
    ))
  }    

  if(which=="uniform.quadratic.disc") {
     return(list(
       f =  function(a) {
                if(a==0) p=rep(1/100, 100)
                else p=diff( (24*0:100/100+a+8*a*((0:100/100)-1/2)^3)/(24+2*a) )
                list(
                   x=c(rmultinom(1, nsample, rep(1/100,100))),
                   y=c(rmultinom(1, nsample, p)), 
                   vals=0:99/100+1/200
                )
            },
            param_alt=round(seq(0, 6, length=25), 2)
    ))
  }          
  if(which=="uniform.bump.cont") {
     return(list(
       f = function(alpha) {
          x=runif(nsample)
          y=c(runif(floor(nsample*(1-alpha))), 
          rnorm(ceiling(nsample*alpha), 0.5, 0.05))
          list(x=x, y=y)
       },
       param_alt=round(seq(0, 0.3, length=25), 3)
    ))
  }     

  if(which=="uniform.bump.disc") {
     return(list(
       f = function(alpha) {
         vals=(1:50)/51
         bins=0:50/50
         px=diff(punif(bins))
         py=diff( (1-alpha)*punif(bins) + alpha*pnorm(bins, 0.5, 0.05))
         x=c(rmultinom(1, nsample, px))
         y=c(rmultinom(1, nsample, py))
         I=c(1:50)[x+y>0]
         list(x=x[I],y=y[I],vals=vals[I])
       },
       param_alt=round(seq(0, 0.3, length=25), 3)
    ))
  }     

  if(which=="uniform.sine.cont") {
     return(list(
       f =  function(a) {
              if(a==0) y = runif(nsample)
              else {
                 y = NULL
                 repeat { 
                    z=runif(nsample)
                    I=ifelse(runif(nsample)<(1+a*sin(4*pi*z))/(1+a), TRUE, FALSE)
                    y = c(y, z[I]) 
                    if(length(y)>nsample) break
                 }
              }   
              list(x=runif(nsample), y=y[1:nsample])  
            },
            param_alt=round(seq(0, 0.7, length=25), 3)
    ))
  }           

  if(which=="uniform.sine.disc") {
     return(list(
       f =  function(a) {
         psin=function(x, a=0) x-a*(cos(4*pi*x)-1)/4/pi
         bins = 0:50/50
         vals = (bins[-1]+bins[-51])/2
         px=diff(punif(bins))
         py=diff(psin(bins, a))
         x=c(rmultinom(1, nsample, px))
         y=c(rmultinom(1, nsample, py))
         I=seq_along(vals)[x+y>0]
         list(x=x[I],y=y[I], vals=vals[I])
       },
       param_alt=round(seq(0, 0.7, length=25), 3)
    ))
  }    
   
  if(which=="beta22.betaaa.cont") {
     return(list(
       f =  function(a) {
         list(x=rbeta(nsample, 2, 2), y=rbeta(nsample, a, a))
       },
       param_alt=round(seq(2, 4, length=25), 3)
    ))
  }  

  if(which=="beta22.betaaa.disc") {
     return(list(
       f =  function(a) {
        bins=qbeta(0:50/50,2,2)
        px=diff(pbeta(bins, 2, 2))
        py=diff(pbeta(bins, a, a))
        x=c(rmultinom(1, nsample, px))
        y=c(rmultinom(1, nsample, py))
        vals=(bins[-1]+bins[-51])/2
        I=c(1:50)[x+y>0]
        list(x=x[I],y=y[I],vals=vals[I])
      },
      param_alt=round(seq(2, 4, length=25), 3)
    ))
  }
    
  if(which=="beta22.beta2a.cont") {
     return(list(
       f =  function(beta) {
              list(x=rbeta(nsample, 2, 2), y=rbeta(nsample, 2, beta))
            },
            param_alt=round(seq(2, 3, length=25), 3)
    ))
  }          

  if(which=="beta22.beta2a.disc") {
     return(list(
       f =  function(beta) {
        bins=qbeta(0:50/50,2,2)
        px=diff(pbeta(bins, 2, 2))
        py=diff(pbeta(bins, 2, beta))
        x=c(rmultinom(1, nsample, px))
        y=c(rmultinom(1, nsample, py))
        vals=(bins[-1]+bins[-51])/2
        I=c(1:50)[x+y>0]
        list(x=x[I],y=y[I],vals=vals[I])
      },
      param_alt=round(seq(2, 3, length=25), 3)
    ))
  } 
   

  if(which=="normal.shift.cont") {
     return(list(
       f =  function(mu) {list(x=rnorm(nsample), y=rnorm(nsample, mu))},
       param_alt=round(seq(0, 0.4, length=25), 3)
     ))
   }    
   
  if(which=="normal.shift.disc") {
     return(list(
       f =  function(mu) {
          bins=seq(-3, 3, length=51)
          vals=(bins[-1]+bins[-51])/2
          px=diff((pnorm(bins)-pnorm(-3))/(1-2*pnorm(-3)))
          py=diff((pnorm(bins,mu)-pnorm(-3,mu))/(1-2*pnorm(-3,mu)))
          x=c(rmultinom(1, nsample, px))
          y=c(rmultinom(1, nsample, py))
          I=c(1:50)[x+y>0]
          list(x=x[I],y=y[I],vals=vals[I])
      },
      param_alt=round(seq(0, 0.4, length=25), 3)
    ))
  }    

  if(which=="normal.stretch.cont") {
     return(list(
       f =  function(sigma) {list(x=rnorm(nsample), y=rnorm(nsample, 0, sigma))},
       param_alt=round(seq(1, 1.5, length=25), 3)
     ))
  }     

  if(which=="normal.stretch.disc") {
     return(list(
       f =  function(sigma) {
          bins=seq(-3, 3, length=51)
          vals=(bins[-1]+bins[-51])/2
          px=diff((pnorm(bins)-pnorm(-3))/(1-2*pnorm(-3)))
          py=diff((pnorm(bins, 0, sigma)-pnorm(-3, 0, sigma))/(1-2*pnorm(-3, 0, sigma)))
          x=c(rmultinom(1, nsample, px))
          y=c(rmultinom(1, nsample, py))
          I=c(1:50)[x+y>0]
          list(x=x[I],y=y[I], vals=vals[I])
       },
       param_alt=round(seq(1, 1.5, length=25), 3)
     ))
  }    
  
  if(which=="normal.t.cont") {
     return(list(
       f =  function(df) {
          y = rt(2*nsample, df)
          y = y[abs(y)<6][1:nsample]
          list(x=rnorm(nsample), y=y) 
       },
       param_alt=1:25
    ))
  }
  
  if(which=="normal.t.disc") {
     return(list(
       f =  function(df) {
          bins=seq(-5, 5, length=51)
          vals=(bins[-1]+bins[-51])/2
          px=diff((pnorm(bins)-pnorm(-5))/(1-2*pnorm(-5)))
          py=diff((pt(bins,df)-pt(-5, df))/(1-2*pt(-5, df)))
          x=c(rmultinom(1, nsample, px))
          y=c(rmultinom(1, nsample, py))
          I=c(1:50)[x+y>0]
          list(x=x[I],y=y[I], vals=vals[I])
       },
       param_alt=1:25
    ))
  }     

  if(which=="normal.outlier1.cont") {
     return(list(
       f =  function(alpha) {
          x = rnorm(nsample) 
          y = c(rnorm(floor((1-alpha)*nsample)), rnorm(ceiling(alpha*nsample), 5))
          list(x=x, y=y)
       },    
        param_alt=round(seq(0, 0.2, length=25), 3)
    ))
  }      

  if(which=="normal.outlier1.disc") {
     return(list(
       f =  function(alpha) {
          pnormtrunc = function(x, mu)
              (pnorm(bins, mu)-pnorm(-3, mu))/diff(pnorm(c(-3,8), mu))
          bins=seq(-3, 8, length=51)
          vals=(bins[-1]+bins[-51])/2
          px=diff(pnormtrunc(bins, 0))
          py=diff((1-alpha)*pnormtrunc(bins, 0)+alpha*pnormtrunc(bins, 5))
          x=c(rmultinom(1, nsample, px))
          y=c(rmultinom(1, nsample, py))
          I=c(1:50)[x+y>0]
          list(x=x[I],y=y[I], vals=vals[I])
      },
      param_alt=round(seq(0, 0.2, length=25), 3)
    ))
  }  
     

  if(which=="normal.outlier2.cont") {
     return(list(
       f =  function(alpha) {
          x=rnorm(nsample)
          m=round(alpha/2*nsample)
          y=c(rnorm(nsample-2*m), runif(m, -5, -3), runif(m, 3, 5))
          list(x=x, y=y) 
        },         
        param_alt=round(seq(0, 0.2, length=25), 3)
    ))
  }  

  if(which=="normal.outlier2.disc") {
     return(list(
       f =  function(alpha) {
          bins=seq(-5, 5, length=51)
          vals=(bins[-1]+bins[-51])/2
          px=diff( (pnorm(bins)-pnorm(-5))/(2*pnorm(5)-1) )
          py=diff( (1-alpha)*(pnorm(bins)-pnorm(-5))/(2*pnorm(5)-1) +
               alpha/2*punif(bins, -5, -3) + 
               alpha/2*punif(bins, 3, 5))
          x=c(rmultinom(1, nsample, px))
          y=c(rmultinom(1, nsample, py))
          I=c(1:50)[x+y>0]
          list(x=x[I],y=y[I], vals=vals[I])
      },
      param_alt=round(seq(0, 0.2, length=25), 3)
    ))
  }  
     

  if(which=="exponential.gamma.cont") {
     return(list(
       f =  function(b) {list(x=rexp(nsample, 1), y=rgamma(nsample, 1, b))},
        param_alt=round(seq(1, 1.5, length=25), 3)
     ))
  }
        
  if(which=="exponential.gamma.disc") {
     return(list(
       f =  function(b) {
          bins=seq(0, 4, length=51)
          vals=(bins[-1]+bins[-51])/2
          px=diff(pexp(bins, 1) )/pexp(4, 1)
          py=diff(pgamma(bins, 1, b) )/pgamma(4, 1, b)
          x=c(rmultinom(1, nsample, px))
          y=c(rmultinom(1, nsample, py))
          I=c(1:50)[x+y>0]
          list(x=x[I],y=y[I], vals=vals[I])
      },
      param_alt=round(seq(1, 1.5, length=25), 3)
    ))
  }
      

  if(which=="exponential.weibull.cont") {
     return(list(
       f =  function(b) {list(x=rexp(nsample, 1), y=rweibull(nsample, 1, b))},
       param_alt=round(seq(1, 1.6, length=25), 3)
     ))
  }
       
  if(which=="exponential.weibull.disc") {
     return(list(
       f =  function(b) {
         bins=seq(0, 6, length=51)
         vals=(bins[-1]+bins[-51])/2
         px=diff(pexp(bins, 1)/pexp(6, 1))
         py=diff(pweibull(bins, 1, b)/pweibull(6,1, b))
         x=c(rmultinom(1, nsample, px))
         y=c(rmultinom(1, nsample, py))
         I=c(1:50)[x+y>0]
         list(x=x[I],y=y[I], vals=vals[I])
      },
      param_alt=round(seq(1, 1.6, length=25), 3)
    ))
  }     

  if(which=="exponential.bump.cont") {
     return(list(
       f = function(alpha) {
          x=rexp(nsample, 1)
          y=c(rexp(floor(nsample*(1-alpha)), 1), 
              rnorm(ceiling(nsample*alpha), 0.5, 0.05))
          list(x=x, y=y)
      },   
      param_alt=round(seq(0, 0.3, length=25), 3)
    ))
  }     

  if(which=="exponential.bump.disc") {
     return(list(
       f = function(alpha) {
          bins=seq(0, 3, length=51)
          vals=(bins[-1]+bins[-51])/2
          px=diff(pexp(bins, 1)/pexp(3, 1))
          py=diff((1-alpha)*pexp(bins, 1)/pexp(3, 1) + 
               alpha*pnorm(bins, 0.5,0.05)/pnorm(3, 0.5, 0.05))
          x=c(rmultinom(1, nsample, px))
          y=c(rmultinom(1, nsample, py))
          I=c(1:50)[x+y>0]
          list(x=x[I],y=y[I],vals=vals[I])
      },
      param_alt=round(seq(0, 0.3, length=25), 3)
    ))
  }     

  if(which=="gamma.normal.cont") {
     return(list(
       f = function(mu) {
         x=rgamma(nsample, mu, 1)
         list(x=x, y=rnorm(nsample, mean(x), sd(x)))
       },
       param_alt=round(seq(0.5, 20,length=25), 2)
    ))
  }     

  if(which=="gamma.normal.disc") {
     return(list(
       f = function(mu) {
        x=rgamma(nsample, mu, 1)
        y=rnorm(nsample, mean(x), sd(x))
        bins=quantile(c(x,y),0:50/50)
        vals=quantile(c(x,y),1:50/51)
        list(
          x=hist(x, bins, plot=FALSE)$counts, 
          y=hist(y, bins, plot=FALSE)$counts,
          vals=vals)
      },
      param_alt=round(seq(0.5, 20,length=25), 2)
    ))
  } 
      
  if(which=="normal.normalmixture.cont") {
     return(list(
       f=function(mu) {
          x = rnorm(nsample)
          y = c(rnorm(250, -mu), rnorm(250, mu))
          list(x=x, y=y)
       },
       param_alt=round(seq(0, 1, length=25), 3)
    ))
  }      

  if(which=="normal.normalmixture.disc") {
     return(list(
       f = function(mu) {
          bins = seq(-4, 4, length=51)
          vals = (bins[-1]+bins[-51])/2
          px=diff((pnorm(bins)-pnorm(-4))/(2*pnorm(4)-1))
          py=diff((pnorm(bins, -mu)-pnorm(-4, -mu))/(2*pnorm(4,-mu)-1)/2 +
              (pnorm(bins, mu)-pnorm(-4, mu))/(2*pnorm(4, mu)-1)/2)
          x=c(rmultinom(1, nsample, px))
          y=c(rmultinom(1, nsample, py))
          I=seq_along(vals)[x+y>0]
          list(x=x[I],y=y[I], vals=vals[I]) 
      },
      param_alt=round(seq(0, 1, length=25), 3)
    ))
  }    
      
  if(which=="uniform.uniformmixture.cont") {
     return(list(
       f=function(alpha=0.5) {
          x = runif(nsample)
          z = runif(nsample)
          y = ifelse(z<1-alpha, z/2/(1-alpha), (z-1+alpha)/2/alpha+0.5)
          list(x=x, y=y)
      },
      param_alt=round(seq(0.5, 0.7, length=25), 3)
    ))
  }    
  
  if(which=="uniform.uniformmixture.disc") {
     return(list(
       f = function(alpha) {
          bins = seq(0, 1, length=51)
          vals = (bins[-1]+bins[-51])/2
          px=diff(punif(bins))
          py=diff(alpha*punif(bins, 0, 1/2) + (1-alpha)*punif(bins, 1/2, 1))
          x=c(rmultinom(1, nsample, px))
          y=c(rmultinom(1, nsample, py))
          I=seq_along(vals)[x+y>0]
          list(x=x[I],y=y[I], vals=vals[I]) 
       },
       param_alt=round(seq(0.5, 0.7, length=25), 3)
    ))
  }  
      
  if(which=="uniform.betamixture.cont") {
     return(list(
       f=function(alpha) {
          x = runif(nsample)
          n=round(nsample*(1-alpha))
          y = c(runif(n), rbeta(nsample-n, 2, 2))
          list(x=x, y=y)
      },
      param_alt=round(seq(0, 0.7, length=25), 3)
    ))
  }     

  if(which=="uniform.betamixture.disc") {
     return(list(
       f = function(alpha) {
          bins = seq(0, 1, length=51)
          vals = (bins[-1]+bins[-51])/2
          px=diff(punif(bins))
          py=diff((1-alpha)*punif(bins) + alpha*pbeta(bins, 2, 2))
          x=c(rmultinom(1, nsample, px))
          y=c(rmultinom(1, nsample, py))
          I=seq_along(vals)[x+y>0]
          list(x=x[I],y=y[I], vals=vals[I]) 
      },
      param_alt=round(seq(0, 0.7, length=25), 3)
    ))
  } 
     
  if(which=="chisquare.noncentral.cont") {
     return(list(
       f = function(ncp) {list(x=rchisq(nsample, 5), y=rchisq(nsample, 5, ncp))},
       param_alt=round(seq(0, 1.5, length=25), 3)
     ))
  }     

  if(which=="chisquare.noncentral.disc") {
     return(list(
       f = function(ncp) {
          bins = seq(0, 20, length=51)
          vals = (bins[-1]+bins[-51])/2
          px=diff(pchisq(bins, 5)/pchisq(20, 5))
          py=diff(pchisq(bins, 5, ncp)/pchisq(20, 5, ncp))
          x=c(rmultinom(1, nsample, px))
          y=c(rmultinom(1, nsample, py))
          I=seq_along(vals)[x+y>0]
          list(x=x[I],y=y[I], vals=vals[I]) 
      },
      param_alt=round(seq(0, 1.5, length=25), 3)
    ))
  }    

  if(which=="uniform.triangular.cont") {
     return(list(
       f = function(a) {
          x = runif(nsample)
          y = runif(nsample)
          if(a>0) y = ifelse(y<1/2, (a-1+sqrt((1-a)^2+8*a*y))/(4*a),
                (1+3*a-sqrt((1+3*a)^2-8*a*(a+y)))/(4*a))
          list(x=x, y=y)
      },
      param_alt=round(seq(0, 0.6, length=25), 3)
    ))
  }    

  if(which=="uniform.triangular.disc") {
     return(list(
       f = function(a) {
          ptriangular = function(x, a) {
            if(a==0) return(punif(x))
            ifelse(x<0.5, (1-a)*x+2*a*x^2, -2*a*x^2+(1+3*a)*x-a)
          }
          bins = seq(0, 1, length=51)
          vals = (bins[-1]+bins[-51])/2
          px=diff(punif(bins))
          py=diff(ptriangular(bins, a))
          x=c(rmultinom(1, nsample, px))
          y=c(rmultinom(1, nsample, py))
          I=seq_along(vals)[x+y>0]
          list(x=x[I],y=y[I], vals=vals[I]) 
      },
      param_alt=round(seq(0, 0.8, length=25), 3)
    ))
  } 
     
}
