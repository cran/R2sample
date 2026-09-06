#' Create case studies
#' 
#' This function creates the functions needed to run the various case studies.
#' That means it generates two data sets the are equal with respect
#' to three out of these four: mean, variance, skewness and kurtosis, and differ
#' in the other. If which_diff=mixed, all are different.
#'
#' @param which name of the case study.
#' @param which_diff ="mean", name of parameter in which data sets differ.
#' @param data_type ="cont", is data continuous or discrete?
#' @param n =500, sample sizes.
#' @param nbins =50, number of bins for discrete data
#' @param funs_list list of alternative functions
#' @param ReturnCaseNames =FALSE, should list of case studies be returned?
#' @return a function that generates data sets
#' @export 
case_studies=function(which, which_diff ="mean", data_type="cont",
                n=500, nbins=50, funs_list, ReturnCaseNames =FALSE) {
# check if old version is assumed
  if(length(which_diff)==1 &&
     is.numeric(which_diff) &&
     !(which_diff %in% 1:6)) {
     stop("the arguments list of the function case.studies has changed.
          Please check the help file or the vignette",
          call.=FALSE)
  }
  list.of.studies=c(
    "uniform", "linear", "quadratic", "uniform-betabump", "sine", "betaaa",
    "beta2a", "triangular","uniformmixture", "betamixture", "normal-pure",
    "normal-mix1",  "normal-mix2", "normal-mix3", "normal-mix4", "gamma1", 
    "gamma2", "gammamixture",  "exponential-mixture", "noncentral",
    "uniform-linear", "uniform-quadratic", "uniform-bump",
    "uniform-sine",  "beta22.betaaa", "beta22.beta2a", 
    "uniform.uniformmixture", "uniform.betamixture",
    "uniform.triangular", "double.exponential", "normal.onestretch",
    "normal.t", "normal.outlier1", 
    "normal.outlier2", "normal.normalmixture",
    "exponential.gamma", "exponential.weibull", 
    "exponential.bump", "chisquare.noncentral",  
    "gamma.gammamixture"
  )
  if(ReturnCaseNames) return(list.of.studies)
  if(is.numeric(which[1])) which=list.of.studies[which]
  alltypes=c("mean", "variance", "skewness", "kurtosis",
             "mixed", "Null")
  if(is.numeric(which_diff[1]))
      which_diff=alltypes[which_diff]
  
  if(which %in% list.of.studies[1:20]) {
    if(missing(funs_list)) {
      funs = R2sample::funs_list[[data_type]][[which]][[which_diff]]
    } else {
      funs = funs_list
    }
    ralt = function(n) funs$random(n)
    palt = function(x) funs$cdf(x)
  }
  if(data_type=="disc") {
     if(which%in% list.of.studies[c(1:10, 21:29)]) {
       bins=0:nbins/nbins
       Ra=c(0, 1)
     }
    if(which%in% list.of.studies[c(11:15, 30:35)]) {
      bins=seq(-5, 5, length=nbins+1)
      Ra=c(-Inf, Inf)
    }
    if(which%in% list.of.studies[c(16:20,36:40)]) {
      bins=seq(0, 3, length=nbins+1)
      Ra=c(0, Inf)
    }
    vals=(bins[-1]+bins[-(nbins+1)])/2 # midpoints of intervals
  }
 
# 1  
  if(which=="uniform") {
    rnull=function(n) runif(n)
    dnull=function(x) dunif(x)
    pnull=function(x) punif(x)
  }
# 2  
  if(which=="linear") {
    rnull=function(n) 10/6*(-0.7+sqrt(0.49+1.2*runif(n)))
    dnull=function(x) 0.6*x+0.7
    pnull=function(x) 0.3*x^2+0.7*x
  }  
# 3
  if(which=="quadratic") {
     rnull=function(n) {
       a=-0.5
       x=rep(0,n)
       for(i in 1:n) {
         repeat {
           z=runif(1)
           if(runif(1)<3*a*(z-1/2)^2+1-a/4) {
             x[i]=z
             break
           }
         } 
       }
       x
     }   
     dnull=function(x) {a=-0.5;3*a*(x-1/2)^2+1-a/4}
     pnull=function(x) {a=-0.5;a*(x-1/2)^3+(1-a/4)*x+a/8}
  }          
#4  
  if(which=="uniform-betabump") {
    rnull=function(n) c(runif(floor(n*0.825)), 
                        rbeta(ceiling(n*0.175), 25, 25))
    dnull=function(x) 0.825*dunif(x)+0.175*dbeta(x, 25, 25)
    pnull=function(x) 0.825*punif(x)+0.175*pbeta(x, 25, 25)
  }     
#5  
  if(which=="sine") {
    rnull=function(n) acos(2*runif(n)-1)/pi   
    dnull=function(x) pi*sin(pi*x)/2
    pnull=function(x) (1-cos(pi*x))/2
  }    
#6   
  if(which=="betaaa") {
    rnull=function(n) rbeta(n, 3, 3)
    dnull=function(x) dbeta(x, 3, 3)
    pnull=function(x) pbeta(x, 3, 3)
  }
#7   
  if(which=="beta2a") {
    rnull=function(n) rbeta(n, 2, 2.5)
    dnull=function(x) dbeta(x, 2, 2.5)
    pnull=function(x) pbeta(x, 2, 2.5)
  }
#8
  if(which=="triangular") {
    rnull=function(n) {
      y=runif(0.35*n)
      y=ifelse(y<1/2, sqrt(y/2), 1-sqrt((1-y)/2))
      c(y, runif(0.65*n))
   }    
   dnull=function(x) 0.65*dunif(x)+0.35*ifelse(x<0.5, 4*x, 4-4*x)
   pnull=function(x) 0.65*punif(x)+0.35*ifelse(x<0.5, 2*x^2, 1-2*(1-x)^2)
  }    
#9
  if(which=="uniformmixture") {
    rnull=function(n) c(runif(n/6, 0, 1/3), runif(n/3,1/3,2/3), runif(n/2, 2/3, 1))
    dnull=function(x) (dunif(x, 0, 1/3)+2*dunif(x,1/3,2/3)+3*dunif(x, 2/3, 1))/6
    pnull=function(x) (punif(x, 0, 1/3)+2*punif(x,1/3,2/3)+3*punif(x, 2/3, 1))/6
  }     
#10
  if(which=="betamixture") {
    rnull=function(n) c(runif(n/2), rbeta(n/2, 2, 2))
    dnull=function(x) (dunif(x)+dbeta(x, 2, 2))/2
    pnull=function(x) (punif(x)+pbeta(x, 2, 2))/2
  }
#11  
  if(which=="normal-pure") {
    rnull=function(n) rnorm(n)
    dnull=function(x) dnorm(x)
    pnull=function(x) pnorm(x)
  }
#12  
  if(which=="normal-mix1") {
    rnull=function(n) c(rnorm(floor(n/2)), rnorm(ceiling(n/2), 0, 3))
    dnull=function(x) (dnorm(x)+dnorm(x, 0, 3))/2
    pnull=function(x) (pnorm(x)+pnorm(x, 0, 3))/2
  }
#13  
  if(which=="normal-mix2") {
    rnull=function(n) c(rnorm(floor(n/2)), rnorm(ceiling(n/2), 1))
    dnull=function(x) (dnorm(x)+dnorm(x, 1))/2
    pnull=function(x) (pnorm(x)+pnorm(x, 1))/2
  }
#14
  if(which=="normal-mix3") {
    rnull=function(n) c(rnorm(n/4, -3), rnorm(n/2), rnorm(n/4,5))
    dnull=function(x) dnorm(x, -3)/4+dnorm(x)/2+dnorm(x, 5)/4
    pnull=function(x) pnorm(x, -3)/4+pnorm(x)/2+pnorm(x, 5)/4
  }
#15  
  if(which=="normal-mix4") {
    rnull=function(n) c(rnorm(n/2), rnorm(n/4, 0, 3), rnorm(n/4, 0, 5))
    dnull=function(x) (2*dnorm(x)+dnorm(x, 0, 3)+dnorm(x, 0, 5))/4
    pnull=function(x) (2*pnorm(x)+pnorm(x, 0, 3)+pnorm(x, 0, 5))/4
  }
#16
  if(which=="gamma1") {
    rnull=function(n) c(rexp(0.2*n, 1), rgamma(0.8*n, 2, 0.5))
    dnull=function(x) 0.2*dexp(x, 1)+0.8*dgamma(x, 2, 0.5)
    pnull=function(x) 0.2*pexp(x, 1)+0.8*pgamma(x, 2, 0.5)
  }
#17
  if(which=="gamma2") {
    rnull=function(n) c(rgamma(0.2*n, 2, 2), rgamma(0.8*n, 3, 5))
    dnull=function(x) 0.2*dgamma(x, 2, 2)+0.8*dgamma(x, 3, 5)
    pnull=function(x) 0.2*pgamma(x, 2, 2)+0.8*pgamma(x, 3, 5)
  }
#18
  if(which=="gammamixture") {
    rnull=function(n) c(rgamma(n/2, 2, 1), rgamma(n/2, 10,1))
    dnull=function(x) (dgamma(x, 2, 1)+dgamma(x, 10,1))/2
    pnull=function(x) (pgamma(x, 2, 1)+pgamma(x, 10,1))/2
  }
#19  
  if(which=="exponential-mixture") {
    rnull=function(n) c(rexp(floor(n/2), 1), 
            rexp(ceiling(n/2), 3))
    dnull=function(x) (dexp(x, 1)+dexp(x, 3))/2
    pnull=function(x) (pexp(x, 1)+pexp(x, 3))/2
  }
#20
  if(which=="noncentral") {
    rnull=function(n) rchisq(n, 5, 0.875)
    dnull=function(x) dchisq(x, 5, 0.875)
    pnull=function(x) pchisq(x, 5, 0.875)
  }
#21
  if(which=="uniform-linear") {
    rnull=function(n) runif(n)
    dnull=function(x) dunif(x)
    pnull=function(x) punif(x)
    ralt=function(n) 10/6*(-0.7+sqrt(0.49+1.2*runif(n)))
    palt=function(x) 0.3*x^2+0.7*x
  }
#22
  if(which=="uniform-quadratic") {
    rnull=function(n) runif(n)
    dnull=function(x) dunif(x)
    pnull=function(x) punif(x)
    a=ifelse(data_type=="cont", -1.2, -1)
    ralt=function(n) {
      x=rep(0,n)
      for(i in 1:n) {
        repeat {
          z=runif(1)
          if(runif(1)<3*a*(z-1/2)^2+1-a/4) {
            x[i]=z
            break
          }
        } 
      }
      x
    }   
    palt=function(x) {
      a*(x-1/2)^3+(1-a/4)*x+a/8
    }
  }
#23
  if(which=="uniform-bump") {
    alpha=0.125
    rnull=function(n) runif(n)
    dnull=function(x) dunif(x)
    pnull=function(x) punif(x)
    ralt=function(n) c(runif(floor(n*(1-alpha))), 
        rnorm(ceiling(n*alpha), 0.5, 0.05))
    palt=function(x) (1-alpha)*punif(x)+alpha*pnorm(x, 0.5, 0.05)
  }
#24
  if(which=="uniform-sine") {
    rnull=function(n) runif(n)
    dnull=function(x) dunif(x)
    pnull=function(x) punif(x)
    a=ifelse(data_type=="cont", 0.5, 0.4)
    ralt=function(n) {
      y = NULL
      repeat { 
        z=runif(n)
        I=ifelse(runif(n)<1+a*sin(4*pi*z), TRUE, FALSE)
        y = c(y, z[I]) 
        if(length(y)>n) break
      }
      y
    }
    palt=function(x) x-a*cos(4*pi*x)/4/pi+a/4/pi
  }  
#25
  if(which=="beta22.betaaa") {
    a=ifelse(data_type=="cont", 2.7, 2.85)
    rnull=function(n) rbeta(n, 2, 2)
    dnull=function(x) dbeta(x, 2, 2)
    pnull=function(x) pbeta(x, 2, 2)
    ralt=function(n) rbeta(n, a, a)
    palt=function(x) pbeta(x, a, a)
  }
#26
  if(which=="beta22.beta2a") {
    a=2.35
    rnull=function(n) rbeta(n, 2, 2)
    dnull=function(x) dbeta(x, 2, 2)
    pnull=function(x) pbeta(x, 2, 2)
    ralt=function(n) rbeta(n, 2, a)
    palt=function(x) pbeta(x, 2, a)
  }  
#27
  if(which=="uniform.uniformmixture") {
    a=ifelse(data_type=="cont", 0.6, 0.55)
    rnull=function(n) runif(n)
    dnull=function(x) dunif(x)
    pnull=function(x) punif(x)
    ralt=function(n) {
      z = runif(n)
      y = ifelse(z<1-a, z/2/(1-a), (z-1+a)/2/a+0.5)
    }
    palt=function(x) {
      ifelse(x<1-a, (1-a)*punif(x, 0, 1/2),
             1-a+a*punif(x, 1/2, 1))
    }
  } 
#28
  if(which=="uniform.betamixture") {
    a=0.5
    rnull=function(n) runif(n)
    dnull=function(x) dunif(x)
    pnull=function(x) punif(x)
    ralt=function(n) c(runif((1-a)*n), rbeta(a*n, 2, 2))
    palt=function(x) (1-a)*punif(x)+a*pbeta(x, 2, 2)
  }  
#29
  if(which=="uniform.triangular") {
    a=0.44
    rnull=function(n) runif(n)
    dnull=function(x) dunif(x)
    pnull=function(x) punif(x)
    ralt=function(n) {
      y=runif(a*n)
      y=ifelse(y<1/2, sqrt(y/2), 1-sqrt((1-y)/2))
      c(y, runif((1-a)*n))
    }    
    palt=function(x) (1-a)*punif(x)+a*ifelse(x<0.5, 2*x^2, 1-2*(1-x)^2)
  }     
#30
  if(which=="double.exponential") {
    a=ifelse(data_type=="cont", 1.2, 1.3)
    rnull=function(n) ifelse(runif(n)<0,-rexp(n, 1), rexp(n, 1))
    dnull=function(x) dexp(abs(x), 1)/2
    pnull=function(x) ifelse(x<0,(1-pexp(-x, 1))/2,1/2+pexp(x,1)/2)
    ralt=function(n) ifelse(runif(n)<0,-rexp(n, a), rexp(n, a))
    palt=function(x) ifelse(x<0,(1-pexp(-x, a))/2,1/2+pexp(x,a)/2)
  }
#31
  if(which=="normal.onestretch") {
    a=ifelse(data_type=="cont", 1.25, 1.35)
    rnull=function(n) rnorm(n)
    dnull=function(x) dnorm(x)
    pnull=function(x) pnorm(x)
    ralt=function(n)  c(-abs(rnorm(n/2)), abs(rnorm(n/2, 0,a)))
    palt=function(x) ifelse(x<0,pnorm(x), pnorm(x,0,a))
  }
#32
  if(which=="normal.t") {
    a=ifelse(data_type=="cont", 6, 5)
    rnull=function(n) rnorm(n)
    dnull=function(x) dnorm(x)
    pnull=function(x) pnorm(x)
    ralt=function(n) rt(n, a)
    palt=function(x) pt(x, a)
  }  
#33
  if(which=="normal.outlier1") {
    a=ifelse(data_type=="cont", 0.015, 0.085)
    rnull=function(n) rnorm(n)
    dnull=function(x) dnorm(x)
    pnull=function(x) pnorm(x)
    ralt=function(n) c(rnorm(floor((1-a)*n)), rnorm(ceiling(a*n), 5))
    palt=function(x) (1-a)*pnorm(x)+a*pnorm(x, 5)
  }    
#34
  if(which=="normal.outlier2") {
    a=ifelse(data_type=="cont", 0.012, 0.03)
    rnull=function(n) rnorm(n)
    dnull=function(x) dnorm(x)
    pnull=function(x) pnorm(x)
    ralt=function(n) c(rnorm((1-2*a)*n), runif(a*n, -5, -3), runif(a*n, 3, 5))
    palt=function(x) (1-2*a)*pnorm(x)+a*punif(x, -5, -3)+a*punif(x, 3, 5)
  }    
#35
  if(which=="normal.normalmixture") {
    a=ifelse(data_type=="cont", 0.62, 0.68)
    rnull=function(n) rnorm(n)
    dnull=function(x) dnorm(x)
    pnull=function(x) pnorm(x)
    ralt=function(n) c(rnorm(n/2, -a), rnorm(n/2, a))
    palt=function(x) (pnorm(x,-a)+pnorm(x, a))
  }    
#36
  if(which=="exponential.gamma") {
    a=ifelse(data_type=="cont", 1.22, 1.28)
    rnull=function(n) rexp(n, 1)
    dnull=function(x) dexp(x, 1)
    pnull=function(x) pexp(x, 1)
    ralt=function(n) rgamma(n, 1, a)
    palt=function(x) pgamma(x, 1, a)
  }
#37
  if(which=="exponential.weibull") {
    a=ifelse(data_type=="cont", 1.22, 1.3)
    rnull=function(n) rexp(n, 1)
    dnull=function(x) dexp(x, 1)
    pnull=function(x) pexp(x, 1)
    ralt=function(n) rweibull(n, 1, a)
    palt=function(x) pweibull(x, 1, a)
  }
#38
  if(which=="exponential.bump") {
    a=0.124
    rnull=function(n) rexp(n, 1)
    dnull=function(x) dexp(x, 1)
    pnull=function(x) pexp(x, 1)
    ralt=function(n) c(rexp((1-a)*n, 1), rnorm(a*n, 0.5, 0.05))
    palt=function(x) (1-a)*pexp(x, 1)+a*pnorm(x, 0.5, 0.05)
  }  
#39
  if(which=="chisquare.noncentral") {
    a=ifelse(data_type=="cont", 0.6, 3)
    rnull=function(n) rchisq(n, 5)
    dnull=function(x) dchisq(x, 5)
    pnull=function(x) pchisq(x, 5)
    ralt=function(n) rchisq(n, 5, a)
    palt=function(x) pchisq(x, 5, a)
  } 
#40
  if(which=="gamma.gammamixture") {
    a=0.275
    rnull=function(n) rgamma(n, 2, 3)
    dnull=function(x) dgamma(x, 2, 3)
    pnull=function(x) pgamma(x, 2, 3)
    ralt=function(n) c(rgamma((1-a)*n, 2, 3), rgamma(a*n, 3, 3))
    palt=function(x) (1-a)*pgamma(x, 2, 3)+a*pgamma(x, 3, 3)
  }
# End of definitions  
#   

# Create output  
  if(data_type=="cont") {
    f=function() {list(x=rnull(n), y=ralt(n))} 
  }
  else {
    prob1=diff(pnull(bins))
    prob2=diff(palt(bins))
    f=function() {
      x=c(rmultinom(1, n, prob1))
      y=c(rmultinom(1, n, prob2)) 
      vals=vals
      I=seq_along(x)[x+y>0]
      list(
        x=x[I],
        y=y[I], 
        vals=vals[I]
      )
    }
  }  
  return(list(f=f, rnull=rnull, dnull=dnull, pnull=pnull,
              ralt=ralt, palt=palt))   
}
