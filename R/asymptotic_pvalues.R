#' This function finds the p values of several tests based on large sample theory
#' @param  x a vector of test statistics
#' @param  n size of sample 1
#' @param  m size of sample 2
#' @return A vector of p values. 

asymptotic_pvalues = function(x, n, m) {
   KS=function(x, n, m) {
    x=x*sqrt(m)*sqrt(n)/sqrt(m+n)
    k=1:10
    1-2*sum((-1)^(k-1)*exp(-2*k^2*x^2))
   }
  Kuiper = function(kp, nx, ny) {
    n <- nx * (ny/(nx + ny))
    if (kp < 0 | kp > 2) {
        stop("The test statistic much be in the range of (0,2) by definition of the Kuiper Test")
    }
    else if (kp < 2/n) {
        PVAL <- 1 - factorial(n) * (kp - 1/n)^(n - 1)
    }
    else if (kp < 3/n) {
        k <- -(n * kp - 1)/2
        r <- sqrt(k^2 - (n * kp - 2)/2)
        a <- -k + r
        b <- -k - r
        PVAL <- 1 - factorial(n - 1) * (b^(n - 1) * (1 - a) - 
            a^(n - 1) * (1 - b))/(n^(n - 2) * (b - a))
    }
    else if ((kp > 0.5 && n%%2 == 0) | (kp > (n - 1)/(2 * n) && 
        n%%2 == 1)) {
        temp <- as.integer(floor(n * (1 - kp))) + 1
        c <- matrix(0, ncol = 1, nrow = temp)
        for (t in 0:temp) {
            y <- kp + t/n
            tt <- y^(t - 3) * (y^3 * n - y^2 * t * (3 - 2/n)/n - 
                t * (t - 1) * (t - 2)/n^2)
            p_temp <- choose(n, t) * (1 - kp - t/n)^(n - t - 
                1) * tt
            c[t] <- p_temp
        }
        PVAL <- sum(c)
    }
    else {
        z <- kp * sqrt(n)
        s1 <- 0
        term <- 1e-12
        abs <- 1e-100
        for (m in 1:1e+08) {
            t1 <- 2 * (4 * m^2 * z^2 - 1) * exp(-2 * m^2 * z^2)
            s <- s1
            s1 <- s1 + t1
            if ((abs(s1 - s)/(abs(s1) + abs(s)) < term) | (abs(s1 - 
                s) < abs)) 
                break
        }
        s2 <- 0
        for (m in 1:1e+08) {
            t2 <- m^2 * (4 * m^2 * z^2 - 3) * exp(-2 * m^2 * 
                z^2)
            s <- s2
            s2 <- s2 + t2
            if ((abs(s2 - s)/(abs(s2) + abs(s))) < term | (abs(s1 - 
                s) < abs)) 
                break
        }    
        PVAL <- s1 - 8 * kp/(3 * sqrt(n)) * s2
    }
    PVAL
  }
  CvMp = function(x) {
    interp=function(z, x, a) {
      slope = diff(a)/diff(x)
      a[1] + slope*(z-x[1])
    }
    z = x 
    if(z<0.3473) out=interp(z, c(0,0.3473), c(1, 0.1))
    if(z>1.16786) out=interp(z, c(0.001,0), c(1.16786, 2))
    if(z>0.3473 & z<0.46136) 
       out=interp(z, c(0.3473, 0.46136), c(0.1, 0.05)) 
    if(z>0.46136 & z<0.74346) 
       out=interp(z, c(0.46136, 0.74346), c(0.05, 0.01)) 
    if(z>0.74346 & z<1.16786) 
       out=interp(z, c(0.74346, 1.16786), c(0.01, 0.001)) 
    out
  }
  ADp <- function(x) {
    if(x < 2)
      x=exp(-1.2337141/x)/sqrt(x)*(2.00012+(.247105-
              (.0649821-(.0347962-(.011672-.00168691*x)*x)*x)*x)*x)
    else
      x=exp(-exp(1.0776-(2.30695-(.43424-
                (.082433-(.008056 -.0003146*x)*x)*x)*x)*x))
    1-x
    
  }
  out=c(1-KS(x[1], n, m), Kuiper(x[2], n, m), CvMp(x[3]), ADp(x[4]))
  names(out) = c("KS", "Kuiper", "CvM", "AD")
  out
}
