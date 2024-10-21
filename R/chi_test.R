#' This function runs the chi-square test for continuous or discrete data
#' @param  dta  a list with two elements for continuous data or three elements for discrete data, Can also include weights for continuous data
#' @param  nbins =c(50, 10) number of desired bins
#' @param  minexpcount =5 smallest number of counts required in each bin
#' @param  typeTS =5  type of problem, continuous/discrete, with/without weights
#' @param  ponly Should the p value alone be returned?
#' @keywords internal
#' @return A list with the test statistics, the p value and the degree of freedom for each test

chi_test = function(dta, nbins=c(50, 10), minexpcount=5, typeTS, ponly=FALSE) {
  combine.bins = function(x, y, minexpcount=5, maxbins=1000, wx, wy) {
    Weights = ifelse(missing(wx), FALSE, TRUE) 
    nx = sum(x)
    ny = sum(y)
    n = nx+ny
    k = length(x)
    while ( (min(x+y)<minexpcount) || (k>maxbins) ){
      i = which.min(x+y)
      if(i==1) {
        x = c(x[1]+x[2], x[3:k])
        y = c(y[1]+y[2], y[3:k])
        if(Weights) {
          wx = c(wx[1]+wx[2], wx[3:k])
          wy = c(wy[1]+wy[2], wy[3:k])        }
      }
      if(i==2) {
         if((x+y)[1]<(x+y)[3]) {
           x = c(x[1]+x[2], x[3:k])
           y = c(y[1]+y[2], y[3:k])
           if(Weights) {
             wx = c(wx[1]+wx[2], wx[3:k])
             wy = c(wy[1]+wy[2], wy[3:k])
           }
         }
         else {
           x = c(x[1], x[2]+x[3], x[4:k])
           y = c(y[1], y[2]+y[3], y[4:k])
           if(Weights) {
             wx = c(wx[1], wx[2]+wx[3], wx[4:k])
             wy = c(wy[1], wy[2]+wy[3], wy[4:k])
           }
         }
      }
      if(i==k-1) {
         if((x+y)[k-2]<(x+y)[k]) {
            x = c(x[1:(k-3)], x[k-2]+x[k-1], x[k])
            y = c(y[1:(k-3)], y[k-2]+y[k-1], y[k])
            if(Weights) {
              wx = c(wx[1:(k-3)], wx[k-2]+wx[k-1], wx[k])
              wy = c(wy[1:(k-3)], wy[k-2]+wy[k-1], wy[k])
            }
         }
         else {
            x = c(x[1:(k-2)], x[k-1]+x[k])
            y = c(y[1:(k-2)], y[k-1]+y[k])
            if(Weights) {
              wx = c(wx[1:(k-2)], wx[k-1]+wx[k])
              wy = c(wy[1:(k-2)], wy[k-1]+wy[k])
            }
         }
      }
      if(i==k) {
         x = c(x[1:(k-2)], x[k-1]+x[k])
         y = c(y[1:(k-2)], y[k-1]+y[k])
         if(Weights) {
           wx = c(wx[1:(k-2)], wx[k-1]+wx[k])
           wy = c(wy[1:(k-2)], wy[k-1]+wy[k])
         }
      }
      if(i>2 && i<k-1) {
        if( (x+y)[i-1]<(x+y)[i+1] ) {
          x = c(x[1:(i-2)], x[i-1]+x[i], x[(i+1):k])
          y = c(y[1:(i-2)], y[i-1]+y[i], y[(i+1):k])
          if(Weights) {
            wx = c(wx[1:(i-2)], wx[i-1]+wx[i], wx[(i+1):k])
            wy = c(wy[1:(i-2)], wy[i-1]+wy[i], wy[(i+1):k])
          }
        }
        else {
          x = c(x[1:(i-1)], x[i]+x[i+1], x[(i+2):k])
          y = c(y[1:(i-1)], y[i]+y[i+1], y[(i+2):k])
          if(Weights) {
            wx = c(wx[1:(i-1)], wx[i]+wx[i+1], wx[(i+2):k])
            wy = c(wy[1:(i-1)], wy[i]+wy[i+1], wy[(i+2):k])  
          }
        }
      }
      k = k-1
    }
    if(!Weights) out = list(x=x, y=y)
    else out = list(wx=cbind(x, wx), wy=cbind(y,wy))
    out
  }  
  dotest = function(dta) {
    x = dta$x
    y = dta$y
    tmp=sqrt(sum(x)/sum(y))
    chi = sum((x/tmp-y*tmp)^2/(x+y))
    p = 1-stats::pchisq(chi, length(x)-1)
    c(chi, p, length(x)-1)
  }
  dotestweightedcontinuous = function(x, y, wx, wy) {
    n = wx[, 1]
    v = wx[, 2]
    o = wy[, 1]
    u = wy[, 2]
    t = length(x)/length(y)
    chi = sum( (n-t*o)^2/(v+t^2*u) )
    p = 1-stats::pchisq(chi, length(n)-1)
    c(chi, p, length(n)-1)
  }
  dotestweighteddiscrete = function(dta) {
    n = dta[[1]][, 1]
    v = dta[[1]][, 2]
    o = dta[[2]][, 1]
    u = dta[[2]][, 2]
    t = sum(n)/sum(o)
    chi = sum( (n-t*o)^2/(v+t^2*u) )
    p = 1-stats::pchisq(chi, length(n)-1)
    c(chi, p, length(n)-1)
  }
  if(typeTS %in% c(1,2)) {
    if(typeTS==2) { # continuous data, with weights
      dta$wx = dta$wx[order(dta$x)]
      dta$x = sort(dta$x)
      dta$wy = dta$wy[order(dta$y)]
      dta$y = sort(dta$y)
    }
    xy = c(dta$x, dta$y)
    out = matrix(0, 4, 3)
    colnames(out) = c("chi", "p.value", "df")
    rownames(out) = c("ES large", "ES small", "EP large", "EP small")
    # Equal size bins 
    for(j in 1:2) {
      bins = seq(min(xy),max(xy), length=nbins[j])
      if(typeTS==1) {
        x = bincounter(dta$x, bins)
        y = bincounter(dta$y, bins)
        dta1 = combine.bins(x, y, minexpcount, nbins[j])
        out[j, ] = dotest(dta1)        
      }
      else {
        wx = wbincounter(dta$x, bins, dta$wx)
        wy = wbincounter(dta$y, bins, dta$wy)
        dta1 = combine.bins(wx[,1], wy[,1],
                    minexpcount, nbins[j], wx[,2], wy[,2])
        out[j, ] = dotestweightedcontinuous(dta$x, dta$y, dta1$wx, dta1$wy)
      }

    }
    # Equal probability bins
    for(j in 1:2) {
      bins = quantile(xy, (0:nbins[j])/nbins[j])
      if(typeTS==1) {
        x = bincounter(dta$x, bins)
        y = bincounter(dta$y, bins)
        dta1 = combine.bins(x, y, minexpcount, nbins[j])
        out[j+2, ] = dotest(dta1)        
      }
      else {
        wx = wbincounter(dta$x, bins, dta$wx)
        wy = wbincounter(dta$y, bins, dta$wy)
        dta1 = combine.bins(wx[,1], wy[,1],
                            minexpcount, nbins[j], wx[,2], wy[,2])
        out[j+2, ] = dotestweightedcontinuous(dta$x, dta$y, dta1$wx, dta1$wy)
      }
    }
  }
  else { # discrete data
    out = matrix(0,2,3)
    colnames(out) = c("chi", "p.value", "df")
    rownames(out) = c("large", "small")
    if(typeTS==6) { #with weights
      dta1 = combine.bins(dta$x, dta$y, 
         minexpcount, nbins[1], dta$TSextra[[1]], dta$TSextra[[2]])          
      out[1, ] = dotestweighteddiscrete(dta1)      
      dta1 = combine.bins(dta$x, dta$y, 
         minexpcount, nbins[2], dta$TSextra[[1]], dta$TSextra[[2]])        
      out[2, ] = dotestweighteddiscrete(dta1)      
    }
    else {
      dta1 = combine.bins(dta$x, dta$y, minexpcount, nbins[1])          
      out[1, ] = dotest(dta1)      
      dta1 = combine.bins(dta$x, dta$y, minexpcount, nbins[2])          
      out[2, ] = dotest(dta1)      
    }
    
  } 
  if(ponly) return(round(out[, 2], 4))
  list(statistics = round(out[, 1], 2), 
       p.values = round(out[, 2], 4), df=round(out[,3]))
}
