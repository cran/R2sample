## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(R2sample)

## -----------------------------------------------------------------------------
set.seed(123)

## -----------------------------------------------------------------------------
x1 = rnorm(100)
y1 = rnorm(120)
twosample_test(x1, y1, B=500)

## -----------------------------------------------------------------------------
twosample_test_adjusted_pvalue(x1, y1)

## -----------------------------------------------------------------------------
myTS1 = function(x, y) {
   out = c(0, 0)
   out[1] = abs(sd(x) - sd(y))
   out[2] = abs(mean(x)/sd(x) - mean(y)/sd(y))
   names(out) = c("std", "std t test")
   out
}

## -----------------------------------------------------------------------------
twosample_test(x1, y1, TS=myTS1, B=500)

## -----------------------------------------------------------------------------
x2 = table(c(0:5,rbinom(1000, 5, 0.5)))-1
y2 = table(c(0:5,rbinom(1200, 5, 0.55)))-1
rbind(x2, y2)
twosample_test(x2, y2, vals=0:5, B=500)$p.values

## -----------------------------------------------------------------------------
twosample_test(x2, y2, vals=0:5, TS=R2sample::myTS2, B=500)

## -----------------------------------------------------------------------------
rnull=function(dta) {
   list(x=dta$x,
        y=rnorm(length(dta$y), mean(dta$x), sd(dta$x)))
}
twosample_test(x1, y1, rnull=rnull, B=500)

## -----------------------------------------------------------------------------
x=rnorm(10)
y=rnorm(12)
twosample_test(x, y, B=500, doMethods=c("KS","AD"))

## ----eval=FALSE---------------------------------------------------------------
# run_shiny()

## ----eval=FALSE---------------------------------------------------------------
# plot_power(twosample_power(
#   f=function(df) list(x=rnorm(100), y=rt(200, df))
# )

## ----eval=FALSE---------------------------------------------------------------
# plot_power(twosample_power(
#   function(p) {
#     vals=0:10
#     x=table(c(vals, rbinom(100, 10, 0.5))) - 1 #Make sure each vector has length 11
#     y=table(c(vals, rbinom(120, 10, p))) - 1  #and names 1-10
#     vals=vals[x+y>0] # only vals with at least one count
#     list(x=x[as.character(vals)],y=y[as.character(vals)], vals=vals)
#   },
#   p=seq(0.5, 0.6, length=5)
# ), "p")

## ----eval=FALSE---------------------------------------------------------------
# plot_power(twosample_power(
#   f=function(n) list(x=rnorm(n), y=rnorm(n, 1)),
#   n=10*1:10), "n")

## ----eval=FALSE---------------------------------------------------------------
# pvals=matrix(0,1000,13)
# for(i in 1:1000)
#   pvals[i, ]=R2sample::twosample_test(runif(100), runif(100), B=1000)$p.values

## ----eval=FALSE---------------------------------------------------------------
# colnames(pvals)=names(R2sample::twosample_test(runif(100), runif(100), B=1000)$p.values)
# p1=apply(pvals[, c("ZK", "ZC", "Wassp1", "Kuiper", "ES small" )], 1, min)
# p2=apply(pvals[, c("KS", "Kuiper", "AD", "CvM", "LR")], 1, min)

## -----------------------------------------------------------------------------
tmp=R2sample::pvaluecdf
Tests=factor(c(rep("Identical Tests", nrow(tmp)),
        rep("Correlated Selection", nrow(tmp)),
        rep("Best Selection", nrow(tmp)),
        rep("Independent Tests", nrow(tmp))),
        levels=c("Identical Tests",  "Correlated Selection", 
                 "Best Selection", "Independent Tests"),
        ordered = TRUE)
dta=data.frame(x=c(tmp[,1],tmp[,1],tmp[,1],tmp[,1]),
          y=c(tmp[,1],tmp[,3],tmp[,2],1-(1-tmp[,1])^4),
          Tests=Tests)
ggplot2::ggplot(data=dta, ggplot2::aes(x=x,y=y,col=Tests))+
  ggplot2::geom_line(linewidth=1.2)+
  ggplot2::labs(x="p value", y="CDF")+
  ggplot2::scale_color_manual(values=c("blue","red", "Orange", "green"))

## -----------------------------------------------------------------------------
chitest = function(x, y, TSextra) {
   nbins=TSextra$nbins
   nx=length(x);ny=length(y);n=nx+ny
   xy=c(x,y)
   bins=quantile(xy, (0:nbins)/nbins)
   Ox=hist(x, bins, plot=FALSE)$counts
   Oy=hist(y, bins, plot=FALSE)$counts
   tmp=sqrt(sum(Ox)/sum(Oy))
   chi = sum((Ox/tmp-Oy*tmp)^2/(Ox+Oy))
   pval=1-pchisq(chi, nbins-1)
   out=ifelse(TSextra$statistic,chi,pval)
   names(out)="ChiSquare"
   out
}
TSextra=list(nbins=5,statistic=FALSE)

## -----------------------------------------------------------------------------
pwr=R2sample::run.studies(chitest, "uniform.linear", TSextra=TSextra, With.p.value = TRUE)
R2sample::plot_power(pwr, "slope")

## ----eval=FALSE---------------------------------------------------------------
# R2sample::run.studies(chitest, TSextra=TSextra, With.p.value=TRUE)

## ----eval=FALSE, message=FALSE------------------------------------------------
# R2sample::run.studies(TRUE, # continuous data/model
#         study=c("uniform.linear", "uniform.quadratic"),
#         param_alt=list(c(0.1, 0.2), 1.3),
#         nsample=2000,
#         alpha=0.1)

