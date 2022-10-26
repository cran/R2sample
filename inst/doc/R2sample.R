## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(R2sample)

## -----------------------------------------------------------------------------
set.seed(123)

## -----------------------------------------------------------------------------
x=rnorm(10)
y=rnorm(12)
twosample_test(x, y, B=500, maxProcessor = 2, doMethod="all")

## -----------------------------------------------------------------------------
x1=table(rbinom(1000, 5, 0.5))
y1=table(rbinom(1200, 5, 0.55))
rbind(0:5,x1, y1)
twosample_test(x1, y1, vals=0:5, B=500, 
               maxProcessor = 2, doMethod="all")$p.values

## -----------------------------------------------------------------------------
set.seed(111)

## -----------------------------------------------------------------------------
set.seed(111)
x=rnorm(10)
y=rnorm(12)
twosample_test(x, y, B=500, maxProcessor = 2, doMethod=c("KS","AD"))

## -----------------------------------------------------------------------------
doMethod = c("chi small", "ZK", "ZC", "Wassp1")

## -----------------------------------------------------------------------------
doMethod = c("chi small", "Kuiper", "ZC", "Wassp1")

## -----------------------------------------------------------------------------
domethod="all"

## ----eval=FALSE---------------------------------------------------------------
#  run_shiny()

## ----eval=FALSE---------------------------------------------------------------
#  plot_power(twosample_power(
#    f=function(mu) list(x=rnorm(100), y=rt(200, mu)),
#    mu=1:10
#  ))

## ----eval=FALSE---------------------------------------------------------------
#  plot_power(twosample_power(
#    function(p) {
#      vals=0:10
#      x=table(c(vals, rbinom(100, 10, 0.5))) - 1 #Make sure each vector has length 11
#      y=table(c(vals, rbinom(120, 10, p))) - 1  #and names 1-10
#      vals=vals[x+y>0] # only vals with at least one count
#      list(x=x[as.character(vals)],y=y[as.character(vals)], vals=vals)
#    },
#    p=seq(0.5, 0.6, length=5)
#  ), "p")

## ----eval=FALSE---------------------------------------------------------------
#  plot_power(twosample_power(
#    f=function(n) list(x=rnorm(n), y=rnorm(n, 1)),
#    n=10*1:10), "n")

