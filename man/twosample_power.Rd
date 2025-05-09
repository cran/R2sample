% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/twosample_power.R
\name{twosample_power}
\alias{twosample_power}
\title{Find the power of various two sample tests using Rcpp and parallel computing.}
\usage{
twosample_power(
  f,
  ...,
  TS,
  TSextra,
  alpha = 0.05,
  B = 1000,
  nbins = c(50, 10),
  minexpcount = 5,
  UseLargeSample,
  samplingmethod = "independence",
  rnull,
  SuppressMessages = FALSE,
  maxProcessor
)
}
\arguments{
\item{f}{function to generate a list with data sets x, y and (optional) vals, weights}

\item{...}{additional arguments passed to f, up to 2}

\item{TS}{routine to calculate test statistics for non-chi-square tests}

\item{TSextra}{additional info passed to TS, if necessary}

\item{alpha}{=0.05, the level of the hypothesis test}

\item{B}{=1000, number of simulation runs.}

\item{nbins}{=c(50,10), number of bins for chi large and chi small.}

\item{minexpcount}{=5 minimum required count for chi square tests}

\item{UseLargeSample}{should p values be found via large sample theory if n,m>10000?}

\item{samplingmethod}{=independence or MCMC in discrete data case}

\item{rnull}{a function that generates data from a model, possibly with parameter estimation.}

\item{SuppressMessages}{= FALSE print informative messages?}

\item{maxProcessor}{maximum number of cores to use. If maxProcessor=1 no parallel computing is used.}
}
\value{
A numeric vector of power values.
}
\description{
Find the power of various two sample tests using Rcpp and parallel computing.
}
\examples{
 f=function(mu) list(x=rnorm(25), y=rnorm(25, mu))
 twosample_power(f, mu=c(0,2), B=100, maxProcessor = 1)
 f=function(n, p) list(x=table(sample(1:5, size=1000, replace=TRUE)), 
       y=table(sample(1:5, size=n, replace=TRUE, 
       prob=c(1, 1, 1, 1, p))), vals=1:5)
 twosample_power(f, n=c(1000, 2000), p=c(1, 1.5), B=100, maxProcessor = 1)
}
