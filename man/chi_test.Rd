% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chi_test.R
\name{chi_test}
\alias{chi_test}
\title{This function runs the chi-square test for continuous or discrete data}
\usage{
chi_test(dta, nbins = c(50, 10), minexpcount = 5, typeTS, ponly = FALSE)
}
\arguments{
\item{dta}{a list with two elements for continuous data or three elements for discrete data, Can also include weights for continuous data}

\item{nbins}{=c(50, 10) number of desired bins}

\item{minexpcount}{=5 smallest number of counts required in each bin}

\item{typeTS}{=5  type of problem, continuous/discrete, with/without weights}

\item{ponly}{Should the p value alone be returned?}
}
\value{
A list with the test statistics, the p value and the degree of freedom for each test
}
\description{
This function runs the chi-square test for continuous or discrete data
}
\keyword{internal}
