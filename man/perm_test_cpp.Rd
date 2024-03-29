% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{perm_test_cpp}
\alias{perm_test_cpp}
\title{run permutation test.}
\usage{
perm_test_cpp(
  x,
  y,
  vals = 0,
  nbins = as.integer(c(100, 10)),
  B = 5000L,
  doMethod = as.character(c("chi large", "chi small", "t test", "KS", "Kuiper", "CvM",
    "AD", "LR", "ZA", "ZK", "ZC", "Wassp1"))
)
}
\arguments{
\item{x}{A numeric vector.}

\item{y}{A numeric vector.}

\item{vals}{A numeric vector. Indicates discrete data.}

\item{nbins}{Two bin numbers for chi square  test.}

\item{B}{Number of simulation runs.}

\item{doMethod}{A character vector of methods to include.}
}
\value{
A list with test statistics and p values
}
\description{
run permutation test.
}
\keyword{internal}
