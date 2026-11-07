#' Settings of the HHG test
#'
#' The HHG test (Heller, Heller, Kaufman, Brill and Gorfine 2016) is run in its
#' equipartition version: the pooled sample is split into \code{natoms} groups of
#' consecutive order statistics (continuous data) or into its distinct values
#' (discrete data), the likelihood ratio statistics of all partitions of these
#' atoms into m contiguous cells are averaged, for m=2,..,\code{mmax}, and the
#' results are combined with the MinP rule. For this the null distribution of each
#' statistic is found once, from \code{Bnull} random relabelings of the pooled sample.
#' @keywords internal
hhg_settings <- list(natoms=40L, mmax=10L, Bnull=2000L, seed=20161L)

# null distributions for continuous data depend only on the sample sizes, so they
# are found once per session for each pair of sample sizes
.hhg_cache <- new.env(parent=emptyenv())

# Run expr with its own random number stream, so that the HHG null distribution
# does not change the random numbers used by anything else (and is the same
# whenever the same data or sample sizes are used).
hhg_with_seed <- function(expr) {
  had <- exists(".Random.seed", envir=globalenv(), inherits=FALSE)
  if(had) old <- get(".Random.seed", envir=globalenv(), inherits=FALSE)
  on.exit({
    if(had) assign(".Random.seed", old, envir=globalenv())
    else if(exists(".Random.seed", envir=globalenv(), inherits=FALSE))
      rm(".Random.seed", envir=globalenv())
  })
  set.seed(hhg_settings$seed)
  expr
}

#' Null distribution of the HHG statistics for continuous data
#'
#' The test is distribution-free for continuous data, so the null distribution
#' depends only on the sample sizes. It is found once per session for each pair
#' of sample sizes.
#' @param nx size of the first sample
#' @param ny size of the second sample
#' @keywords internal
#' @return a list with the number of atoms A, mmax and the null distribution, or
#'   NULL if the samples are too small
#' @noRd
hhg_null_cont <- function(nx, ny) {
  N <- nx+ny
  A <- min(N, hhg_settings$natoms)
  mmax <- min(hhg_settings$mmax, A)
  if(mmax<2) return(NULL)
  key <- paste(nx, ny, A, mmax)
  if(is.null(.hhg_cache[[key]]))
    .hhg_cache[[key]] <- hhg_with_seed(
      HHG_null(HHG_atoms(N, A), nx, mmax, hhg_settings$Bnull))
  list(A=A, mmax=mmax, null=.hhg_cache[[key]])
}

#' p value of the HHG test for continuous data without permutation
#'
#' Because the test is distribution-free, the null distribution of the HHG
#' statistic is found by simulation once for given sample sizes, and the p value
#' is the proportion of null values at least as large as the observed one. This is
#' used for large samples.
#' @param x first data set
#' @param y second data set
#' @keywords internal
#' @return p value
#' @noRd
hhg_pvalue_cont <- function(x, y) {
  nul <- hhg_null_cont(length(x), length(y))
  if(is.null(nul)) return(NA)
  s <- HHG_cont(x, y, nul$A, nul$mmax, nul$null)
  ms <- attr(nul$null, "minp")
  (1+sum(ms >= s-1e-9*abs(s)))/(length(ms)+1)
}

#' Add the HHG test to the built-in test statistics
#'
#' @param TS the routine for the built-in test statistics (TS_cont or TS_disc)
#' @param x first data set (continuous) or its counts (discrete)
#' @param y second data set (continuous) or its counts (discrete)
#' @param Continuous is the data continuous?
#' @keywords internal
#' @return a function with the same arguments as TS that also returns HHG
#' @noRd
add_HHG <- function(TS, x, y, Continuous) {
  force(TS)  # the caller replaces its TS by the result, so TS has to be evaluated now
  mmax <- hhg_settings$mmax
  if(Continuous) {
    nul <- hhg_null_cont(length(x), length(y))
    if(is.null(nul)) return(TS)
    A <- nul$A; mmax <- nul$mmax; nullHHG <- nul$null
    function(x, y) c(TS(x, y), HHG=HHG_cont(x, y, A, mmax, nullHHG))
  } else {
    if(length(x)<2) return(TS)
    mmax <- min(mmax, length(x))
    # null distribution for the atoms of this data set (its distinct values)
    nullHHG <- hhg_with_seed(
      HHG_null(as.integer(round(x+y)), as.integer(round(sum(x))), mmax,
               hhg_settings$Bnull))
    function(x, y, vals, adw) c(TS(x, y, vals, adw), HHG=HHG_disc(x, y, mmax, nullHHG))
  }
}
