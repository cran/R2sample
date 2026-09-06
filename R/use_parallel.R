#' Determine whether parallel processing should be used
#'
#' Determines whether a computation is expected to take long enough to
#' benefit from parallel processing. If the estimated computation time is
#' below the specified threshold, a single processor is returned.
#' Otherwise, the requested number of processors is returned.
#'
#' @param tm Estimated time, in seconds, for one computation.
#' @param nparam Number of parameter values to be processed.
#' @param B Number of simulation runs.
#' @param maxProcessor Maximum number of processors requested.
#' @param threshold Minimum estimated computation time, in seconds, for
#'   parallel processing to be used. The default is 20.
#'
#' @return An integer giving the number of processors to use.
#'
#' @keywords internal
#'
use_parallel <- function(tm, nparam, B, maxProcessor, threshold = 20) {
  
  # Always respect an explicit request for one processor
  if(maxProcessor <= 1)
    return(1L)
  
  # For short computations, parallel overhead is not worthwhile
  if(tm * nparam * 2 * B < threshold)
    return(1L)
  
  as.integer(maxProcessor)
}