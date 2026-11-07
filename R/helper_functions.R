#' makeTSextra
#'
#' Create the TSextra list required by the two-sample routines.
#'
#' @param TSextra a user-supplied list
#' @param samplingmethod sampling scheme for discrete data
#' @param rnull optional function that generates data under a fitted null model
#' @return a list
#' @keywords internal
#' @noRd
makeTSextra <- function(TSextra, samplingmethod, rnull) {
  if(missing(TSextra)) TSextra <- list()
  TSextra$samplingmethod <- samplingmethod
  if(!missing(rnull)) TSextra$rnull <- rnull
  TSextra
}

#' maketypeTS
#'
#' Find typeTS and TS.
#'
#' @param TS a function
#' @param Continuous is data continuous?
#' @param WithWeights are continuous observations weighted?
#' @return a list
#' @keywords internal
#' @noRd
maketypeTS <- function(TS, Continuous, WithWeights) {
  useSingleProcessor <- FALSE
  if(missing(TS)) {
    if(Continuous) {
      if(!WithWeights) {
        typeTS <- 1L
        TS <- TS_cont
      } else {
        typeTS <- 2L
        TS <- TSw_cont
      }
    } else {
      typeTS <- 5L
      TS <- TS_disc
    }
  } else {
    if(any(grepl(".Call", deparse(TS), fixed=TRUE))) {
      message("Parallel Programming is not possible if custom TS is written in C++. Switching to single processor")
      useSingleProcessor <- TRUE
    }
    nargs <- length(formals(TS))
    if(Continuous) {
      if(!(nargs %in% c(2L,3L)))
        stop("TS for continuous data should have either 2 or 3 arguments", call.=FALSE)
      typeTS <- nargs+1L
    } else {
      if(!(nargs %in% c(3L,4L)))
        stop("TS for discrete data should have either 3 or 4 arguments", call.=FALSE)
      typeTS <- nargs+3L
    }
  }
  list(TS=TS, typeTS=typeTS, useSingleProcessor=useSingleProcessor)
}

#' set_seed
#'
#' Validate and set an optional simulation seed.
#'
#' @param seed NULL or a single finite numeric value
#' @return invisibly returns seed
#' @keywords internal
#' @noRd
set_seed <- function(seed=NULL) {
  if(!is.null(seed)) {
    if(!is.numeric(seed) || length(seed)!=1L || is.na(seed) || !is.finite(seed))
      stop("seed must be NULL or a single finite number", call.=FALSE)
    set.seed(seed)
  }
  invisible(seed)
}

#' validate_integer_scalar
#'
#' Validate a positive integer-valued scalar.
#'
#' @param x value to validate
#' @param name argument name used in error messages
#' @param lower smallest permitted value
#' @return invisibly returns the validated integer
#' @keywords internal
#' @noRd
validate_integer_scalar <- function(x, name, lower=1L) {
  if(!is.numeric(x) || length(x)!=1L || is.na(x) || !is.finite(x) ||
     x != floor(x) || x < lower)
    stop(name, " must be a single integer >= ", lower, call.=FALSE)
  invisible(as.integer(x))
}

#' makemaxProcessor
#'
#' Choose the number of processors for a simulation.
#'
#' @param maxProcessor maximum number of processors requested
#' @param dta data set used for a timing check
#' @param TS test-statistic routine
#' @param typeTS test-statistic interface code
#' @param TSextra additional information passed to TS
#' @param B number of simulation runs
#' @param SuppressMessages suppress informative messages?
#' @param useSingleProcessor force a single processor, e.g. for a C++ custom TS
#' @return an integer giving the number of processors to use
#' @keywords internal
#' @noRd
makemaxProcessor <- function(maxProcessor, dta, TS, typeTS, TSextra, B,
                             SuppressMessages=FALSE,
                             useSingleProcessor=FALSE) {
  if(useSingleProcessor) return(1L)

  if(!missing(maxProcessor)) {
    validate_integer_scalar(maxProcessor, "maxProcessor", 1L)
    return(as.integer(maxProcessor))
  }

  m <- parallel::detectCores(logical=FALSE)
  if(is.na(m) || m < 1L) m <- 1L
  ans <- max(1L, as.integer(m)-1L)

  if(ans > 1L) {
    tm <- timecheck(dta, TS, typeTS, TSextra)
    if(length(tm) > 1L) tm <- max(tm, na.rm=TRUE)

    if(!is.finite(tm) || 2*tm*B < 20 || B < 2*ans) {
      ans <- 1L
      if(!SuppressMessages)
        message("maxProcessor set to 1 because parallel overhead would dominate this short computation")
    }
    else if(!SuppressMessages)
      message("Using ", ans, " cores.")
  }

  ans
}

