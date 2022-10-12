#' run_shiny
#'
#' Runs the shiny app associated with R2sample package
#' @return No return value, called for side effect of opening a shiny app
#' @export
run_shiny <- function() {
  appDir <- system.file("shiny-example", "r2sample_shiny", package = "R2sample")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `R2sample`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}
