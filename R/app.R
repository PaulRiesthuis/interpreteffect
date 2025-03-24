#' Launch Shiny App
#'
#' @param name The name of the app to run
#' @param ... arguments to pass to shiny::runApp
#'
#' @export
#'
app <- function(name = "app", ...) {
  baseDir <- system.file("apps", package = "interpreteffect")

  if (baseDir == "") {
    stop("The 'apps' directory does not exist in the interpreteffect package.")
  }

  appPath <- file.path(baseDir, paste0(name, ".R"))  # Check for a single file

  if (file.exists(appPath)) {
    shiny::runApp(appPath, ...)
  } else {
    stop("The shiny app '", name, "' does not exist in interpreteffect.")
  }
}
