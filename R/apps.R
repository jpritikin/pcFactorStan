#' A Shiny app to experiment with the item response model
#'
#' @export
#' @examples
#' \dontrun{
#' itemModelExplorer()  # will launch a browser in RStudio
#' }
itemModelExplorer <- function() {
  library(shiny)
  shiny::runApp(system.file('itemModelExplorer', package='pcFactorStan'))
}

.onAttach <- function(libname, pkgname) {
	packageStartupMessage(paste(
          "Use itemModelExplorer() to explore the item model."))
}
