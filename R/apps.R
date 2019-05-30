#' A Shiny app to experiment with the item response model
#'
#' @export
#' @import shiny
#' @examples
#' \dontrun{
#' itemModelExplorer()  # will launch a browser in RStudio
#' }
itemModelExplorer <- function() {
	shiny::runApp(system.file('itemModelExplorer', package='pcFactorStan'))
}

.onAttach <- function(libname, pkgname) {
	packageStartupMessage(paste(
          "Use itemModelExplorer() to explore the item model."))
}
