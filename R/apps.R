#' A Shiny app to experiment with the item response model
#'
#' @export
#' @examples
#' \dontrun{
#' itemModelExplorer()  # will launch a browser in RStudio
#' }
itemModelExplorer <- function() {
  want <- c("shiny","ggplot2")
  if (requireNamespace(want, quietly = TRUE)) {
    shiny::runApp(system.file('itemModelExplorer', package='pcFactorStan'))
  } else {
    stop(paste0("Please install.packages(",deparse(want),") and try again"))
  }
}
