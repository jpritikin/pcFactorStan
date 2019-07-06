#' A Shiny app to experiment with the item response model
#'
#' @template args-dl
#' @template args-fit
#' @param item name of the item to visualize
#' @description
#'
#' When data \code{dl} and fitted model \code{fit} are provided, the
#' item parameters associated with \code{item} are loaded for
#' inspection.
#'
#' @template detail-response
#' @template ref-masters1982
#' 
#' @importFrom rstan summary
#' @export
#' @examples
#' \donttest{
#' itemModelExplorer()  # will launch a browser in RStudio
#' }
itemModelExplorer <- function(dl=NULL, fit=NULL, item=NULL) { # nocov start
  want <- c("shiny","ggplot2")
  if (requireNamespace(want, quietly = TRUE)) {
    if (!is.null(dl)) {
      assertDataFitCompat(dl, fit)
      ii <- match(item, dl$nameInfo$item)
      if (is.na(ii)) stop(paste("Cannot find item", item))
      thrInd <- dl$TOFFSET[ii] + 1:dl$NTHRESH[ii] - 1L

      if ("alpha" %in% names(fit@par_dims)) {
        alphaName <- 'alpha'
        if (length(fit@par_dims$alpha)) alphaName <- paste0("alpha[",ii,"]")
        alpha <- summary(fit, pars=alphaName, .5)$summary[,'50%']
      } else if (!is.null(dl[['alpha']])) {
        alpha <- dl[['alpha']][ii]
      } else {
        stop("Where to find alpha?")
      }

      itemModelExplorer.default <-
        list(scaleValue = dl$scale[ii],
             discrimination = alpha,
             thresholds = summary(fit, pars=paste0("threshold[",thrInd,"]"))$summary[,'mean'])
    }

    tryCatch(shiny::runApp(system.file('itemModelExplorer', package='pcFactorStan')),
             finally=rm("itemModelExplorer.default", envir=globalenv()))

  } else {
    stop(paste0("Please install.packages(",deparse(want),") and try again"))
  }
} # nocov end
