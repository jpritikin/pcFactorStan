#' @section Response model:
#' 
#' The paired comparison item response model has thresholds and a
#' scale parameter similar to the partial credit model (Masters,
#' 1982). The model is cumbersome to describe in traditional
#' mathematical notation, but the R code is fairly straightforward,
#'
#' \preformatted{
#' softmax <- function(y) exp(y) / sum(exp(y))
#' 
#' cmp_probs <- function(scale, alpha, pa1, pa2, thRaw) {
#'   th <- cumsum(thRaw)
#'   diff <- scale * (pa2 - pa1)
#'   unsummed <- c(0, diff + rev(th), diff - th, use.names = FALSE)
#'   softmax(cumsum(alpha * unsummed))
#' }
#' }
#'
#' The function \code{cmp_probs} takes a \code{scale} constant,
#' \code{alpha} discrimination, the latent scores for two objects
#' \code{pa1} and \code{pa2}, and a vector of thresholds
#' \code{thRaw}. The thresholds are parameterized as the difference
#' from the previous threshold. For example, thresholds c(0.5, 0.5)
#' are not at the same location but are at locations c(0.5,
#' 1.0). Thresholds are symmetric. If there is one threshold then the
#' model admits three possible response outcomes (e.g. \emph{win}, \emph{tie}, and
#' \emph{lose}). Responses are always stored centered with zero representing
#' a tie. Therefore, it is necessary to add one plus the number of
#' thresholds to response data to index into the vector returned by
#' \code{cmp_probs}. For example, if our response data is (-1, 0, 1)
#' and has one threshold then we would add 2 (1 + 1 threshold) to
#' obtain the indices (1, 2, 3).
#' 
#' Use \code{\link{itemModelExplorer}} to explore the item model. In
#' this \pkg{shiny} app, the \emph{discrimination} parameter does what
#' is customary in item response models. However, it is not difficult
#' to show that discrimination is a function of thresholds and
#' scale. That is, discrimination is not an independent parameter.  In
#' paired comparison models, discrimination and measurement error are
#' confounded.
