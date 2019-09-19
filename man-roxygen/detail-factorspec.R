#' @details
#' For each factor, you need to specify its name, which items it predicts, and its scale prior.
#' The connections from factors to items is specified by the `path` argument.
#' The scale priors are given in the `factorScalePrior` argument.
#' Both factors and items are specified by name (not index).
#' The example shows how everything fits together.
#' Paths are ordered as given in the `path` argument.
#' 
#' The units of `factorScalePrior` is a standard deviation of the
#' normal prior for the logit transformed factor proportion.
