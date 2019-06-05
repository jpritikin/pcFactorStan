#' Stan Models for the Pairwise Comparison Factor Model
#'
#' @docType package
#' @name pcFactorStan-package
#' @aliases pcFactorStan
#'
#' @description
#' \pkg{pcFactorStan} makes it easy to fit the pairwise comparison
#' factor model using \pkg{rstan}.
#' 
#' A user will generally want to use \code{\link{prepData}} and
#' \code{\link{pcStan}} to fit a model.
#'
#' The package includes a number of Stan models (see
#' \code{\link{findModel}} for a list) and an example dataset
#' \code{\link{phyActFlowPropensity}}.
#'
#' After gaining some experience with the pre-defined models, we
#' anticipate that users may write their own Stan models and fit them
#' with \code{\link[rstan]{stan}}, for which \code{\link{pcStan}} is a
#' wrapper.
#'
#' @useDynLib pcFactorStan, .registration = TRUE
#' @importFrom Rcpp loadModule
#' @import methods
#' 
NULL
