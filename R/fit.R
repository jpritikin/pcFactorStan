#' Normalize data according to a canonical order
#'
#' @template args-df
#' @param ...  Not used.  Forces remaining arguments to be specified by name.
#' @param .palist a character vector giving an order to use instead of the default
#' @param .sortRows logical. Using the same order, sort rows in addition to vertex pairs.
#'
#' @description
#'
#' Pairwise comparison data are not commutative.
#' Alice beating Bob in chess is equivalent to Bob losing to
#' Alice. \code{normalizeData} assigns an arbitrary order to all
#' vertices and reorders vertices column-wise to match,
#' flipping signs as needed.
#'
#' @examples
#' df <- data.frame(pa1=NA, pa2=NA, i1=c(1, -1))
#' df[1,paste0('pa',1:2)] <- c('a','b')
#' df[2,paste0('pa',1:2)] <- c('b','a')
#' normalizeData(df)
#' @export
normalizeData <- function(df, ..., .palist=NULL, .sortRows=TRUE) {
  palist <- verifyIsData(df)
  if (!is.null(.palist)) {
    if (length(palist) != length(.palist)) {
      stop(paste(".palist must be length", length(palist)))
    }
    if (any(is.na(match(palist, .palist)))) {
      stop(".palist must contain the names of all vertices")
    }
    palist <- .palist
  }
  flip <- match(df$pa1, palist) > match(df$pa2, palist)
  dataCols <- -match(paste0('pa',1:2), colnames(df))
  df[flip, dataCols] <- -df[flip, dataCols]
  tmp <- df[flip, 'pa1']
  df[flip, 'pa1'] <- df[flip, 'pa2']
  df[flip, 'pa2'] <- tmp
  if (.sortRows) {
    df <- df[order(df$pa1, df$pa2),]
  }
  df
}

#' Transforms data into a form tailored for efficient evaluation by Stan
#'
#' Vertex names, if not already factors, are converted to
#' factors.  The number of thresholds per item is determined by the
#' largest absolute response value.  Missing responses are filtered
#' out.  Responses on the same pair of vertices on the same item are
#' grouped together.  Within a vertex pair and item, responses
#' are ordered from negative to positive.
#'
#' @details
#' Note: Reordering of responses is likely unless something like
#' \code{\link{normalizeData}} has been used with \code{.sortRows=TRUE}.
#'
#' @template args-df
#'
#' @return a data list suitable for passing as the \code{data}
#'   argument to \code{\link{pcStan}} or \code{\link[rstan]{stan}}
#' @family data preppers
#' @examples
#' df <- prepCleanData(phyActFlowPropensity)
#' str(df)
#' @importFrom reshape2 melt
#' @export
prepCleanData <- function(df) {
  palist <- verifyIsData(df)

  if (!is.factor(df$pa1)) {
    df$pa1 <- factor(df$pa1, levels=palist)
    df$pa2 <- factor(df$pa2, levels=palist)
  }

  dataCols <- -match(paste0('pa',1:2), colnames(df))
  nthr <- apply(df[,dataCols,drop=FALSE], 2, function(x) max(abs(x), na.rm=TRUE))
  names(nthr) <- NULL

  tall <- melt(df, id.vars = paste0('pa',1:2),
               variable.name = "item", value.name="pick")

  base <- 1L+max(length(nthr), length(palist))
  likeId <- (base * base * unclass(tall$pa1) +
             base * unclass(tall$pa2) + unclass(tall$item))

  l <- split(tall, likeId)
  pa1 <- c()
  pa2 <- c()
  item <- c()
  weight <- c()
  pick <- c()
  refresh <- c()
  for (lx in 1:length(l)) {
    x <- l[[lx]]
    pt <- table(x$pick, useNA="no")
    if (length(pt) == 0) next
    pa1 <- c(pa1, rep(unclass(x$pa1)[1], length(pt)))
    pa2 <- c(pa2, rep(unclass(x$pa2)[1], length(pt)))
    item <- c(item, rep(unclass(x$item)[1], length(pt)))
    weight <- c(weight, pt, use.names=FALSE)
    pick <- c(pick, as.integer(names(pt)))
    refresh <- c(refresh, c(1, rep(0, length(pt)-1)))
  }

  dl <- list(
    # all models
    NPA=length(palist),
    NCMP=length(pick),
    N=sum(weight),
    # multivariate models
    NITEMS=length(nthr),
    NTHRESH=nthr,
    TOFFSET=c(1L, 1L + cumsum(nthr)[-length(nthr)]),
    # data (all models)
    pa1=pa1,
    pa2=pa2,
    weight=weight,
    pick=pick,
    refresh=refresh,
    # multivariate data
    item=item
  )
  if (any(is.na(match(preppedDataFields, names(dl))))) {
    stop("Bug in prepCleanData(); contact developers")
  }
  dl
}

#' Transforms data into a form tailored for efficient evaluation by Stan
#'
#' @template args-df
#'
#' @description
#' Invokes \code{\link{filterGraph}} and \code{\link{normalizeData}}.
#' Vertex names, if not already factors, are converted to
#' factors.  The number of thresholds per item is determined by the
#' largest absolute response value.  Missing responses are filtered
#' out.  Responses on the same pair of vertices on the same item are
#' grouped together.  Within a vertex pair and item, responses
#' are ordered from negative to positive.
#'
#' @return a data list suitable for passing as the \code{data}
#'   argument to \code{\link{pcStan}} or \code{\link[rstan]{stan}}
#' @family data preppers
#' @examples
#' df <- prepData(phyActFlowPropensity)
#' str(df)
#'
#' @export
prepData <- function(df) {
  df <- filterGraph(df)
  df <- normalizeData(df)
  prepCleanData(df)
}

preppedDataFields <- c('NPA','NCMP','N','pa1','pa1','weight','pick','refresh')

verifyIsPreppedData <- function(data) {
  if (is.data.frame(data)) stop("Data must processed by prepData. Try prepData(data)")
  if (is.list(data) && all(!is.na(match(preppedDataFields, names(data))))) return()
  stop("Is data an object returned by prepData()? Seems not")
}

#' Given a model name and data, return the path to a Stan model
#'
#' @template args-locate
#' @description
#'
#' This is a convenience function to help you look up the path to an
#' appropriate model for your data.
#'
#' @details
#'
#' There are essentially three models: \sQuote{unidim}, \sQuote{covariance}, and \sQuote{factor}.
#' \sQuote{unidim} analyzes a single item. \sQuote{covariance} is suitable for two or more items.
#' Once you have vetted your items with the \sQuote{unidim} and \sQuote{covariance} models,
#' then you can try the \sQuote{factor} model.
#' For each model, there is a \sQuote{+ll} variation. This model
#' includes row-wise log likelihoods suitable for feeding to \pkg{loo}
#' for efficient approximate leave-one-out cross-validation.
#'
#' There is also a special model \sQuote{unidim+adapt}.  Except for
#' this model, the other models require a scaling constant.  To find
#' an appropriate scaling constant, we recommend fitting
#' \sQuote{unidim+adapt} to each item separately and then take the
#' median of median point estimates to set the scale. \sQuote{unidim+adapt} requires a
#' varCorrection constant. In general, a varCorrection of 2.0 or 3.0
#' should provide optimal results.
#'
#' @return a file system path to the selected model
#' @seealso \code{\link{pcStan}}
#' @examples
#' locateModel()  # shows available models
#'
#' df1 <- prepData(phyActFlowPropensity[,1:3])
#' df1$varCorrection <- 2.0
#' locateModel(data=df1)
#'
#' df2 <- prepData(phyActFlowPropensity[,1:5])
#' df2$scale <- 1.5
#' locateModel(data=df2)
#' @export
locateModel <- function(model=NULL, data=NULL) {
  if (!is.null(data)) {
    verifyIsPreppedData(data)
  }

  extdata <- system.file("extdata", package = "pcFactorStan")
  avail <- sub(".stan$", "", dir(extdata, pattern=".stan$"), perl=TRUE)

  if (is.null(model)) {
    if (is.null(data)) {
      message(paste("Models available:", paste0(deparse(avail), collapse="")))
      return(invisible())
    }
    model <- ifelse(data$NITEMS == 1, "unidim+adapt", "covariance")
  }

  stan_path <- file.path(extdata, ifelse(grepl("\\.stan$", model), model, paste0(model, ".stan")))
  if(!file.exists(stan_path)) {
    stop(paste("Stan model not found:", stan_path))
  }

  if (!is.null(data)) {
    if (model == 'unidim+adapt') {
      if (is.null(data[['varCorrection']])) {
        stop("You must choose a varCorrection. For example, data$varCorrection <- 2.0")
      }
    } else if (is.null(data[['scale']])) {
      stop(paste0("You must choose a scale. For example, data$scale <- 1.8\n",
                  "The 'unidim+adapt' model can help determine an optimal scale."))
    }
  }

  stan_path
}

#' Fit a pairwise comparison Stan model
#' @template args-locate
#' @template args-stan
#' @description Uses \code{\link{locateModel}} to find the appropriate
#'   model and then invokes \link[rstan]{stan}.
#' @return A \code{\link[rstan]{stanfit-class}} object.
#' @seealso See \code{\link[rstan]{stan}}, for which this function is
#'   a wrapper, for additional options. See \code{\link{prepData}} to
#'   create a suitable data list.  See
#'   \code{\link[rstan]{print.stanfit}} for ways of getting tables
#'   summarizing parameter posteriors.
#'
#' @return An object of S4 class \code{\link[rstan]{stanfit}}.
#' @seealso \code{\link{calibrateItems}}
#' @examples
#' dl <- prepData(phyActFlowPropensity[,c(1,2,3)])
#' dl$varCorrection <- 2.0
#' \dontrun{pcStan(data=dl)}
#' @importFrom rstan stan
#' @export
pcStan <- function(model=NULL, data, ...) {
  verifyIsPreppedData(data)
  stan_path <- locateModel(model, data)
  message("Using ", file.path(stan_path))
  rstan::stan(stan_path, data = data, ...)
}

#' Determine the optimal scale constant for a set of items
#' @template args-df
#' @template args-stan
#' @param iter A positive integer specifying the number of iterations for each chain (including warmup).
#' @param chains A positive integer specifying the number of Markov chains.
#' @param varCorrection A correction factor greater than or equal to 1.0
#' @param maxAttempts How many times to try re-running a model with more iterations.
#'
#' @description
#'
#' Data are passed through \code{\link{filterGraph}} and \code{\link{normalizeData}}.
#' Then the \sQuote{unidim+adapt} model is fit to each item individually.
#' A larger \code{varCorrection} will obtain a more accurate
#' \code{scale}, but is also more likely to produce an intractable
#' model. A good compromise is between 2.0 and 4.0.
#'
#' @seealso \code{\link[rstan]{check_hmc_diagnostics}}
#' @template ref-vehtari2019
#' @examples
#' \dontrun{
#' result <- calibrateItems(phyActFlowPropensity)
#' print(result)
#' }
#' @importMethodsFrom rstan summary
#' @export
calibrateItems <- function(df, iter=2000L, chains=4L, varCorrection=3.0, maxAttempts=5L, ...) {
  df <- filterGraph(df)
  df <- normalizeData(df)
  vCol <- match(paste0('pa',1:2), colnames(df))
  chains <- 4L
  varCorrection <- 3.0
  result <- expand.grid(item=colnames(df[,-vCol]),
                        iter=iter,
                        divergent=NA, treedepth=NA, low_bfmi=NA, n_eff=NA, Rhat=NA, scale=NA)
  for (attempt in 1:maxAttempts) {
    for (rx in 1:nrow(result)) {
      if (!is.na(result[rx,'divergent']) &&
          (result[rx,'divergent'] || result[rx,'treedepth'] || result[rx,'low_bfmi'])) next
      if (!is.na(result[rx, 'Rhat']) &&
          result[rx, 'Rhat'] < 1.015 && result[rx, 'n_eff'] > 100 * chains) next

      itemCol <- match(result[rx,'item'], colnames(df))
      dl <- prepCleanData(df[,c(vCol, itemCol)])
      dl$varCorrection <- varCorrection
      fit1 <- pcStan(data=dl, chains=chains, iter=result[rx,'iter'])
      result[rx,'divergent'] <- get_num_divergent(fit1)
      result[rx,'treedepth'] <- sum(get_max_treedepth_iterations(fit1))
      result[rx,'low_bfmi'] <- length(get_low_bfmi_chains(fit1))
      allPars <- summary(fit1, probs=0.5)$summary
      result[rx,'n_eff'] <- min(allPars[,'n_eff'])
      result[rx,'Rhat'] <- max(allPars[,'Rhat'])
      result[rx,'scale'] <- allPars['scale','50%']
      result[rx,'thetaVar'] <- var(summary(fit1, pars="theta", probs=c())$summary[,'mean'])
      result[rx,'iter'] <- round(result[rx,'iter'] * 1.5)
    }
  }
  result
}
