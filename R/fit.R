#' Normalize data according to a canonical order
#'
#' @template args-df
#' @template args-dots-barrier
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
  if (length(list(...)) > 0) stop("Rejected are any values passed in the '...' argument")

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
#' @template return-datalist
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
  itemNames <- names(nthr)
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
    nameInfo=list(pa=palist, item=itemNames),
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
    stop("Bug in prepCleanData(); contact developers") # nocov
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
#' @template return-datalist
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

preppedDataFields <- c('nameInfo','NPA','NCMP','N','pa1','pa1','weight',
                       'pick','refresh')

verifyIsPreppedData <- function(data) {
  if (is.data.frame(data)) stop("Data must be processed by prepData. Try prepData(data)")
  if (is.list(data) && all(!is.na(match(preppedDataFields, names(data))))) return()
  stop("Is data an object returned by prepData()? Seems not")
}

#' Given a model name, return stanmodel object
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
#' For each model, there is a \sQuote{_ll} variation. This model
#' includes row-wise log likelihoods suitable for feeding to \pkg{loo}
#' for efficient approximate leave-one-out cross-validation (Vehtari, Gelman, & Gabry, 2017).
#'
#' There is also a special model \sQuote{unidim_adapt}.  Except for
#' this model, the other models require a scaling constant.  To find
#' an appropriate scaling constant, we recommend fitting
#' \sQuote{unidim_adapt} to each item separately and then take the
#' median of median point estimates to set the scale. \sQuote{unidim_adapt} requires a
#' varCorrection constant. In general, a varCorrection of 2.0 or 3.0
#' should provide optimal results.
#'
#' TODO discuss factor model
#'
#' @return An instance of S4 class \code{\link[rstan:stanmodel-class]{stanmodel}} that can be passed to \code{\link{pcStan}}.
#' @seealso \code{\link{toLoo}}
#' @template ref-vehtari2017
#' @examples
#' findModel()  # shows available models
#' findModel('unidim')
#' @export
findModel <- function(model=NULL) {
  avail <- names(stanmodels)

  if (is.null(model)) {
    message(paste("Models available:", paste0(deparse(avail), collapse="")))
    return(invisible())
  }

  obj <- stanmodels[[model]]
  if(is.null(obj)) {
    stop(paste("Stan model not found:", model))
  }

  obj
}

#' Specify a single factor model
#'
#' Specify a single latent factor with a path to each item.
#'
#' @param factorScalePrior standard deviation of the normal prior for the logit transformed factor proportion
#' @template args-data
#' @template return-datalist
#' @examples
#' dl <- prepData(phyActFlowPropensity)
#' dl <- prepSingleFactorModel(dl, 0.9)
#' str(dl)
#' @family factor model
#' @family data preppers
#' @export
prepSingleFactorModel <- function(data, factorScalePrior) {
  verifyIsPreppedData(data)
  data$factorScalePrior <- as.array(factorScalePrior)
  ni <- data$NITEMS
  data$factorItemPath <- matrix(c(rep(1,ni), 1:ni), nrow=2, byrow=TRUE)
  data$NFACTORS <- 1L
  data$NPATHS <- ni
  data
}

#' Specify a factor model
#'
#' Specify a factor model with an arbitrary number of factors and
#' arbitrary factor-to-item structure.
#'
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
#'
#' @param path a named list of item names
#' @param factorScalePrior named numeric vector
#' @template args-data
#' @template return-datalist
#' @examples
#' pa <- phyActFlowPropensity[,setdiff(colnames(phyActFlowPropensity),
#'                                     c('goal1','feedback1'))]
#' dl <- prepData(pa)
#' dl <- prepFactorModel(dl,
#'                       list(flow=c('complex','skill','predict',
#'                                   'creative', 'novelty', 'stakes',
#'                                   'present', 'reward', 'chatter',
#'                                   'body'),
#'                            f2=c('waiting','control','evaluated','spont'),
#'                            rc=c('novelty', 'waiting')),
#'                       c(flow=0.9, f2=0.5, rc=0.2))
#' str(dl)
#' @family factor model
#' @family data preppers
#' @export
prepFactorModel <- function(data, path, factorScalePrior) {
  verifyIsPreppedData(data)
  if (length(path) == 0) stop("No paths specified")
  n1 <- names(path)
  n2 <- names(factorScalePrior)
  if (length(n1) == 0) {
    stop("paths must be named, the name is the name of the factor")
  }
  if (length(n2) == 0) {
    stop("factorScalePrior must be named, the name is the name of the factor")
  }
  if (length(n1) != length(n2)) {
    stop(paste("Number of factors mismatch between path",
               length(path),"and factorScalePrior",
               length(factorScalePrior)))
  }
  if (!setequal(n1, n2)) {
    stop("path and factorScalePrior specify different factor names")
  }
  itemsPerFactor <- sapply(path, length)
  if (any(itemsPerFactor < 2)) {
    stop(paste("Some factors have less than 2 indicators:",
               paste(names(itemsPerFactor)[itemsPerFactor<2],
                     collapse=", ")))
  }
  items <- data$nameInfo$item
  indicators <- Reduce(union, path, c())
  noItem <- is.na(match(indicators, items))
  if (any(noItem)) {
    stop(paste("No matching item for factor indicator(s):",
               paste(indicators[noItem], collapse=", ")))
  }
  noIndicator <- is.na(match(items, indicators))
  if (any(noIndicator)) {
    stop(paste("No factor predicts item(s):",
               paste(items[noIndicator], collapse=", ")))
  }
  data$factorScalePrior <- as.array(unname(factorScalePrior[n1]))
  data$factorItemPath <- matrix(c(rep(1:length(itemsPerFactor), itemsPerFactor),
                                  unlist(lapply(path, function(x) match(x, items)))),
                                nrow=2, byrow=TRUE)
  data$NFACTORS <- length(factorScalePrior)
  data$NPATHS <- sum(itemsPerFactor)
  data$nameInfo[['factor']] <- names(path)
  data
}

#' Fit a paired comparison Stan model
#' @template args-locate
#' @template args-data
#' @template args-stan
#' @description Uses \code{\link{findModel}} to find the appropriate
#'   model and then invokes \link[rstan:sampling]{sampling}.
#' @return A \code{\link[rstan:stanfit-class]{stanfit}} object.
#' @seealso See \code{\link[rstan:sampling]{sampling}}, for which this function is
#'   a wrapper, for additional options. See \code{\link{prepData}} to
#'   create a suitable data list.  See
#'   \code{\link[rstan:print.stanfit]{print.stanfit}} for ways of getting tables
#'   summarizing parameter posteriors.
#'
#' @return An object of S4 class \code{\link[rstan:stanfit-class]{stanfit}}.
#' @seealso \code{\link{calibrateItems}}, \code{\link{outlierTable}}
#' @examples
#' dl <- prepData(phyActFlowPropensity[,c(1,2,3)])
#' dl$varCorrection <- 2.0
#' \donttest{pcStan('unidim_adapt', data=dl)}  # takes more than 5 seconds
#' @importFrom rstan sampling
#' @export
pcStan <- function(model, data, ...) {
  verifyIsPreppedData(data)
  obj <- findModel(model)
  rstan::sampling(obj, data = data, ...)
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
#' Then the \sQuote{unidim_adapt} model is fit to each item individually.
#' A larger \code{varCorrection} will obtain a more accurate
#' \code{scale}, but is also more likely to produce an intractable
#' model. A good compromise is between 2.0 and 4.0.
#'
#' @return
#' A data.frame (one row per item) with the following columns:
#' \describe{
#'  \item{item}{Name of the item}
#'  \item{iter}{Number of iterations per chain}
#'  \item{divergent}{Number of divergent transitions observed after warmup}
#'  \item{treedepth}{Number of times the treedepth was exceeded}
#'  \item{low_bfmi}{Number of chains with low E-BFMI}
#'  \item{n_eff}{Minimum effective number of samples across all parameters}
#'  \item{Rhat}{Maximum Rhat across all parameters}
#'  \item{scale}{Median marginal posterior of \code{scale}}
#'  \item{thetaVar}{Median variance of theta (latent scores)}
#'  }
#' @seealso \code{\link[rstan:check_hmc_diagnostics]{check_hmc_diagnostics}}
#' @template ref-vehtari2019
#' @examples
#' \donttest{
#' result <- calibrateItems(phyActFlowPropensity)  # takes more than 5 seconds
#' print(result)
#' }
#' @importMethodsFrom rstan summary
#' @importFrom rstan get_num_divergent get_max_treedepth_iterations get_low_bfmi_chains
#' @export
calibrateItems <- function(df, iter=2000L, chains=4L, varCorrection=3.0, maxAttempts=5L, ...) {
  df <- filterGraph(df)
  df <- normalizeData(df)
  vCol <- match(paste0('pa',1:2), colnames(df))
  result <- expand.grid(item=colnames(df[,-vCol]),
                        iter=NA,
                        divergent=NA, treedepth=NA, low_bfmi=NA, n_eff=NA, Rhat=NA,
                        scale=NA, thetaVar=NA)
  for (attempt in 1:maxAttempts) {
    for (rx in 1:nrow(result)) {
      if (!is.na(result[rx,'divergent']) &&
          (result[rx,'divergent'] || result[rx,'treedepth'] || result[rx,'low_bfmi'])) next
      if (!is.na(result[rx, 'Rhat']) &&
          result[rx, 'Rhat'] < 1.015 && result[rx, 'n_eff'] > 100 * chains) next

      itemCol <- match(result[rx,'item'], colnames(df))
      dl <- prepCleanData(df[,c(vCol, itemCol)])
      dl$varCorrection <- varCorrection
      result[rx,'iter'] <- ifelse(is.na(result[rx,'iter']), iter, result[rx,'iter'] * 1.5)
      fit1 <- suppressWarnings(pcStan("unidim_adapt", data=dl, chains=chains, iter=result[rx,'iter'], ...))
      result[rx,'divergent'] <- get_num_divergent(fit1)
      result[rx,'treedepth'] <- sum(get_max_treedepth_iterations(fit1))
      result[rx,'low_bfmi'] <- length(get_low_bfmi_chains(fit1))
      allPars <- summary(fit1, probs=0.5)$summary
      result[rx,'n_eff'] <- min(allPars[,'n_eff'])
      result[rx,'Rhat'] <- max(allPars[,'Rhat'])
      result[rx,'scale'] <- allPars['scale','50%']
      result[rx,'thetaVar'] <- allPars['thetaVar','50%']
    }
  }
  result
}

#' Compute approximate leave-one-out (LOO) cross-validation for Bayesian
#' models using Pareto smoothed importance sampling (PSIS)
#'
#' @template args-fit
#' @param ... Additional options passed to \code{\link[loo:loo]{loo}}.
#'
#' @description
#'
#' You must use an \sQuote{_ll} model variation (see \code{\link{findModel}}).
#'
#' @return a loo object
#' @seealso \code{\link{outlierTable}}, \code{\link[loo:loo]{loo}}
#' @importFrom loo loo extract_log_lik relative_eff
#' @export
#' @examples
#' palist <- letters[1:10]
#' df <- twoLevelGraph(palist, 300)
#' theta <- rnorm(length(palist))
#' names(theta) <- palist
#' df <- generateItem(df, theta, th=rep(0.5, 4))
#'
#' df <- filterGraph(df)
#' df <- normalizeData(df)
#' dl <- prepCleanData(df)
#' dl$scale <- 1.5
#'
#' \donttest{
#' m1 <- pcStan("unidim_ll", dl)
#'
#' loo1 <- toLoo(m1, cores=1)
#' print(loo1)
#' }
toLoo <- function(fit, ...) {
  ll <- extract_log_lik(fit, merge_chains = FALSE)
  loo(ll, r_eff=relative_eff(exp(ll)), ...)
}

#' List observations with Pareto values larger than a given threshold
#'
#' @template args-data
#' @param x An object created by \code{\link[loo:loo]{loo}}
#' @param threshold threshold is the minimum k value to include
#'
#' @description
#'
#' The function \code{\link{prepCleanData}} compresses observations
#' into the most efficient format for evaluation by Stan. This function
#' maps indices of observations back to the actual observations,
#' filtering by the largest Pareto k values. It is assumed that
#' \code{data} was processed by \code{\link{normalizeData}} or is in
#' the same order as seen by \code{\link{prepCleanData}}.
#'
#' @return
#' A data.frame (one row per observation) with the following columns:
#' \describe{
#' \item{pa1}{Name of object 1}
#' \item{pa2}{Name of object 2}
#' \item{item}{Name of item}
#' \item{pick}{Observed response}
#' \item{k}{Associated Pareto k value}
#' }
#' @seealso \code{\link{toLoo}}, \code{\link[loo:pareto-k-diagnostic]{pareto_k_ids}}
#' @importFrom loo pareto_k_ids pareto_k_values
#' @export
#' @examples
#' palist <- letters[1:10]
#' df <- twoLevelGraph(palist, 300)
#' theta <- rnorm(length(palist))
#' names(theta) <- palist
#' df <- generateItem(df, theta, th=rep(0.5, 4))
#'
#' df <- filterGraph(df)
#' df <- normalizeData(df)
#' dl <- prepCleanData(df)
#' dl$scale <- 1.5
#'
#' \donttest{
#' m1 <- pcStan("unidim_ll", dl)
#'
#' loo1 <- toLoo(m1, cores=1)
#' ot <- outlierTable(dl, loo1, threshold=.2)
#' df[df$pa1==ot[1,'pa1'] & df$pa2==ot[1,'pa2'], 'i1']
#' }
outlierTable <- function(data, x, threshold=0.5) {
  verifyIsPreppedData(data)
  ids <- pareto_k_ids(x, threshold)
  loc <- c(1L, 1L+cumsum(data$weight))
  offset <- findInterval(ids, loc)
  palist <- data$nameInfo$pa
  itemName <- data$nameInfo$item
  df <- data.frame(pa1=data$pa1[offset],
             pa2=data$pa2[offset],
             item=data$item[offset],
             pick=data$pick[offset],
             k=pareto_k_values(x)[ids])
  df <- df[!duplicated(offset),]
  for (k in paste0('pa',1:2)) {
    levels(df[[k]]) <- palist
    class(df[[k]]) <- 'factor'
  }
  levels(df[['item']]) <- itemName
  class(df[['item']]) <- 'factor'
  df <- df[order(-df$k),]
  rownames(df) <- c()  # original order is meaningless
  df
}

assertDataFitCompat <- function(dl, fit) {
  if (!is(dl, 'list')) stop("dl must be a list of data")
  if (!is(fit, 'stanfit')) stop("fit must be a stanfit object")
  pd <- fit@par_dims
  if (pd$theta[1] != dl$NPA) {
    stop(paste0("dl has ",dl$NPA," objects but fit has ",pd$theta[1]," objects"))
  }
  fitItems <- ifelse(length(pd$theta)==1, 1, pd$theta[2])
  if (!fitItems == dl$NITEMS) {
    stop(paste0("dl has ",dl$NITEMS," items but fit has ",fitItems," items"))
  }
  if (pd$threshold != sum(dl$NTHRESH)) {
    stop(paste0("dl has a total of ",dl$NTHRESH," thresholds across all items but fit has ",
                pd$threshold," thresholds"))
  }
}

#' Produce data suitable for plotting item response curves
#'
#' @template args-dl
#' @template args-fit
#' @param item a vector of item names
#' @param responseNames a vector of labels for the possible responses
#' @template args-samples
#' @param from the starting latent difference value
#' @param to the ending latent difference value
#' @param by the grid increment
#'
#' @description
#' Selects \code{samples} random draws from the posterior and evaluates the item
#' response curve on the grid given by \code{seq(from,to,by)}.
#' All items use the same \code{responseNames}. If you have some items
#' with a different number of thresholds or different response names
#' then you can call \code{responseCurve} for each item separately
#' and \code{rbind} the results together.
#'
#' @template detail-response
#' @template ref-masters1982
#' @return
#' A data.frame with the following columns:
#' \describe{
#' \item{response}{Which response}
#' \item{worthDiff}{Difference in worth}
#' \item{item}{Which item}
#' \item{sample}{Which sample}
#' \item{prob}{Associated probability}
#' \item{responseSample}{A grouping index for independent item response samples}
#' }
#' @export
#' @importFrom rstan extract
#' @family data extractor
#' @examples
#' \donttest{ vignette('manual', 'pcFactorStan') }
responseCurve <- function(dl, fit, responseNames, item=dl$nameInfo$item,
                          samples=100, from=-6, to=-from, by=.1) {
  assertDataFitCompat(dl, fit)
  pd <- fit@par_dims
  itemIndex <- match(item, dl$nameInfo$item)
  if (any(is.na(itemIndex))) {
    stop(paste0("Item not found: ",
                paste(item[is.na(itemIndex)], collapse=', ')))
  }
  mismatch <- (1 + 2 * dl$NTHRESH[itemIndex]) != length(responseNames)
  if (any(mismatch)) {
    stop(paste0("Item with different number of responseNames: ",
                paste(item[mismatch], collapse=', ')))
  }
  grid <- seq(from,to,by)
  df <- expand.grid(response=responseNames, worthDiff=grid,
                    item=item, sample=1:samples, prob=NA)
  pick <- c()
  for (i1 in item) {
    ii <- match(i1, dl$nameInfo$item)
    thrInd <- dl$TOFFSET[ii] + 1:dl$NTHRESH[ii] - 1L
    thrData <- extract(fit, pars=paste0("threshold[",thrInd,"]"))
    if ('alpha' %in% names(pd)) {
      if (length(pd$alpha) == 0) {
        alphaData <- extract(fit, pars="alpha")[[1]]
      } else {
        alphaData <- extract(fit, pars=paste0("alpha[",ii,"]"))[[1]]
      }
    } else {
      alphaData <- dl$alpha[ii]
    }
    if (length(pick)==0) {
      samples <- min(length(thrData[[1]]), samples)
      pick <- sample.int(length(thrData[[1]]), samples)
    }
    for (sx in 1:samples) {
      mask <- df$item == i1 & df$sample == sx
      p1 <- sapply(grid, function(gx) {
        if (length(alphaData) > 1) {
          alpha <- alphaData[pick[sx]]
        } else {
          alpha <- alphaData
        }
        cmp_probs(dl$scale[ii], alpha, 0, gx,
                  sapply(thrData, function(x) x[pick[sx]]))
      })
      df[mask, 'prob'] <- c(p1)
    }
  }
  df$responseSample <-
    (match(df$item, dl$nameInfo$item) * length(responseNames) * samples +
     df$sample * length(responseNames) +
     unfactor(df$response)-1L)
  df
}

#' Remove the array indexing from a parameter name
#' @param name a parameter name
#' @return the name without the square bracket parameter indexing
#' @examples
#' withoutIndex("foo[1,2")
#' @export
withoutIndex <- function(name) {
  sub("\\[[\\d,]+\\]$", "", name, perl=TRUE)
}

#' Produce data suitable for plotting parameter estimates
#' 
#' @template args-fit
#' @param pars a vector of parameter names
#' @param label column name for \code{nameVec}
#' @param nameVec a vector of explanatory parameters names
#' @param width a width in probability units for the uncertainty interval
#'
#' @return
#' A data.frame with the following columns:
#' \describe{
#' \item{L}{Lower quantile}
#' \item{M}{Median}
#' \item{U}{Upper quantile}
#' \item{\emph{label}}{\emph{nameVec}}
#' }
#' @export
#' @importFrom rstan summary
#' @family data extractor
#' @examples
#' \donttest{ vignette('manual', 'pcFactorStan') }
parInterval <- function(fit, pars, nameVec,
                        label=withoutIndex(pars[1]), width=0.8) {
  probs <- 0.5 + c(-width/2, 0, width/2)
  interval <- summary(fit, pars=pars, probs=probs)$summary[,4:6,drop=FALSE]
  colnames(interval) <- c("L","M","U")
  interval <- as.data.frame(interval)
  interval[[label]] <- factor(nameVec, levels=nameVec)
  interval
}

#' Produce data suitable for plotting parameter distributions
#' 
#' @template args-fit
#' @param pars a vector of parameter names
#' @param label column name for \code{nameVec}
#' @param nameVec a vector of explanatory parameters names
#' @template args-samples
#' @return
#' A data.frame with the following columns:
#' \describe{
#' \item{sample}{Sample index}
#' \item{\emph{label}}{A name from \emph{nameVec}}
#' \item{value}{A single sample of the associated parameter}
#' }
#' @export
#' @importFrom rstan extract
#' @importFrom reshape2 melt
#' @family data extractor
#' @examples
#' \donttest{ vignette('manual', 'pcFactorStan') }
parDistributionCustom <- function(fit, pars, nameVec,
                                  label=withoutIndex(pars[1]), samples=500) {
  colSet <- extract(fit, pars)
  nextCol <- 1L
  pick <- c()
  tall <- NULL
  for (c1 in colSet) {
    if (length(dim(c1)) == 1) c1 <- as.matrix(c1)
    colnames(c1) <- nameVec[nextCol:(nextCol + ncol(c1) - 1L)]
    nextCol <- nextCol + ncol(c1)
    if (length(pick) == 0) {
      samples <- min(nrow(c1), samples)
      pick <- sample.int(nrow(c1), samples)
    }
    c1 <- c1[pick,,drop=FALSE]
    tall1 <- melt(c1)
    colnames(tall1)[1:2] <- c('sample',label)
    tall <- rbind(tall, tall1)
  }
  tall
}

#' @param pi a data.frame returned by \code{\link{parInterval}}
#' @rdname parDistributionCustom
#' @export
parDistributionFor <- function(fit, pi, samples=500) {
  parDistributionCustom(fit, rownames(pi), nameVec=pi[[4]], label=colnames(pi)[4])
}
