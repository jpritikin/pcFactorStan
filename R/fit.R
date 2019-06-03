#' Normalize data according to a canonical order
#'
#' @param df a data frame with pairs of vertices given in columns \code{pa1} and \code{pa2}
#' @param .palist a character vector giving an order to use instead of the default
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
normalizeData <- function(df, ..., .palist=NULL) {
  palist <- verifyIsData(df)
  if (!missing(.palist)) {
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
  df
}

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

  list(
    # all models
    NPA=length(palist),
    NCMP=length(pick),
    N=sum(weight),
    # multivariate models
    NITEMS=length(nthr),
    NTHRESH=nthr,
    TOFFSET=cumsum(nthr-1L) + 1L,
    # data (all models)
    pa1=pa1,
    pa2=pa2,
    weight=weight,
    pick=pick,
    refresh=refresh,
    # multivariate data
    item=item
  )
}

#' @export
prepData <- function(df) {
  df <- filterGraph(df)
  df <- normalizeData(df)
  prepCleanData(df)
}

#' @export
locateModel <- function(model, data) {
  extdata <- system.file("extdata", package = "pcFactorStan")

  if (missing(model)) {
    model <- ifelse(data$NITEMS == 1, "unidim+adapt", "covariance")
  }

  stan_path <- file.path(extdata, ifelse(grepl("\\.stan$", model), model, paste0(model, ".stan")))
  if(!file.exists(stan_path)) {
    avail <- sub(".stan$", "", dir(extdata, pattern=".stan$"), perl=TRUE)
    message("Models available: ", deparse(avail))
    stop(paste("Stan model not found:", stan_path))
  }

  if (model == 'unidim+adapt' && is.null(data[['varCorrection']])) {
    stop("You must choose a varCorrection. For example, data$varCorrection <- 2.0")
  } else if (is.null(data[['scale']])) {
    stop(paste0("You must choose a scale. For example, data$scale <- 1.8\n",
                "The 'unidim+adapt' model can help determine an optimal scale."))
  }

  stan_path
}

#' @export
pcStan <- function(model = "", data, ... ) {
  stan_path <- locateModel(model, data)
  message("Using ", file.path(stan_path))
  rstan::stan(stan_path, data = data, ...)
}