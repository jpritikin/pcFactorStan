#' Item response function for pairwise comparisons
#'
#' Use \code{\link{itemModelExplorer}} to explore the item model. In
#' this \pkg{shiny} app, the \emph{discrimination} parameter does what
#' is customary in item response models. However, it is not difficult
#' to show that discrimination is a function of thresholds and
#' scale. That is, discrimination is not an independent parameter.  In
#' paired comparison models, discrimination and measurement error are
#' confounded.
#'
#' @details
#'
#' The thresholds are parameterized as the difference
#' from the previous threshold. For example, thresholds c(0.5, 0.6)
#' are not at the same location but are at locations c(0.5,
#' 1.1). Thresholds are symmetric. If there is one threshold then the
#' model admits three possible response outcomes (e.g. \emph{win}, \emph{tie}, and
#' \emph{lose}). Responses are always stored centered with zero representing
#' a tie. Therefore, it is necessary to add one plus the number of
#' thresholds to response data to index into the vector returned by
#' \code{cmp_probs}. For example, if our response data is (-1, 0, 1)
#' and has one threshold then we would add 2 (1 + 1 threshold) to
#' obtain the indices (1, 2, 3).
#'
#' @section Math:
#' Up until version 1.4, the item response model was based on the
#' partial credit model (Masters, 1982). In version 1.5,
#' the graded response model is used instead (Samejima, 1969).
#' The advantage of the graded response model is greater
#' independence among threshold parameters and the ability to
#' compute only the parts of the model that are actually needed
#' given particular observations. The curves predicted by both
#' models are similar and should obtain similar results in data
#' analyses.
#'
#' @param alpha discrimination parameter
#' @param scale scale correction factor
#' @param pa1 first latent worth
#' @param pa2 second latent worth
#' @param thRaw vector of positive thresholds
#' @template ref-samejima1969
#' @template ref-masters1982
#' @export
#' @return A vector of probabilities of observing each outcome
#' @examples
#' # Returns probabilities of
#' # c(pa1 > pa2, pa1 = pa2, pa1 < pa2)
#' cmp_probs(1,1,0,1,.8)
#'
#' # Add another threshold for a symmtric 3 point Likert scale
#' cmp_probs(1,1,0,.5,c(.8, 1.6))
cmp_probs <- function(alpha, scale, pa1, pa2, thRaw) {
  th <- cumsum(abs(thRaw))
  diff <- scale * (pa1 - pa2)
  at <- c(diff - rev(th), diff + th)
  #  pr <- plogis(at, scale=1.0/alpha)
  pr <- 1/(1+exp(-at*alpha))
  pr <- c(0, pr, 1)
  diff(pr)
}

pairMap <- function(n) {
  result <- matrix(0, 2, n*(n-1)/2)
  col <- 1
  for (rx in 2:n) {
    for (cx in 1:(rx-1)) {
      result[,col] <- c(cx,rx)
      col <- col + 1L
    }
  }
  result
}

verifyIsData <- function(df) {
  if (!is.data.frame(df)) stop("df is not a data.frame")
  if (any(is.na(match(paste0('pa',1:2), colnames(df))))) {
    stop("Expected columns pa1, pa2 to contain the vertex names")
  }
  if (is.factor(df$pa1) && is.factor(df$pa2)) {
    if (length(levels(df$pa1)) != length(levels(df$pa2)) ||
        any(levels(df$pa1) != levels(df$pa2))) {
      stop("Some levels(pa1) don't match levels(pa2)")
    }
    levels(df$pa1)
  } else {
    unique(c(as.character(df$pa1), as.character(df$pa2)))
  }
}

assertNameUnused <- function(df, name) {
  collision <- !is.na(match(name, colnames(df)))
  if (any(collision)) {
    stop(paste("Colname", paste(name[collision], collapse=' '),
               "is already taken. Please provide name="))
  }
}

#' Generate paired comparison data with a common factor that
#' accounts for some proportion of the variance
#'
#' @param prop the number of items or a vector of signed proportions of variance
#' @inheritParams generateItem
#' @description
#'
#' Imagine that there are people that play in tournaments of more than
#' one board game. For example, the computer player AlphaZero (Silver
#' et al. 2018) has trained to play chess, shogi, and Go. We can take
#' the tournament match outcome data and find rankings among the
#' players for each of these games. We may also suspect that there is
#' a latent board game skill that accounts for some proportion of the
#' variance in the per-board game rankings.
#'
#' @template detail-response
#' @template detail-genfactor
#'
#' @template ref-silver2018
#' @return
#' The given data.frame \code{df} with additional columns for each item.
#' @family item generators
#' @examples
#' df <- twoLevelGraph(letters[1:10], 100)
#' df <- generateSingleFactorItems(df, 3)
#' @export
#' @importFrom stats rnorm sd rbinom rbeta
generateSingleFactorItems <- function(df, prop, th=0.5, name, ..., scale=1, alpha=1) {
  palist <- verifyIsData(df)
  if (length(prop) == 1) {
    if (prop < 3) stop(paste0("At least 3 indicators are required (", prop," requested)"))
    prop <- rbeta(prop, 4, 3)
    for (cx in 2:length(prop)) {
      if (rbinom(1, 1, .5)) prop[cx] <- -prop[cx]
    }
  }
  if (length(prop) < 3) stop(paste0("At least 3 indicators are required (", length(prop)," given)"))
  if (any(prop <= -1 | prop >= 1)) stop("Signed proportions must be between -1 and 1")
  if (missing(name)) {
    num <- ncol(df)-1
    name <- paste0('i', num:(num+length(prop)-1))
  }
  assertNameUnused(df, name)
  posProp <- abs(prop)
  factorScore <- c(scale(rnorm(length(palist))))
  thetaF <- factorScore %*% t(sqrt(posProp/(1-posProp)))
  thetaF <- t(t(thetaF) * sign(prop))
  thetaN <- matrix(rnorm(length(palist)*length(prop)),
                   length(palist), length(prop))
  thetaN <- scale(thetaN)

  # to double check
  ## varF <- apply(thetaF, 2, var)
  ## varN <- apply(thetaN, 2, var)
  ## varF / (varF + varN) - posProp

  theta <- scale(thetaN + thetaF)
  dimnames(theta) <- list(palist, name)
  generateItem(df, theta, th, scale=scale, alpha=alpha)
}

validateFactorModel <- function(items, path) {
  if (length(path) == 0) stop("No paths specified")
  n1 <- names(path)
  if (length(n1) == 0) {
    stop("paths must be named, the name is the name of the factor")
  }
  itemsPerFactor <- sapply(path, length)
  if (any(itemsPerFactor < 2)) {
    stop(paste("Some factors have less than 2 indicators:",
               paste(names(itemsPerFactor)[itemsPerFactor<2],
                     collapse=", ")))
  }
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
}

ssqrt <- function(v) sign(v)*sqrt(abs(v))

#' Generate paired comparison data for a factor model
#'
#' @template args-path
#' @template args-factorScalePrior
#' @inheritParams generateItem
#' @description
#' Generate paired comparison data given a mapping from factors to items.
#'
#' @template detail-factorspec
#'
#' @details Path proportions (factor-to-item loadings) are sampled
#'   from a logistic transformed normal distribution with scale
#'   0.6. A few attempts are made to resample path
#'   proportions if any of the item proportions sum to more than
#'   1.0. An exception will be raised if repeated attempts fail to
#'   produce viable proportion assignments.
#'
#' @template detail-response
#' @template detail-genfactor
#'
#' @template ref-silver2018
#' @family item generators
#' @return
#' The given data.frame \code{df} with additional columns for each item.
#' In addition, you can obtain path proportions (factor-to-item loadings) from \code{attr(df, "pathProp")},
#' the factor scores from \code{attr(df, "score")},
#' and latent worths from \code{attr(df, "worth")}.
#'
#' @seealso To fit a factor model: \link{prepFactorModel}
#' @examples
#' df <- twoLevelGraph(letters[1:10], 100)
#' df <- generateFactorItems(df, list(f1=paste0('i',1:4),
#'                            f2=paste0('i',2:4)),
#'                       c(f1=0.9, f2=0.5))
#' head(df)
#' attr(df, "pathProp")
#' attr(df, "score")
#' attr(df, "worth")
#' @export
#' @importFrom stats rnorm sd plogis rbinom
generateFactorItems <- function(df, path, factorScalePrior=deprecated(), th=0.5, name, ..., scale=1, alpha=1) {
  if (lifecycle::is_present(factorScalePrior)) {
    deprecate_warn("1.5.0", "prepSingleFactorModel(factorScalePrior = )")
  }
  palist <- verifyIsData(df)
  items <- Reduce(union, path, c())
  numItems <- length(items)
  if (missing(name)) {
    num <- ncol(df)-1
    name <- paste0('i', num:(num+numItems-1))
  }
  assertNameUnused(df, name)
  validateFactorModel(items, path)

  numFactors <- length(path)
  itemsPerFactor <- sapply(path, length)
  factorItemPath <- matrix(c(rep(1:length(itemsPerFactor), itemsPerFactor),
                             unlist(lapply(path, function(x) match(x, items)))),
                           nrow=2, byrow=TRUE)
  prop <- matrix(0, ncol=numItems, nrow=numFactors+1)
  propTry <- 1L
  while(1) {
    for (fx in 1:numFactors) {
      sel <- factorItemPath[1,]==fx
      prop[fx+1, factorItemPath[2,sel]] <-
        2*(plogis(abs(rnorm(sum(factorItemPath[1,]==fx), sd=.6))) - 0.5)
    }
    if (all(colSums(prop) < .99)) break
    propTry <- propTry + 1L
    if (propTry > 10L) stop("factorScalePrior is too large")
  }
  prop[1,] <- 1 - colSums(prop)

  # random sign
  prop <- prop * matrix((2*rbinom(prod(dim(prop)), 1, .5)-1),
                        nrow(prop), ncol(prop))
  for (fx in 1:numFactors) {
    c <- which(prop[1+fx,]!=0)[1]
    prop[1+fx, c] <- abs(prop[1+fx, c])
  }

  theta <- matrix(0, nrow=length(palist), ncol=numItems)
  for (ix in 1:numItems) {
    theta[,ix] <- ssqrt(prop[1,ix]) * c(scale(rnorm(length(palist))))
  }

  ftheta <- matrix(0, nrow=length(palist), ncol=numFactors)
  for (fx in 1:numFactors) {
    ftheta[,fx] <- c(scale(rnorm(length(palist))))
  }
  for (px in 1:ncol(factorItemPath)) {
    fx <- factorItemPath[1,px]
    ix <- factorItemPath[2,px]
    theta[,ix] <- theta[,ix] + ftheta[,fx] * ssqrt(prop[1+fx,ix])
  }

  theta <- scale(theta)
  dimnames(theta) <- list(palist, name)
  df <- generateItem(df, theta, th, scale=scale, alpha=alpha)
  pathProp <- mapply(function(r,c,x) x[r,c],
                     1+factorItemPath[1,], factorItemPath[2,], MoreArgs = list(prop))
  attr(df, "pathProp") <- pathProp
  attr(df, "score") <- ftheta
  attr(df, "worth") <- theta
  df
}

#' Generate paired comparison data with random correlations between items
#'
#' @inheritParams generateItem
#' @param numItems how many items to create
#' @description
#'
#' If you need access to the correlation matrix used to generate the
#' absolute latent scores then you will need to generate them yourself.
#' This is not difficult. See how in the example.
#'
#' @template detail-response
#' @family item generators
#' @return
#' The given data.frame \code{df} with additional columns for each item.
#' In addition, you can obtain the correlation matrix used
#' to generate the latent worths from \code{attr(df, "cor")} and
#' and latent worths from \code{attr(df, "worth")}.
#'
#' @examples
#' library(mvtnorm)
#' df <- twoLevelGraph(letters[1:10], 100)
#' df <- generateCovItems(df, 3)
#'
#' # generateCovItems essentially does the same thing as:
#' numItems <- 3
#' palist <- letters[1:10]
#' trueCor <- cov2cor(rWishart(1, numItems, diag(numItems))[,,1])
#' theta <- rmvnorm(length(palist), sigma=trueCor)
#' dimnames(theta) <- list(palist, paste0('i', 3 + 1:numItems))
#' df <- generateItem(df, theta)
#' attr(df, "cor")
#'
#' @export
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats cov2cor rWishart
generateCovItems <- function(df, numItems, th=0.5, name, ..., scale=1, alpha=1) {
  if (numItems < 2) stop("numItems must be 2 or greater")
  palist <- verifyIsData(df)
  if (missing(name)) {
    num <- ncol(df)-1
    name <- paste0('i', num:(num+numItems-1))
  }
  assertNameUnused(df, name)
  trueCor <- cov2cor(rWishart(1, numItems, diag(numItems))[,,1])
  theta <- rmvnorm(length(palist), sigma=trueCor)
  dimnames(theta) <- list(palist, name)
  df <- generateItem(df, theta, th, scale=scale, alpha=alpha)
  attr(df, "cor") <- trueCor
  attr(df, "worth") <- theta
  df
}

#' Generate paired comparison data for one or more items given
#' absolute latent scores
#'
#' @template args-df
#' @param theta a vector or matrix of absolute latent scores. See details below.
#' @param th a vector of thresholds
#' @param name a vector of item names
#' @template args-dots-barrier
#' @param scale a vector of scaling constants
#' @param alpha a vector of item discriminations
#'
#' @description
#' To add a single item, \code{theta} should be a vector of latent
#' scores. To add multiple items at a time, \code{theta} should be a
#' matrix with one item in each column. Item names can be given as
#' the colnames of \code{theta}.
#'
#' The interpretation of \code{theta} depends on the context where the
#' data were generated. For example, in chess, \code{theta} represents
#' unobserved chess skill that is partially revealed by match
#' outcomes.
#'
#' The graph can be regarded as undirected, but data are generated
#' relative to the order of vertices within each row. Vertices do not
#' commute. For example, a \code{-1} for vertices \sQuote{a} and
#' \sQuote{b} is the same as \code{1} for vertices \sQuote{b} and
#' \sQuote{a}.
#'
#' @template detail-response
#' @family item generators
#' @return
#' The given data.frame \code{df} with additional columns for each item.
#' @examples
#' df <- roundRobinGraph(letters[1:5], 40)
#' df <- generateItem(df)
#' @export
generateItem <- function(df, theta, th=0.5, name, ..., scale=1, alpha=1) {
  if (length(list(...)) > 0) stop("Rejected are any values passed in the '...' argument")
  if (!missing(theta) && is.matrix(theta)) {
    if (length(colnames(theta))) {
      if (!missing(name)) {
        if (!all(name == colnames(theta))) {
          stop("Mismatch between name and colnames(theta)")
        }
      } else {
        name <- colnames(theta)
      }
    } else if (missing(name)) {
      num <- ncol(df)-1
      name <- paste0('i', num:(num+ncol(theta)-1))
    }
    scale <- c(matrix(scale, nrow=ncol(theta)))
    alpha <- c(matrix(alpha, nrow=ncol(theta)))
    for (ix in 1:ncol(theta)) {
      df <- generateItem(df, theta[,ix], th, name[ix],
                         scale=scale[ix], alpha=alpha[ix])
    }
    return(df)
  }

  palist <- verifyIsData(df)
  if (missing(theta)) {
    theta <- rnorm(length(palist))
    theta <- c(scale(theta))
    names(theta) <- palist
  }
  if (missing(name)) {
    name <- paste0('i', ncol(df)-1)
  }
  assertNameUnused(df, name)
  if (length(theta) != length(palist)) {
    stop(paste("length(theta)",length(theta),"!=", length(palist)))
  }
  missing <- is.na(match(palist, names(theta)))
  if (any(missing)) {
    stop(paste("No latent score for:",
      paste(palist[missing], collapse=' ')))
  }
  pick <- 1:length(th)
  pick <- c(-rev(pick),0,pick)
  df[[name]] <- NA
  for (rx in 1:nrow(df)) {
    p1 <- match(df[rx,'pa1'], palist)
    p2 <- match(df[rx,'pa2'], palist)
    prob <- cmp_probs(scale, alpha, theta[p1], theta[p2], th)
    df[rx,name] <- sample(pick, 1, prob=prob)
  }
  df
}

objectNamesToFactor <- function(name, df) {
  for (v in 1:2) {
    df[[ paste0('pa',v) ]] <- factor(df[[ paste0('pa',v) ]], name)
  }
  df
}

#' Create an edge list with round-robin connectivity
#'
#' @param name vector of vertex names
#' @param N number of comparisons
#' @return An undirected graph represented as a data frame with each row describing an edge.
#' @family graph generators
#' @examples
#' roundRobinGraph(letters[1:5], 10)
#' @export
roundRobinGraph <- function(name, N) {
  if (length(name) < 2) stop("Must provide at least 2 names")
  pmap <- pairMap(length(name))
  if (N < ncol(pmap)) {
    warning(paste("Sample size too small", N,
                  "to round robin connect all the vertices", ncol(pmap)))
  }
  df <- data.frame(pa1=rep(NA,N), pa2=NA)
  for (rx in 1:N) {
    pick <- pmap[,1 + (rx-1) %% ncol(pmap)]
    df[rx,'pa1'] <- name[pick[1]]
    df[rx,'pa2'] <- name[pick[2]]
  }
  df <- objectNamesToFactor(name, df)
  df
}

#' Create an edge list with a random two level connectivity
#'
#' @param shape1 beta distribution parameter for first edge
#' @param shape2 beta distribution parameter for second edge
#'
#' @description
#' Initially, edges are added from the first vertex to all the other
#' vertices. Thereafter, the first vertex is drawn from a Beta(shape1,
#' 1.0) distribution and the second vertex is drawn from a
#' Beta(shape2, 1.0) distribution. The idea is that the edges will
#' tend to connect a small subset of vertices from the top of the tree
#' to leaf vertices. These vertex connections are similar to the pairs
#' that you might observe in an elimination tournament. The selected
#' vertices are sorted so it doesn't matter whether \code{shape1 >
#' shape2} or \code{shape1 < shape2}.
#'
#' @inheritParams roundRobinGraph
#' @inherit roundRobinGraph return
#' @family graph generators
#' @examples
#' twoLevelGraph(letters[1:5], 20)
#' @export
twoLevelGraph <- function(name, N, shape1=0.8, shape2=0.5) {
  if (length(name) < 2) stop("Must provide at least 2 names")
  if (N < length(name)) {
    warning(paste("Sample size too small", N,
                  "to connect all the vertices", length(name)))
  }
  if (shape1 > 1 || shape1 < 0) stop("shape1 must be between 0 and 1")
  if (shape2 > 1 || shape2 < 0) stop("shape2 must be between 0 and 1")

  df <- data.frame(pa1=rep(NA,N), pa2=NA)

  # Generate all comparisons between 'a' and others to ensure
  # that there are no disconnected vertices.
  for (rx in 1:min(N,(length(name)-1))) {
    df[rx,'pa1'] <- name[1]
    ox <- 1 + (rx %% length(name))
    df[rx,'pa2'] <- name[ox]
  }
  if (N >= length(name)) {
    for (rx in length(name):nrow(df)) {
      pick <- 1+floor(sort(c(rbeta(1, shape1, 1),
                             rbeta(1, shape2, 1))) * length(name))
      if (pick[1] == pick[2]) {
        if (pick[1] > 1) pick[1] <- pick[1] - 1
        else pick[2] <- pick[2] + 1
      }
      df[rx,'pa1'] <- name[pick[1]]
      df[rx,'pa2'] <- name[pick[2]]
    }
  }
  df <- objectNamesToFactor(name, df)
  df
}

#' Filter graph to remove vertices that are not well connected
#'
#' @template args-df
#' @param minAny the minimum number of edges
#' @param minDifferent the minimum number of vertices
#'
#' @description
#'
#' Vertices not part of the largest connected component are excluded (Hopcroft & Tarjan, 1973).
#' Vertices that have fewer than \code{minAny} edges and are not
#' connected to \code{minDifferent} or more different vertices are
#' excluded. For example, vertex \sQuote{a} connected to vertices
#' \sQuote{b} and \sQuote{c} will be include so long as these vertices
#' are part of the largest connected component.
#'
#' @details
#' Given that \code{minDifferent} defaults to 2,
#' if activity \eqn{A} was compared to at least
#' two other activities, \eqn{B} and \eqn{C}, then \eqn{A} is retained.
#' The rationale is that,
#' although little may be learned about \eqn{A},
#' there may be a transitive relationship,
#' such as \eqn{B < A < C}, by which the model can infer that \eqn{B < C}.
#' Therefore, per-activity sample size is less of a concern
#' when the graph is densely connected.
#'
#' A young novice asked the wise master, "Why is 11 the default \code{minAny} instead of 10?"
#' The master answered, "Because 11 is a prime number."
#'
#' @return The same graph excluding some
#'   vertices.
#'
#' @importFrom igraph graph_from_edgelist components incident
#' @references
#' Hopcroft, J., & Tarjan, R. (1973). Algorithm 447: Efficient algorithms for graph
#' manipulation. \emph{Communications of the ACM, 16}(6), 372–378.
#' doi:10.1145/362248.362272
#' @examples
#' df <- filterGraph(phyActFlowPropensity[,c(paste0('pa',1:2),'predict')])
#' head(df)
#'
#' @export
filterGraph <- function(df, minAny=11L, minDifferent=2L) {
  verifyIsData(df)
  el <- as.matrix(df[,c(paste0('pa',1:2))])
  gr <- graph_from_edgelist(el, directed = FALSE)
  c1 <- components(gr)
  largest <- which(c1$csize == max(c1$csize))
  disconnected <- names(c1$membership[c1$membership != largest])
  good <- names(c1$membership[c1$membership == largest])
  weak <- c()
  for (g1 in good) {
    ie <- incident(gr, g1)
    vn <- attr(ie, 'vnames')
    if (length(vn) >= minAny || sum(!duplicated(vn)) >= minDifferent) next
    weak <- c(weak, g1)
  }
  good <- setdiff(good, weak)

  df <- df[df$pa1 %in% good & df$pa2 %in% good,]
  for (v in paste0('pa',1:2)) {
    # remove weak from factor levels
    df[[v]] <- factor(as.character(df[[v]]), good)
  }
  attr(df, 'disconnected') <- disconnected
  attr(df, 'weak') <- weak
  class(df) <- c("filteredGraph", class(df))
  df
}

#' @export
print.filteredGraph <- function(x, ...) {
  df <- x
  class(df) <- 'data.frame'
  print(df, ...)
  dis <- attr(x, 'disconnected')
  weak <- attr(x, 'weak')
  if (length(c(dis,weak))) {
    message("Some vertices were excluded,")
    if (length(dis)) {
      message("  not part of the largest connected component: ", deparse(dis))
    }
    if (length(weak)) {
      message("  too weakly connected: ", deparse(weak))
    }
  }
}

#' Turn a factor back into a vector of integers
#'
#' @param f a factor
#'
#' @description
#'
#' Factors store values as integers and use a 'levels' attribute to
#' map the integers to labels. This function removes the 'factor'
#' class and levels attribute, leaving the vector of integers.
#'
#' @examples
#' f <- factor(letters[1:3])
#' print(f)
#' print(unfactor(f))
#' @export
unfactor <- function(f) {
  if (!is.factor(f)) return(f)
  f <- unclass(f)
  levels(f) <- c()
  f
}
