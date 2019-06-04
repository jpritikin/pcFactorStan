softmax <- function(y) exp(y) / sum(exp(y))

cmp_probs <- function(scale, pa1, pa2, thRaw) {
  th <- cumsum(thRaw)
  diff <- scale * (pa1 - pa2)
  unsummed <- c(0, c(diff - rev(th)), c(diff + th), use.names = FALSE)
  softmax(cumsum(unsummed))
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
    stop("Expected pa1, pa2 as the first two columns")
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

#' Generate pairwise comparison data with a common factor that
#' accounts for some proportion of the variance
#'
#' @param prop the number of items or a vector of proportions of variance
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
#' @template ref-masters1982
#' @template ref-silver2018
#' @family item generators
#' @examples
#' df <- twoLevelGraph(letters[1:10], 100)
#' df <- generateFactorItems(df, 3)
#' @export
#' @importFrom stats rnorm sd rbinom rbeta
generateFactorItems <- function(df, prop, th=-0.5, scale=1, name) {
  palist <- verifyIsData(df)
  if (length(prop) == 1) {
    if (prop < 3) stop(paste0("At least 3 indicators are required (", prop," requested)"))
    prop <- rbeta(prop, 4, 3)
  }
  if (length(prop) < 3) stop(paste0("At least 3 indicators are required (", length(prop)," given)"))
  if (any(prop < 0 | prop >= 1)) stop("Proportions must be between 0 and 1")
  if (missing(name)) {
    num <- ncol(df)-1
    name <- paste0('i', num:(num+length(prop)-1))
  }
  assertNameUnused(df, name)
  factorScore <- rnorm(length(palist))
  factorScore <- c(scale(factorScore))
  thetaF <- factorScore %*% t(sqrt(prop/(1-prop)))
  thetaN <- matrix(rnorm(length(palist)*length(prop)),
                   length(palist), length(prop))
  sdN <- apply(thetaN, 2, sd)
  for (cx in 1:ncol(thetaN)) {
    thetaN[,cx] <- thetaN[,cx] / sdN[cx]
  }

  # to double check
  ## varF <- apply(thetaF, 2, var)
  ## varN <- apply(thetaN, 2, var)
  ## genProp <- varF / (varF + varN)

  for (cx in 2:ncol(thetaF)) {
    if (rbinom(1, 1, .5)) thetaF[,cx] <- -thetaF[,cx]
  }
  theta <- thetaF + thetaN
  theta <- scale(theta)
  dimnames(theta) <- list(palist, name)
  generateItem(df, theta, th, scale)
}

#' Generate pairwise comparison data with random correlations between items
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
#' @template ref-masters1982
#' @family item generators
#' @examples
#' df <- twoLevelGraph(letters[1:10], 100)
#' df <- generateCovItems(df, 3)
#'
#' # generateCovItems essentially does the same thing as:
#' numItems <- 3
#' palist <- unique(c(df$pa1,df$pa2))
#' trueCor <- cov2cor(rWishart(1, numItems, diag(numItems))[,,1])
#' theta <- rmvnorm(length(palist), sigma=trueCor)
#' dimnames(theta) <- list(palist, paste0('i', 3 + 1:numItems))
#' df <- generateItem(df, theta)
#' 
#' @export
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats cov2cor rWishart
generateCovItems <- function(df, numItems, th=-0.5, scale=1, name) {
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
  generateItem(df, theta, th, scale)
}

#' Generate pairwise comparison data for one or more items given
#' absolute latent scores
#' 
#' @template args-df
#' @param theta a vector or matrix of absolute latent scores. See details below.
#' @param th a vector of thresholds
#' @param scale the scaling constant
#' @param name a vector of item names
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
#' relative to the order of vertices in the row. Vertices do not
#' commute. For example, a \code{-1} for vertices \sQuote{a} and
#' \sQuote{b} is the same as \code{1} for vertices \sQuote{b} and
#' \sQuote{a}.
#'
#' @template detail-response
#' @template ref-masters1982
#' @family item generators
#' @examples
#' df <- roundRobinGraph(letters[1:5], 40)
#' df <- generateItem(df)
#' @export
generateItem <- function(df, theta, th=-0.5, scale=1, name) {
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
    for (ix in 1:ncol(theta)) {
      df <- generateItem(df, theta[,ix], th, scale, name[ix])
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
    prob <- cmp_probs(scale, theta[p1], theta[p2], th)
    df[rx,name] <- sample(pick, 1, prob=prob)
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
#' shape2} or \code{shape1 < shape2} as long as \code{shape1 !=
#' shape2}.
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
  for (rx in 1:(length(name)-1)) {
    df[rx,'pa1'] <- name[1]
    ox <- 1 + (rx %% length(name))
    df[rx,'pa2'] <- name[ox]
  }
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
#' Vertices not part of the largest connected component are excluded.
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
