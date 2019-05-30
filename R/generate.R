softmax <- function(y) exp(y) / sum(exp(y))

cmp_probs <- function(scale, pa1, pa2, thRaw) {
  th <- cumsum(thRaw)
  diff <- scale * (pa1 - pa2);
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

#generateCovItems <- function(df)
#generateFactorItems <- function(df)

verifyIsData <- function(df) {
  if (!is.data.frame(df)) stop("df is not a data.frame")
  if (all(match(paste0('pa',1:2), colnames(df)) != c(1,2))) {
    stop("Expected pa1, pa2 as the first two columns")
  }
}

#' @export
generateItem <- function(df, theta, th=-0.5, scale=1, name) {
  verifyIsData(df)
  palist <- unique(c(df$pa1,df$pa2))
  if (missing(theta)) {
    theta <- rnorm(length(palist))
    theta <- c(scale(theta))
    names(theta) <- palist
  }
  if (missing(name)) {
    name <- paste0('i', ncol(df)-1)
  }
  if (name %in% colnames(df)) {
    stop(paste("Colname",name,"is already taken. Please provide name="))
  }
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

genData <- function(df, theta, th=-0.5, scale=1) {
  numIndicators <- ncol(theta)
  for (ix in 1:numIndicators) df[[paste0('i',ix)]] <- NA
  palist <- rownames(theta)
  for (rx in 1:nrow(df)) {
    p1 <- match(df[rx,'pa1'], palist)
    p2 <- match(df[rx,'pa2'], palist)
    for (ix in 1:numIndicators) {
      prob <- cmp_probs(scale, theta[p1,ix], theta[p2,ix], th)
      df[rx,paste0('i',ix)] <- sample(c(-2,-1,0,1,2), 1, prob=prob)
    }
  }
  df
}

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

#' @export
twoLevelGraph <- function(name, N, shape1, shape2) {
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
