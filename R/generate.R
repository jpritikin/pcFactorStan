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

verifyIsData <- function(df) {
  if (!is.data.frame(df)) stop("df is not a data.frame")
  if (all(match(paste0('pa',1:2), colnames(df)) != c(1,2))) {
    stop("Expected pa1, pa2 as the first two columns")
  }
  unique(c(df$pa1,df$pa2))
}

assertNameUnused <- function(df, name) {
  collision <- !is.na(match(name, colnames(df)))
  if (any(collision)) {
    stop(paste("Colname", paste(name[collision], collapse=' '),
               "is already taken. Please provide name="))
  }
}

#' @export
generateFactorItems <- function(df, prop, th=-0.5, scale=1, name) {
  palist <- verifyIsData(df)
  if (length(prop) == 1) {
    if (prop < 3) stop(paste0("At least 3 indicators are required (", prop," requested)"))
    prop <- rbeta(prop, 4, 3)
  }
  if (length(prop) < 3) stop(paste0("At least 3 indicators are required (", length(prop)," given)"))
  if (missing(name)) {
    num <- ncol(df)-1
    name <- paste0('i', num:(num+length(prop)-1))
  }
  assertNameUnused(df, name)
  factorScore <- rnorm(length(palist))
  factorScore <- factorScore / sd(factorScore)
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

#' @export
#' @importFrom mvtnorm rmvnorm
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
