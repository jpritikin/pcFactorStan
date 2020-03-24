library(testthat)
library(pcFactorStan)
context("test-5loo")

skip_on_cran()
options(mc.cores=2)

suppressWarnings(RNGversion("3.6"))
library(rstan)  # for get_logposterior
library(mvtnorm)  # rmvnorm

test_that("loo univariate", {
  set.seed(1)
  palist <- letters[1:10]

  df <- twoLevelGraph(palist, 300)
  theta <- rnorm(length(palist))
  names(theta) <- palist
  df <- generateItem(df, theta, th=rep(0.5, 4))

  df <- filterGraph(df)
  df <- normalizeData(df)
  dl <- prepCleanData(df)
  dl$scale <- 1.5
  m1 <- suppressWarnings(pcStan("unidim_ll", dl, chains=2L, iter=300))

  loo1 <- toLoo(m1, cores=1)
  expect_error(outlierTable(df, loo1),
               "Data must be processed by prepData")
})

test_that("loo", {
  library(testthat)
  library(pcFactorStan)
  options(mc.cores=2)
  suppressWarnings(RNGversion("3.5"))
  library(rstan)  # for get_logposterior
  library(mvtnorm)  # rmvnorm

  set.seed(1)
  palist <- letters[1:10]
  numItems <- 3

  boundary <- 300
  df1 <- twoLevelGraph(palist, boundary)
  trueCor <- cov2cor(rWishart(1, numItems, diag(numItems))[,,1])
  theta <- rmvnorm(length(palist), sigma=trueCor)
  dimnames(theta) <- list(palist, paste0('i', 1:numItems))
  df1 <- generateItem(df1, theta, th=rep(0.5, 4), alpha=1.25)

  df2 <- twoLevelGraph(palist, 30)
  df2 <- generateItem(df2, -theta, th=rep(0.5, 4), alpha=1.25)

  df <- rbind(df1, df2)

  for (k in paste0('pa',1:2)) df[[k]] <- factor(df[[k]], levels=palist)
  df <- filterGraph(df)
  df <- normalizeData(df)
  dl <- prepCleanData(df)
  dl$scale <- rep(1.5, dl$NITEMS)

  # this *sometimes* generates warnings about divergent transitions
  m1 <- suppressWarnings(pcStan("correlation_ll", dl, chains=2L, iter=400L))

  loo1 <- suppressWarnings(toLoo(m1, cores=1))
  ot <- outlierTable(dl, loo1, .4)
  expect_true(nrow(ot) > 0)

  bad <- df[df$pa1==ot[1,'pa1'] & df$pa2==ot[1,'pa2'],]
  # no idea how to write tests for this TODO
})
