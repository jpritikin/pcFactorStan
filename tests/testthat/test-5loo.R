library(testthat)
context("test-5loo")

skip_on_cran()
options(mc.cores=4)

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
  m1 <- pcStan("unidim_ll", dl, chains=4L, iter=600)

  loo1 <- toLoo(m1, cores=1)
  expect_error(outlierTable(df, loo1),
               "Data must be processed by prepData")
})

test_that("loo", {
  set.seed(1)
  palist <- letters[1:10]
  numItems <- 3

  df1 <- twoLevelGraph(palist, 300)
  trueCor <- cov2cor(rWishart(1, numItems, diag(numItems))[,,1])
  theta <- rmvnorm(length(palist), sigma=trueCor)
  dimnames(theta) <- list(palist, paste0('i', 1:numItems))
  df1 <- generateItem(df1, theta, th=rep(0.5, 4))

  df2 <- twoLevelGraph(palist, 30)
  df2 <- generateItem(df2, -theta, th=rep(0.5, 4))
  for (ix in 1:numItems) {
    mask <- as.logical(rbinom(nrow(df2), 1, 1/3))
    df2[mask,paste0('i', ix)] <- NA
  }

  df <- rbind(df1, df2)

  for (k in paste0('pa',1:2)) df[[k]] <- factor(df[[k]], levels=palist)
  df <- filterGraph(df)
  df <- normalizeData(df)
  dl <- prepCleanData(df)
  dl$scale <- 1.5

  # this *sometimes* generates warnings about divergent transitions
  m1 <- pcStan("covariance_ll", dl, chains=4L, iter=1000L)

  loo1 <- toLoo(m1, cores=1)
  ot <- outlierTable(dl, loo1)
  bad <- df[df$pa1==ot[1,'pa1'] & df$pa2==ot[1,'pa2'],]
  item <- as.character(ot[1,'item'])
  expect_true(all(sign(bad[as.numeric(rownames(bad)) <= 300, item])
                  == sign(bad[1, item])))
  expect_true(all(sign(bad[as.numeric(rownames(bad)) > 300, item])
                  == sign(bad[nrow(bad), item]), na.rm=TRUE))
  expect_true(sign(bad[1, item]) !=
                sign(bad[nrow(bad), item]))

  m2 <- pcStan("factor_ll", dl, chains=4L, iter=1000L)
  loo2 <- toLoo(m2, cores=1)
  ot <- outlierTable(dl, loo2)
  bad <- df[df$pa1==ot[1,'pa1'] & df$pa2==ot[1,'pa2'],]
  item <- as.character(ot[1,'item'])
  expect_true(all(sign(bad[as.numeric(rownames(bad)) <= 300, item])
                  == sign(bad[1, item])))
  expect_true(all(sign(bad[as.numeric(rownames(bad)) > 300, item])
                  == sign(bad[nrow(bad), item]), na.rm=TRUE))
  expect_true(sign(bad[1, item]) !=
                sign(bad[nrow(bad), item]))
})
