library(testthat)
library(pcFactorStan)

context("test-4likelihood")

skip_on_cran()
options(mc.cores=4)

suppressWarnings(RNGversion("3.5"))
library(rstan)  # for get_logposterior

test_that("unidim", {
  # As of rstan 2.19, cores=0 suppresses the warnings about chain convergence.

  expect_error(pcStan('unidim', data=phyActFlowPropensity[,1:3]),
               "Data must be processed by prepData")

  expect_error(pcStan('unidim', data=matrix(0, 3, 3)),
               "Is data an object returned by prepData")

  dl <- prepData(phyActFlowPropensity[,c(1,2,3)])

  dl$varCorrection <- 2.0
  m1 <- findModel("unidim_adapt")
  f1 <- sampling(m1, dl, chains=1, cores=0, iter=1, seed=1,warmup=0, refresh=0)
  expect_equal(get_logposterior(f1)[[1]], -3060.868, tolerance=1e-2, scale=1)

  dl$scale <- 1.0
  m2 <- findModel("unidim_ll")
  f2 <- sampling(m2, dl, chains=1, cores=0, iter=1, seed=1,warmup=0, refresh=0)
  expect_equal(get_logposterior(f2)[[1]], -6324.085, tolerance=1e-2, scale=1)
  #cat(deparse(round(fivenum(extract(f2)$log_lik[1,]), 3)))
  expect_equal(fivenum(extract(f2)$log_lik[1,]),
               c(-31.897, -10.49, -4.678, -0.609, -0.023), tolerance=1e-2, scale=1)
})

test_that("correlation", {
  set.seed(1)
  dl <- prepData(phyActFlowPropensity)
  dl$scale <- rnorm(dl$NITEMS, sd=.2)
  m2 <- findModel("correlation_ll")
  f2 <- sampling(m2, dl, chains=1, cores=0, iter=1, seed=1,warmup=0, refresh=0)
  expect_equal(get_logposterior(f2)[[1]], -59564.49, tolerance=1e-2, scale=1)
  #cat(deparse(round(fivenum(extract(f2)$log_lik[1,]), 3)))
  expect_equal(fivenum(extract(f2)$log_lik[1,]),
               c(-23.025, -4.753, -1.724, -1.386, 0), tolerance=1e-2, scale=1)
})

test_that("factor", {
  set.seed(1)
  dl <- prepData(phyActFlowPropensity)
  dl$scale <- rep(1.5, dl$NITEMS)
  dl$alpha <- rnorm(dl$NITEMS, 1, .1)
  dl <- prepSingleFactorModel(dl, .2)
  m2 <- findModel("factor1_ll")
  f2 <- sampling(m2, dl, chains=1, cores=0, iter=1, seed=1,warmup=0, refresh=0)
  expect_equal(get_logposterior(f2)[[1]], -108935.665, tolerance=1e-2, scale=1)
  #cat(deparse(round(fivenum(extract(f2)$log_lik[1,]), 3)))
  expect_equal(fivenum(extract(f2)$log_lik[1,]),
               c(-54.219, -9.448, -4.866, -2.03, 0), tolerance=1e-2, scale=1)
})

test_that("mixed thresholds", {
  library(mvtnorm)
  set.seed(1)
  palist <- letters[1:10]
  df <- twoLevelGraph(palist, 300)
  for (k in paste0('pa',1:2)) df[[k]] <- factor(df[[k]], levels=palist)
  numItems <- 5
  trueCor <- cov2cor(rWishart(1, numItems, diag(numItems))[,,1])
  theta <- rmvnorm(length(palist), sigma=trueCor)
  dimnames(theta) <- list(palist, paste0('i', 1:numItems))
  for (ix in 1:numItems) {
    df <- generateItem(df, theta[,ix], th=rep(0.5, ix))#
  }

  df <- filterGraph(df)
  dl <- prepCleanData(df)
  scaleSave <- rnorm(numItems, .9, .2)
  dl$scale <- scaleSave
  m2 <- findModel("correlation_ll")
  f2 <- sampling(m2, dl, chains=1, cores=0, iter=1, seed=1,warmup=0, refresh=0)
  expect_equal(get_logposterior(f2)[[1]], -10701.834, tolerance=1e-2, scale=1)
  #cat(deparse(round(fivenum(extract(f2)$log_lik[1,]), 3)))
  expect_equal(fivenum(extract(f2)$log_lik[1,]),
               c(-73.863, -10.18, -2.316, -1.524, -0.004), tolerance=1e-2, scale=1)

  df <- normalizeData(df, .palist=sample(palist, 10))
  dl <- prepCleanData(df)
  dl$scale <- scaleSave
  f3 <- sampling(m2, dl, chains=1, cores=0, iter=1, seed=1,warmup=0, refresh=0)
  expect_equal(get_logposterior(f3)[[1]],
               get_logposterior(f2)[[1]], tolerance=1e-2, scale=1)
  expect_equal(fivenum(extract(f2)$log_lik[1,]),
               fivenum(extract(f3)$log_lik[1,]), tolerance=1e-2, scale=1)
})

test_that("calibrateItems", {
  pafp <- phyActFlowPropensity[,c(1,2,6:8)]
  result <- calibrateItems(pafp, iter=1000L, chains=2, seed=1)
  expect_equal(nrow(result), 3)
  expect_true(all(result[1:2,'n_eff'] > 200))
  expect_true(all(result[1:2,'Rhat'] < 1.015))
  expect_true(result[3,'low_bfmi'] > 0 ||
                result[3,'divergent'] > 0)
  expect_true(result[3,'n_eff'] < 200)
  expect_equal(result[1:2,'scale'], c(.469, .535),
               tolerance=.01, scale=1)
})
