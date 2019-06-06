context("test-4likelihood")

skip_on_cran()

RNGversion("3.5")
library(rstan)  # for get_logposterior

test_that("unidim", {
  # As of rstan 2.19, cores=0 suppresses the warnings about chain convergence.

  expect_error(pcStan('unidim', data=phyActFlowPropensity[,1:3]),
               "Data must be processed by prepData")

  dl <- prepData(phyActFlowPropensity[,c(1,2,3)])

  dl$varCorrection <- 2.0
  m1 <- findModel("unidim_adapt")
  f1 <- sampling(m1, dl, chains=1, cores=0, iter=1, seed=1,warmup=0)
  expect_equal(get_logposterior(f1)[[1]], -53935.71, tolerance=1e-2, scale=1)

  dl$scale <- 1.0
  m2 <- findModel("unidim_ll")
  f2 <- sampling(m2, dl, chains=1, cores=0, iter=1, seed=1,warmup=0)
  expect_equal(get_logposterior(f2)[[1]], -3174.317, tolerance=1e-2, scale=1)
  #cat(deparse(round(fivenum(extract(f2)$log_lik[1,]), 3)))
  expect_equal(fivenum(extract(f2)$log_lik[1,]),
               c(-14.4, -4.84, -2.343, -0.792, -0.26), tolerance=1e-2, scale=1)
})

test_that("covariance", {
  dl <- prepData(phyActFlowPropensity)
  dl$scale <- 1.5
  m2 <- findModel("covariance_ll")
  f2 <- sampling(m2, dl, chains=1, cores=0, iter=1, seed=1,warmup=0)
  expect_equal(get_logposterior(f2)[[1]], -89848.91, tolerance=1e-2, scale=1)
  #cat(deparse(round(fivenum(extract(f2)$log_lik[1,]), 3)))
  expect_equal(fivenum(extract(f2)$log_lik[1,]),
               c(-103.731, -5.783, -2.97, -1.486, 0), tolerance=1e-2, scale=1)
})

test_that("factor", {
  dl <- prepData(phyActFlowPropensity)
  dl$scale <- 1.5
  m2 <- findModel("factor_ll")
  f2 <- sampling(m2, dl, chains=1, cores=0, iter=1, seed=1,warmup=0)
  expect_equal(get_logposterior(f2)[[1]], -100692.33, tolerance=1e-2, scale=1)
  #cat(deparse(round(fivenum(extract(f2)$log_lik[1,]), 3)))
  expect_equal(fivenum(extract(f2)$log_lik[1,]),
               c(-48.981, -8.852, -4.65, -1.873, 0), tolerance=1e-2, scale=1)
})

test_that("mixed thresholds", {
  set.seed(1)
  palist <- letters[1:10]
  df <- twoLevelGraph(palist, 300)
  for (k in paste0('pa',1:2)) df[[k]] <- factor(df[[k]], levels=palist)
  numItems <- 5
  trueCor <- cov2cor(rWishart(1, numItems, diag(numItems))[,,1])
  theta <- rmvnorm(length(palist), sigma=trueCor)
  dimnames(theta) <- list(palist, paste0('i', 1:numItems))
  for (ix in 1:numItems) {
    df <- generateItem(df, theta[,ix], th=rep(-0.5, ix))
  }

  df <- filterGraph(df)
  dl <- prepCleanData(df)
  dl$scale <- 1.5
  m2 <- findModel("covariance_ll")
  f2 <- sampling(m2, dl, chains=1, cores=0, iter=1, seed=1,warmup=0)
  expect_equal(get_logposterior(f2)[[1]], -5633.51, tolerance=1e-2, scale=1)
  #cat(deparse(round(fivenum(extract(f2)$log_lik[1,]), 3)))
  expect_equal(fivenum(extract(f2)$log_lik[1,]),
               c(-41.785, -4.048, -1.776, -0.457, 0), tolerance=1e-2, scale=1)

  df <- normalizeData(df, .palist=sample(palist, 10))
  dl <- prepCleanData(df)
  dl$scale <- 1.5
  f3 <- sampling(m2, dl, chains=1, cores=0, iter=1, seed=1,warmup=0)
  expect_equal(get_logposterior(f3)[[1]],
               get_logposterior(f2)[[1]], tolerance=1e-2, scale=1)
  expect_equal(fivenum(extract(f2)$log_lik[1,]),
               fivenum(extract(f3)$log_lik[1,]), tolerance=1e-2, scale=1)
})

test_that("calibrateItems", {
  pafp <- phyActFlowPropensity[,1:5]
  result <- calibrateItems(pafp, iter=500L)
  expect_equal(nrow(result), 3)
  expect_true(all(result[,'n_eff'] > 400))
  expect_true(all(result[,'Rhat'] < 1.015))
  expect_equal(result[,'scale'], c(.64, .61, .50),
               tolerance=1e-2, scale=1)
})
