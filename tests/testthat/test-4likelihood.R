context("test-4likelihood")

skip_on_cran()

RNGversion("3.5")

test_that("unidim", {
  # As of rstan 2.19, cores=0 suppresses the warnings about chain convergence.

  dl <- prepData(phyActFlowPropensity[,c(1,2,3)])

  expect_error(locateModel("unidim+adapt", data=dl),
               "You must choose a varCorrection.")
  expect_error(locateModel("unidim+ll", data=dl),
               "You must choose a scale.")

  dl$varCorrection <- 2.0
  m1 <- stan_model(locateModel("unidim+adapt", data=dl))
  f1 <- sampling(m1, dl, chains=1, cores=0, iter=1, seed=1,warmup=0)
  expect_equal(get_logposterior(f1)[[1]], -3846.102, tolerance=1e-2, scale=1)

  dl$scale <- 1.0
  m2 <- stan_model(locateModel("unidim+ll", data=dl))
  f2 <- sampling(m2, dl, chains=1, cores=0, iter=1, seed=1,warmup=0)
  expect_equal(get_logposterior(f2)[[1]], -2921.243, tolerance=1e-2, scale=1)
  #cat(deparse(round(fivenum(extract(f2)$log_lik[1,]), 3)))
  expect_equal(fivenum(extract(f2)$log_lik[1,]),
               c(-14.382, -3.49, -2.01, -1.338, -0.032), tolerance=1e-2, scale=1)
})

test_that("covariance", {
  dl <- prepData(phyActFlowPropensity)
  dl$scale <- 1.5
  m2 <- stan_model(locateModel("covariance+ll", data=dl))
  f2 <- sampling(m2, dl, chains=1, cores=0, iter=1, seed=1,warmup=0)
  expect_equal(get_logposterior(f2)[[1]], -114592.3, tolerance=1e-2, scale=1)
  #cat(deparse(round(fivenum(extract(f2)$log_lik[1,]), 3)))
  expect_equal(fivenum(extract(f2)$log_lik[1,]),
               c(-191.855, -6.121, -2.951, -1.731, 0), tolerance=1e-2, scale=1)
})

test_that("factor", {
  dl <- prepData(phyActFlowPropensity)
  dl$scale <- 1.5
  m2 <- stan_model(locateModel("factor+ll", data=dl))
  f2 <- sampling(m2, dl, chains=1, cores=0, iter=1, seed=1,warmup=0)
  expect_equal(get_logposterior(f2)[[1]], -114116.3, tolerance=1e-2, scale=1)
  #cat(deparse(round(fivenum(extract(f2)$log_lik[1,]), 3)))
  expect_equal(fivenum(extract(f2)$log_lik[1,]),
               c(-60.806, -10.001, -4.726, -1.691, 0), tolerance=1e-2, scale=1)
})

test_that("mixed thresholds", {
  set.seed(1)
  palist <- letters[1:10]
  df <- twoLevelGraph(palist, 300)
  numItems <- 5
  trueCor <- cov2cor(rWishart(1, numItems, diag(numItems))[,,1])
  theta <- rmvnorm(length(palist), sigma=trueCor)
  dimnames(theta) <- list(palist, paste0('i', 1:numItems))
  for (ix in 1:numItems) {
    df <- generateItem(df, theta[,ix], th=rep(-0.5, ix))
  }

  dl <- prepData(df)
  dl$scale <- 1.5
  m2 <- stan_model(locateModel("covariance+ll", data=dl))
  f2 <- sampling(m2, dl, chains=1, cores=0, iter=1, seed=1,warmup=0)
  expect_equal(get_logposterior(f2)[[1]], -4611.471, tolerance=1e-2, scale=1)
  #cat(deparse(round(fivenum(extract(f2)$log_lik[1,]), 3)))
  expect_equal(fivenum(extract(f2)$log_lik[1,]),
               c(-20.86, -3.784, -2.167, -0.838, -0.002), tolerance=1e-2, scale=1)
})