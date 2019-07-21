library(testthat)
library(pcFactorStan)
library(rstan)

context("test-7report")

test_that("responseCurve", {
  dl1 <- prepData(phyActFlowPropensity[,c(1,2,3)])
  dl1$scale <- 1.0
  m1 <- findModel("unidim")
  f1 <- sampling(m1, dl1, chains=1, cores=0, iter=100, refresh=0)
  rc <- responseCurve(dl1, f1, letters[1:5], samples=2, by=1)
  expect_equal(nrow(rc), 130)
  expect_true(all(diff(subset(rc, response=='a' & sample == 1)$prob) < 0))
  expect_true(all(diff(subset(rc, response=='e' & sample == 2)$prob) > 0))

  dl2 <- prepData(phyActFlowPropensity[,c(1:5)])
  dl2$scale <- rnorm(dl2$NITEMS, .8, .1)
  m2 <- findModel("correlation")
  f2 <- sampling(m2, dl2, chains=1, cores=0, iter=50, refresh=0)
  rc <- responseCurve(dl2, f2, letters[1:5], 'predict', samples=2, by=1)
  expect_equal(nrow(rc), 130)
  expect_true(all(diff(subset(rc, response=='a' & sample == 1)$prob) < 0))
  expect_true(all(diff(subset(rc, response=='e' & sample == 2)$prob) > 0))

  expect_error(responseCurve(dl2, f1, letters[1:5]),
               "dl has 3 items but fit has 1 items")
  expect_error(responseCurve(dl1, f2, letters[1:5]),
               "dl has 1 items but fit has 3 items")

  dl2$alpha <- rnorm(dl2$NITEMS, .8, .1)
  f3 <- sampling(findModel("factor"), dl2, chains=1,
                 cores=0, iter=50, refresh=0)
  rc <- responseCurve(dl2, f3, letters[1:5], 'predict', samples=2, by=1)
  expect_equal(nrow(rc), 130)
  expect_true(all(diff(subset(rc, response=='a' & sample == 1)$prob) < 0))
  expect_true(all(diff(subset(rc, response=='e' & sample == 2)$prob) > 0))
})

test_that("parInterval+parDistributionFor", {
  set.seed(1)
  dl1 <- prepData(phyActFlowPropensity[,c(1,2,3)])
  dl1$scale <- 1.0
  m1 <- findModel("unidim")
  f1 <- sampling(m1, dl1, chains=1, cores=0, iter=100, refresh=0)
  label <- "discrimination"
  pi <- parInterval(f1, "alpha", label, "alpha")
  expect_equal(nrow(pi), 1)
  expect_equal(ncol(pi), 4)
  expect_equal(rownames(pi), "alpha")
  expect_equal(colnames(pi)[1:3], c('L','M','U'))
  expect_equal(colnames(pi)[4], label)
  pd <- parDistributionFor(f1, pi)
  expect_equivalent(c(table(pd$value < pi$M)), c(25,25))

  label <- "activity"
  pi <- parInterval(f1, "theta", label, dl1$nameInfo$pa)
  expect_equal(nrow(pi), 61)
  expect_equal(ncol(pi), 4)
  expect_equal(rownames(pi), paste0('theta[',1:61,']'))
  expect_equal(colnames(pi)[1:3], c('L','M','U'))
  expect_equal(colnames(pi)[4], label)
  pd <- parDistributionFor(f1, pi)
  tbl <- table(subset(pd, activity=='mountain biking')$value <
    subset(pi, activity=='mountain biking')$M)
  expect_equivalent(c(tbl), c(25,25))
})
