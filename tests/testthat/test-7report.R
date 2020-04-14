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
  expect_equal(nrow(rc), 30)
  expect_true(all(diff(subset(rc, response=='a' & sample == 1)$prob) < 0))
  expect_true(all(diff(subset(rc, response=='e' & sample == 2)$prob) > 0))

  squash <- phyActFlowPropensity[,c(1,2,3)]
  squash <- squash[!is.na(squash$skill),]
  squash[squash$skill == 2, 'skill'] <- 1
  squash[squash$skill == -2, 'skill'] <- -1
  dl4 <- prepData(squash)
  expect_error(responseCurve(dl4, f1, letters[1:3]),
               "1 thresholds across all items but fit has 2 thresholds")
  expect_error(responseCurve(letters, f1),
               "dl must be a list of data")
  expect_error(responseCurve(dl1, letters[1:5]),
               "fit must be a stanfit object")
  expect_error(responseCurve(dl1, f1, letters[1:5], "zorg"),
               "Item not found: zorg")

  dl2 <- prepData(phyActFlowPropensity[,c(1:5)])
  dl2$scale <- rnorm(dl2$NITEMS, .8, .1)
  m2 <- findModel("correlation")
  f2 <- sampling(m2, dl2, chains=1, cores=0, iter=50, refresh=0)
  rc <- responseCurve(dl2, f2, letters[1:5], 'predict', samples=2, by=1)
  expect_equal(nrow(rc), 30)
  expect_true(all(diff(subset(rc, response=='a' & sample == 1)$prob) < 0))
  expect_true(all(diff(subset(rc, response=='e' & sample == 2)$prob) > 0))

  dl3 <- prepCleanData(filterGraph(phyActFlowPropensity[,c(1:5)],
                                   minAny = 15, minDifferent = 3))
  expect_error(responseCurve(dl3, f2, letters[1:5]),
               "dl has 49 objects but fit has 61 objects")
  expect_error(responseCurve(dl2, f2, letters[1:4]),
               "different number of responseNames")
  expect_error(responseCurve(dl2, f1, letters[1:5]),
               "dl has 3 items but fit has 1 items")
  expect_error(responseCurve(dl1, f2, letters[1:5]),
               "dl has 1 items but fit has 3 items")

  dl2 <- prepSingleFactorModel(dl2)
  f3 <- sampling(findModel("factor1"), dl2, chains=1,
                 cores=0, iter=50, refresh=0)
  rc <- responseCurve(dl2, f3, letters[1:5], 'predict', samples=2, by=1)
  expect_equal(nrow(rc), 30)
  expect_true(all(diff(subset(rc, response=='a' & sample == 1)$prob) < 0))
  expect_true(all(diff(subset(rc, response=='e' & sample == 2)$prob) > 0))
})

test_that("parInterval+parDistributionFor", {
  set.seed(1)
  dl1 <- prepData(phyActFlowPropensity[,c(1,2,3)])
  dl1$scale <- 1.0
  expect_error(parInterval(dl1, "alpha", nameVec="alpha"),
               "must be a stanfit object")
  m1 <- findModel("unidim")
  f1 <- sampling(m1, dl1, chains=1, cores=0, iter=100, refresh=0)
  label <- "discrimination"
  expect_error(parInterval(f1, "alpha", nameVec=paste0("alpha",1:3)),
               "pars and nameVec must be the same length.")
  pi <- parInterval(f1, "alpha", nameVec="alpha")
  expect_equal(colnames(pi)[4], "alpha")
  pi <- parInterval(f1, "alpha", "alpha", label)
  expect_equal(nrow(pi), 1)
  expect_equal(ncol(pi), 4)
  expect_equal(rownames(pi), "alpha")
  expect_equal(colnames(pi)[1:3], c('L','M','U'))
  expect_equal(colnames(pi)[4], label)
  expect_error(parDistributionFor(dl1, pi),
               "must be a stanfit object")
  pd <- parDistributionFor(f1, pi)
  expect_equivalent(c(table(pd$value < pi$M)), c(25,25))
  expect_equal(nrow(pd), 50)
  pd <- parDistributionFor(f1, pi, samples = 5)
  expect_equal(nrow(pd), 5)
  expect_error(parDistributionCustom(f1, pars = "alpha", nameVec = paste0("alpha",1:2)),
               "pars and nameVec must be the same length.")

  label <- "activity"
  pi <- parInterval(f1, "theta", dl1$nameInfo$pa)
  expect_equal(colnames(pi)[4], "theta")
  pi <- parInterval(f1, "theta", dl1$nameInfo$pa, label)
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
