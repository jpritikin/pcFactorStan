library(testthat)
library(pcFactorStan)

test_that("unidim", {
  skip("Must be checked manually")

  df <- twoLevelGraph(letters[1:10], 300)
  df <- generateItem(df)
  dl <- prepData(df)
  dl$scale <- 1
  fit <- pcStan("unidim", dl, iter=100, chains=1)
  itemModelExplorer(dl, fit, "i1")
})

test_that("correlation", {
  skip("Must be checked manually")

  numItems <- 3
  df <- twoLevelGraph(letters[1:10], 300)
  df <- generateCovItems(df, numItems)
  dl <- prepData(df)
  dl$scale <- rep(1, numItems)
  fit <- pcStan("correlation", dl, iter=100, chains=1)
  itemModelExplorer(dl, fit, "i3")
})

test_that("factor", {
  skip("Must be checked manually")

  numItems <- 4
  df <- twoLevelGraph(letters[1:10], 300)
  df <- generateSingleFactorItems(df, numItems)
  dl <- prepData(df)
  dl$scale <- rep(1, numItems)
  dl$alpha <- rnorm(numItems, .8, .15)
  dl <- prepSingleFactorModel(dl, 0.2)
  fit <- suppressWarnings(pcStan("factor1", dl, iter=100, chains=1,
                include=FALSE,
                pars=c('rawUnique', 'rawUniqueTheta', 'rawPerComponentVar',
                       'rawFactor', 'rawLoadings', 'rawFactorProp', 'rawNegateFactor', 'rawSeenFactor',
                       'unique', 'uniqueTheta')))
  itemModelExplorer(dl, fit, "i3")
})
