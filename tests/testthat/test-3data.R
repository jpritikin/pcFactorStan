library(testthat)
context("test-3data")

suppressWarnings(RNGversion("3.5"))

test_that("normalizeData", {
  df <- data.frame(pa1=NA, pa2=NA, i1=c(1, -1), i2=c(-2, 2))
  df[1,paste0('pa',1:2)] <- c('a','b')
  df[2,paste0('pa',1:2)] <- c('b','a')
  expect_error(normalizeData(df, bob="whatever"),
               "Rejected are any values passed")
  df <- normalizeData(df)
  expect_equal(df$pa1, c('a','a'))
  expect_equal(df$i1, c(1,1))
  expect_equal(df$i2, c(-2,-2))

  df <- rbind(df, data.frame(pa1='c', pa2='b', i1=2, i2=-1))
  df <- normalizeData(df, .palist = c('c','b','a'))
  expect_equal(df$pa1, c( 'b','b','c'))
  expect_equal(df$i1, c(-1,-1,2))

  expect_error(normalizeData(df, .palist = c('b','a')),
               ".palist must be length 3")
  expect_error(normalizeData(df, .palist = c('b','a','q')),
               ".palist must contain the names of all vertices")

  d1 <- phyActFlowPropensity
  d2 <- normalizeData(d1, .sortRows = FALSE)
  # due to truncation from "running;solo" and similar
  expect_equal(which(d1$pa1 != d2$pa1), c(151, 153))
})

test_that("prepCleanData", {
  df <- data.frame(pa1=NA, pa2=NA, i1=c(1, -1))
  df[1,paste0('pa',1:2)] <- c('a','b')
  df[2,paste0('pa',1:2)] <- c('b','a')
  dl <- prepCleanData(df)
  expect_equal(dl$pa1, c(1,2))
  expect_equal(dl$pa2, c(2,1))

  df <- phyActFlowPropensity[,c(1,2,7)]
  dl <- prepCleanData(df)
  expect_equivalent(c(table(dl$refresh)), c(216L, 70L, 25L, 21L, 10L))
  expect_equivalent(c(table(dl$weight)),
                    c(405L, 79L, 33L, 10L, 12L, 6L, 6L, 4L, 1L, 1L, 1L, 2L, 1L, 1L,  2L, 1L))
})

test_that("findModel", {
  expect_message(findModel(), "Models available")
  expect_error(findModel("asdg"), "Stan model not found")
})

test_that("prepSingleFactorModel", {
  dl <- prepData(phyActFlowPropensity)
  dl <- prepSingleFactorModel(dl)

  expect_equal(dl$NFACTORS, 1)
  expect_equal(dl$NPATHS, dl$NITEMS)
  expect_equal(dl$factorScalePrior, as.array(1.2))
  expect_equal(dl$factorItemPath[1,],
               rep(1, dl$NITEMS))
  expect_equal(dl$factorItemPath[2,],
               1:dl$NITEMS)
})

test_that("prepFactorModel", {
  pa <- phyActFlowPropensity[,setdiff(colnames(phyActFlowPropensity),
                                      c('goal1','feedback1'))]
  dl <- prepData(pa)
  expect_error(prepFactorModel(dl, list()),
               "No paths specified")
  expect_error(prepFactorModel(dl,
                        list(c('complex','skill','predict',
                                    'creative', 'novelty', 'stakes',
                                    'present', 'reward', 'chatter',
                                    'body'),
                             c('waiting','control','evaluated','spont'),
                             c('novelty', 'waiting'))),
               "paths must be named")
  expect_error(prepFactorModel(dl,
                               list(flow=c('complex','skill','predict',
                                           'creative', 'novelty', 'stakes',
                                           'present', 'reward', 'chatter',
                                           'body'),
                                    f2=c('waiting','control','evaluated','spont'),
                                    rc=c('novelty'))),
               "less than 2 indicators")
  expect_error(prepFactorModel(dl,
                               list(flow=c('complex','skill','predict',
                                           'creative', 'novelty', 'stakes',
                                           'present', 'reward', 'chatter',
                                           'bodies'),
                                    f2=c('waiting','control','evaluated','spont'),
                                    rc=c('novelty', 'waiting'))),
               "No matching item for factor")
  expect_error(prepFactorModel(dl,
                               list(flow=c('complex','skill','predict',
                                           'creative', 'novelty', 'stakes',
                                           'present', 'reward', 'chatter'),
                                    f2=c('waiting','control','evaluated','spont'),
                                    rc=c('novelty', 'waiting'))),
               "No factor predicts item")
  dl <- prepFactorModel(dl,
                        list(flow=c('complex','skill','predict',
                                    'creative', 'novelty', 'stakes',
                                    'present', 'reward', 'chatter',
                                    'body'),
                             f2=c('waiting','control','evaluated','spont'),
                             rc=c('novelty', 'waiting')))
  expect_equal(dl$nameInfo$factor,
               c('flow','f2','rc'))
  expect_equal(dl$NFACTORS, 3)
  expect_equal(dl$NPATHS, 16)
  expect_equal(dl$factorScalePrior, as.array(rep(1.2,3)))
})
