context("test-3data")

suppressWarnings(RNGversion("3.5"))

test_that("normalizeData", {
  df <- data.frame(pa1=NA, pa2=NA, i1=c(1, -1), i2=c(-2, 2))
  df[1,paste0('pa',1:2)] <- c('a','b')
  df[2,paste0('pa',1:2)] <- c('b','a')
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
  expect_equivalent(c(table(dl$refresh)), c(223, 342))
  expect_equivalent(c(table(dl$weight)),
                    c(405L, 79L, 33L, 10L, 12L, 6L, 6L, 4L, 1L, 1L, 1L, 2L, 1L, 1L,  2L, 1L))
})

test_that("findModel", {
  expect_message(findModel(), "Models available")
  expect_error(findModel("asdg"), "Stan model not found")
})
