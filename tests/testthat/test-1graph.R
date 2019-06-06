context("test-1graph")

RNGversion("3.5")

test_that("roundRobinGraph", {
  df <- roundRobinGraph(letters[1:5], 15)
  expect_identical(paste0(df$pa1, collapse = ''),
                   "aababcabcdaabab")
  expect_identical(paste0(df$pa2, collapse = ''),
                   "bccdddeeeebccdd")

  expect_error(roundRobinGraph("a", 15),
               "Must provide at least 2 names")
  expect_warning(roundRobinGraph(letters, 5),
               "Sample size too small 5 to round robin connect all the vertices 325")
})

test_that("twoLevelGraph", {
  set.seed(1)
  df <- twoLevelGraph(letters[1:8], 20, .8, .5)
  expect_identical(paste0(df$pa1, collapse = ''),
                   "aaaaaaabaabaacaabbea")
  expect_identical(paste0(df$pa2, collapse = ''),
                   "bcdefghfbgebfdbbcchf")

  expect_error(twoLevelGraph("a", 15, .8, .5),
               "Must provide at least 2 names")
  expect_warning(twoLevelGraph(letters, 5, .8, .5),
                 "Sample size too small 5 to connect all the vertices 26")
  expect_error(twoLevelGraph(letters, 50, -1, .5),
               "must be between 0 and 1")
  expect_error(twoLevelGraph(letters, 50, .5, 2),
               "must be between 0 and 1")
})

test_that("filterGraph", {
  junk <- data.frame(bob=rep(1,3), alice=rep(2,3))
  expect_error(filterGraph(junk),
               "Expected columns pa1, pa2 to contain the vertex names")
  junk <- data.frame(pa1=factor(c('a','b')),
                     pa2=factor(c('b','c')))
  expect_error(filterGraph(junk),
               "Some levels(pa1) don't match levels(pa2)", fixed=TRUE)

  df <- data.frame(pa1=rep(NA,15), pa2=NA)
  for (rx in 1:11) {
    df[rx,] <- c('a','b')
  }
  for (rx in 12:15) {
    df[rx,] <- c('a','c')
  }
  df <- filterGraph(df)
  expect_equal(nrow(df), 11)
  expect_equal(attr(df, 'weak'), "c")
  expect_equal(levels(df$pa1), c('a','b'))
  expect_message(print(df),
                 'too weakly connected: "c"')

  df <- data.frame(pa1=rep(NA,15), pa2=NA)
  for (rx in 1:11) {
    df[rx,] <- c('a','b')
  }
  for (rx in 12:14) {
    df[rx,] <- c('b','c')
  }
  df[15,] <- c('a','c')
  df <- filterGraph(df)
  expect_equal(nrow(df), 15)
  expect_equal(length(levels(df$pa1)), 3)

  df <- filterGraph(phyActFlowPropensity[,c(1,2)])
  expect_equal(length(attr(df, 'disconnected')), 8)
  expect_equal(length(attr(df, 'weak')), 18)
  expect_equal(length(levels(df$pa1)), 61)
  expect_message(print(head(df)),
                 "not part of the largest connected component")
  expect_message(print(head(df)),
                 "too weakly connected")
})
