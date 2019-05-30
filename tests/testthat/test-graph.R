context("test-graph")

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
