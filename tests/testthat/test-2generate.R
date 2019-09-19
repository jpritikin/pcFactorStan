library(testthat)
library(pcFactorStan)
context("test-2generate")

suppressWarnings(RNGversion("3.5"))

test_that("generateItem", {
  set.seed(1)
  df <- roundRobinGraph(letters[1:5], 40)
  expect_error(generateItem(df, bob="whatever"),
               "Rejected are any values passed")
  df <- generateItem(df)
  expect_equivalent(c(table(df$i1)), c(7,17,16))

  expect_error(generateItem(df, name="i1"),
               "Colname i1 is already taken")
  expect_error(generateItem(df, rnorm(6)),
               "length(theta) 6 != 5", fixed=TRUE)
  expect_error(generateItem(df, rnorm(5)),
               "No latent score for: a b c d e", fixed=TRUE)

  theta <- rnorm(5)
  names(theta) <- letters[2:6]
  expect_error(generateItem(df, theta),
               "No latent score for: a", fixed=TRUE)

  expect_error(generateItem('i1'), "df is not a data.frame")

  theta <- matrix(rnorm(10), nrow=5)
  expect_error(generateItem(df, theta),
               "No latent score")
  dimnames(theta) <- list(letters[1:5], c('i2','i3'))
  expect_error(generateItem(df, theta, name=c('apple', 'banana')),
               "Mismatch between name and colnames(theta)", fixed=TRUE)
  colnames(theta) <- NULL
  df <- generateItem(df, theta)
  expect_equal(colnames(df), c(paste0('pa',1:2), paste0('i',1:3)))
})

test_that("generateCovItems", {
  set.seed(1)
  df <- roundRobinGraph(letters[1:5], 100)
  df <- generateCovItems(df, 3)

  # This is a nonsensical way to look at the data.
  # Just ensure that nothing has changed.
  c1 <- cov(df[,paste0('i',1:3)])
  expect_equal(c1[lower.tri(c1, diag = TRUE)],
               c(0.543, 0.232, -0.134, 0.77, -0.266, 0.652),
               tolerance=.001, scale=1)

  expect_error(generateCovItems(df, 1),
               "numItems must be 2 or greater")
  colnames(df)[3] <- 'i4'
  expect_error(generateCovItems(df, 3),
               "Colname i4 is already taken.")
})

test_that("generateSingleFactorItems", {
  set.seed(1)
  df <- twoLevelGraph(letters[1:10], 100)
  df <- generateSingleFactorItems(df, 3)

  # This is a nonsensical way to look at the data.
  # Just ensure that nothing has changed.
  c1 <- cov(df[,paste0('i',1:3)])
  expect_equal(c1[lower.tri(c1, diag = TRUE)],
               c(0.624, 0.142, 0.009, 0.644, -0.221, 0.579),
               tolerance=1e-3, scale=1)

  expect_error(generateSingleFactorItems(df, 1),
               "At least 3 indicators are required")
  expect_error(generateSingleFactorItems(df, c(.3,.4)),
               "At least 3 indicators are required")
  expect_error(generateSingleFactorItems(df, c(1.3,.4,.4)),
               "Signed proportions must be between -1 and 1")
})

test_that("generateFactorItems", {
  set.seed(1)
  df <- twoLevelGraph(letters[1:10], 100)
  df <- generateFactorItems(df, list(f1=paste0('i',1:4),
                                     f2=paste0('i',2:4)),
                            c(f1=0.9, f2=0.5))
  # This is a nonsensical way to look at the data.
  # Just ensure that nothing has changed.
  c1 <- cov(df[,paste0('i',1:4)])
  expect_equal(c1[lower.tri(c1, diag = TRUE)],
               c(0.657, -0.031, -0.055, -0.042, 0.652, 0.073, 0.043,
                 0.602,  -0.021, 0.563),
               tolerance=1e-3, scale=1)

  expect_error(generateFactorItems(df, list(f1=paste0('i',1:4),
                                     f2=paste0('i',2:4)),
                            c(f1=100, f2=100)),
               "factorScalePrior is too large")
})
