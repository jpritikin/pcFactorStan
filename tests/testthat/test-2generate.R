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
  expect_equivalent(c(table(df$i1)), c(10,11,19))

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
  # cat(deparse(round(c1[lower.tri(c1, diag = TRUE)],3)))
  expect_equal(c1[lower.tri(c1, diag = TRUE)],
               c(0.794, 0.101, -0.232, 0.829, -0.219, 0.822),
               tolerance=.001, scale=1)

  expect_error(generateCovItems(df, 1),
               "numItems must be 2 or greater")
  colnames(df)[3] <- 'i4'
  expect_error(generateCovItems(df, 3),
               "Colname i4 is already taken.")
})

test_that("generateSingleFactorItems", {
  set.seed(1)
  df <- expect_warning(twoLevelGraph(letters[1:10], 4),
                 "Sample size too small")
  expect_equal(nrow(df), 4)

  df <- twoLevelGraph(letters[1:10], 100)
  df <- generateSingleFactorItems(df, 3)

  # This is a nonsensical way to look at the data.
  # Just ensure that nothing has changed.
  c1 <- cov(df[,paste0('i',1:3)])
  #  cat(deparse(round(c1[lower.tri(c1, diag = TRUE)],3)))
  expect_equal(c1[lower.tri(c1, diag = TRUE)],
               c(0.782, 0.26, -0.049, 0.829, -0.101, 0.825),
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
  # cat(deparse(round(c1[lower.tri(c1, diag = TRUE)],3)))
  expect_equal(c1[lower.tri(c1, diag = TRUE)],
               c(0.768, -0.024, -0.03, -0.051, 0.854, 0.101, -0.038,
                 0.808,  0.121, 0.735),
               tolerance=1e-3, scale=1)

  expect_error(generateFactorItems(df, list(f1=paste0('i',1:4),
                                     f2=paste0('i',2:4)),
                            c(f1=100, f2=100)),
               "factorScalePrior is too large")
})
