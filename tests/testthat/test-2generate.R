context("test-2generate")

RNGversion("3.5")

test_that("generateItem", {
  set.seed(1)
  df <- roundRobinGraph(letters[1:5], 40)
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

test_that("generateFactorItems", {
  set.seed(1)
  df <- twoLevelGraph(letters[1:10], 100)
  df <- generateFactorItems(df, 3)

  # This is a nonsensical way to look at the data.
  # Just ensure that nothing has changed.
  c1 <- cov(df[,paste0('i',1:3)])
  expect_equal(c1[lower.tri(c1, diag = TRUE)],
               c(0.674, -0.14, -0.107, 0.672, 0.229, 0.559),
               tolerance=1e-3, scale=1)

  expect_error(generateFactorItems(df, 1),
               "At least 3 indicators are required")
  expect_error(generateFactorItems(df, c(.3,.4)),
               "At least 3 indicators are required")
  expect_error(generateFactorItems(df, c(1.3,.4,.4)),
               "Proportions must be between 0 and 1")
})
