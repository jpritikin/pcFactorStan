context("test-generate")

RNGversion("3.5")

test_that("generateItem", {
  set.seed(1)
  df <- roundRobinGraph(letters[1:5], 40)
  df <- generateItem(df)
  expect_equivalent(c(table(df$i1)), c(16, 17, 7))

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
