context("Calculating age distributions")

library("socialmixr")

mixing <- contact_matrix(survey = polymod, age.limits = c(0, 5, 10))

test_that("epidemic age distributions can be calculated", {
  expect_equal(round(unname(
    epidemic_age_dist(
      mixing$matrix, r_0 = 5, immunity = 0.50
    )[1, ]
  ), 2), c(0.19, 0.50, 0.96))
})
