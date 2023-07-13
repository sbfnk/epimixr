context("Adjusting immunity levels")

library("socialmixr")

mixing <- contact_matrix(survey = polymod, age.limits = c(0, 5, 10))

test_that("immunity levels can be projected", {
  baseline_immunity <- c(`2` = 0.85, `5` = 0.9, `10` = 0.95)
  coverage <- matrix(rep(0.9, 10), nrow = 2)
  colnames(coverage) <- as.character(seq(2015, 2019))
  expect_equal(unname(round(
    project_immunity(
      baseline_immunity, 2018, 2019,
      coverage = coverage, schedule = c(1, 2),
      0.5, 0.95
    ), 2
  )), c(0.68, 0.89, 0.89, 0.95))
  expect_equal(unname(round(
    project_immunity(
      baseline_immunity, 2018, 2019,
      maternal_immunity = 0.5
    ), 2
  )), c(0.50, 0.85, 0.90, 0.95))
})

test_that("errors are thrown correctly", {
  expect_error(project_immunity(), "baseline immunity")
  expect_error(project_immunity(0.9), "baseline year")
  expect_error(project_immunity(0.9, 2000), "'year'")
  expect_error(project_immunity(0.9, 2000, 1998), "must be greater")
  expect_error(project_immunity(0.9, 2000, 2004), "maternal immunity")
  expect_error(project_immunity(
    0.9, 2000, 2004, matrix(c(0.9, 0.9), nrow = 1)
  ), "'schedule' must be given")
  expect_error(project_immunity(
    0.9, 2000, 2004, matrix(c(0.9, 0.9), nrow = 1), c(1, 5)
  ), "a row for each element")
  expect_error(project_immunity(
    0.9, 2000, 2004, matrix(c(0.9, 0.9), nrow = 1), c(1), 0.5
  ), "efficacy")
})

test_that("adjusted immunity levels can be calculated", {
  expect_equal(round(adjust_immunity(
    mixing$matrix,
    immunity = c(0, 0.5, 0.8)
  ), 1), 0.7)
  expect_equal(round(adjust_immunity(
    mixing$matrix,
    immunity = c(0, 0.5, 0.8)
  ), 1), 0.7)
})
