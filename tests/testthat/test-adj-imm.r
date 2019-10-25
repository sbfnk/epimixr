context("Adjusting immunity levels")

mixing <- socialmixr::contact_matrix(survey=socialmixr::polymod, age.limits = c(0, 5, 10))
mixing2 <- socialmixr::contact_matrix(n=2, survey=socialmixr::polymod, age.limits = c(0, 5, 10))

test_that("immunity levels can be projected",
{
  baseline.immunity <- c(`2`=0.85, `5`=0.9, `10`=0.95)
  coverage <- matrix(rep(0.9, 10), nrow=2)
  colnames(coverage) <- as.character(seq(2015, 2019))
  expect_equal(unname(round(project_immunity(baseline.immunity, 2018, 2019, coverage = coverage, schedule=c(1, 2), 0.5, 0.95), 1)), c(0.7, 0.9, 0.9, 1.0))
  expect_equal(unname(round(project_immunity(baseline.immunity, 2018, 2019, maternal.immunity = 0.5), 1)), c(0.5, 0.8, 0.9, 1.0))
})

test_that("errors are thrown correctly",
{
  expect_error(project_immunity(), "baseline immunity")
  expect_error(project_immunity(0.9), "baseline year")
  expect_error(project_immunity(0.9, 2000), "'year'")
  expect_error(project_immunity(0.9, 2000, 1998), "must be greater")
  expect_error(project_immunity(0.9, 2000, 2004), "maternal immunity")
  expect_error(project_immunity(0.9, 2000, 2004, matrix(c(0.9, 0.9), nrow=1)), "'schedule' must be given")
  expect_error(project_immunity(0.9, 2000, 2004, matrix(c(0.9, 0.9), nrow=1), c(1, 5)), "a row for each element")
  expect_error(project_immunity(0.9, 2000, 2004, matrix(c(0.9, 0.9), nrow=1), c(1), 0.5), "efficacy")
})

test_that("adjusted immunity levels can be calculated",
{
  expect_equal(round(adjust_immunity(mixing, immunity = c(0, 0.5, 0.8)), 1), 0.7)
  expect_equal(round(adjust_immunity(mixing$matrix, immunity = c(0, 0.5, 0.8)), 1), 0.7)
  expect_length(adjust_immunity(mixing2, immunity = c(0, 0.5, 0.8)), 2)
})

