context("Calculating age distributions")

mixing <- socialmixr::contact_matrix(survey=socialmixr::polymod, age.limits = c(0, 5, 10))
mixing2 <- socialmixr::contact_matrix(n=2, survey=socialmixr::polymod, age.limits = c(0, 5, 10))

test_that("epidemic age distributions can be calculated",
{
  expect_equal(round(unname(epidemic_age_dist(mixing, R0=5, immunity=0.50)[1, ]), 1), c(0.2, 0.5, 1))
  expect_equal(round(unname(epidemic_age_dist(mixing$matrix, R0=5, immunity=0.50)[1, ]), 1), c(0.2, 0.5, 1))
  expect_equal(nrow(epidemic_age_dist(mixing2, R0=5, immunity=0.50)), 2)
})

