##' Calculates the age distribution of an epidemic
##'
##' calculates the age distribution in an epidemic setting using the iterative
##' method of: J Wallinga, P Teunis, M Kretschmar (2006)  Using Data on Social
##' Contacts to Estimate Age-specific Transmission Parameters for
##' Respiratory-spread Infectious Agents. Am J Epidemiol 164(10), 945-946.
##' @param mixing_matrix A mixing matrix or set of mixing matrices, as returned
##'   by \code{\link{socialmixr::contact_matrix}}
##' @param r_0 basic reproduction number
##' @param immunity proportion immune before the epidemic
##' @param final_size_start starting value for inidence
##' @param tol tolerance for stopping the iteration
##' @export
##' @return A matrix of the final size(s) (proportion of susceptibles infected)
##'   in each age group (one row per matrix contained in \code{mixing})
##' @examples
##' library("socialmixr")
##' mixing <- contact_matrix(survey = polymod, age.limits = c(0, 5, 10))
##' epidemic_age_dist(mixing, r_0 = 5, immunity = 0.50)
epidemic_age_dist <- function(mixing_matrix, r_0, immunity = 0,
                              final_size_start = 0.01, tol = 1e-5) {
  ## initialise variables
  z <- rep(final_size_start, nrow(mixing_matrix))
  first_run <- TRUE

  ## calculate next-generation matrix
  ngm <- r_0 * mixing_matrix / eigen(mixing_matrix)$values[1]

  ## distribute immunity across age groups
  if (length(immunity) == 1) {
    immunity <- rep(immunity, nrow(mixing_matrix))
  }

  ## initialise variables
  ## correct for immunity
  ingm <- ngm %*% diag(1 - immunity)

  ## set to greater than tol for first time the loop is run
  current_diff <- tol + 1

  ## loop until difference between estimates is smaller than tolerance
  while (current_diff > tol) {
    rhs <- 1 - exp(-z %*% (ingm))
    if (first_run == TRUE) {
      ## run loop at least two times
      current_diff <- tol + 1
      first_run <- FALSE
    } else {
      current_diff <- sum(abs(rhs - z))
    }
    z <- rhs
  }
  colnames(z) <- colnames(mixing_matrix)

  return(z)
}
