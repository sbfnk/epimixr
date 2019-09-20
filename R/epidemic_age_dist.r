##' calculates the age distribution in an epidemic setting using the iterative method
##' of Wallinga (2006)
##'
##' calculates the age distribution in an epidemic setting using the iterative method
##' of Wallinga (2006)
##' @param contact.matrix: social contact matrix
##' @param q: probability of transmission per contact
##' @param child.mixing: first entry of contact matrix
##' @param immunity: proportion immune before the epidemic
##' @param final.size.start.start: starting value for inidence
##' @param tol: tolerance for stopping the iteration
##' @export
##' @return A vector of the final size (proportion of susceptibles infected) in each age group
epidemic_age_dist <- function(contact.matrix, q,
                              child.mixing = NA,
                              immunity = 0,
                              final.size.start = 0.01,
                              tol = 1e-5) {
  ## initialise variables
  z <- rep(final.size.start, nrow(contact.matrix))
  last.z <- rep(0, nrow(contact.matrix))
  first.run <- TRUE
  if (!is.na(child.mixing)) {
    contact.matrix[1,1] <- child.mixing
  }
  if (length(immunity) == 1) {
    immunity <- rep(immunity, nrow(contact.matrix))
  }
  contact.matrix <- t(apply(contact.matrix, 1, function(x) {
    x * (1 - immunity)
  }))

  ## set to greater than tol for first time the loop is run
  current.diff <- tol + 1

  ## loop until difference between estimates is smaller than tolerance
  while (current.diff > tol) {
    rhs <- 1 - exp(- z %*% (q * contact.matrix))
    if (first.run == TRUE) {
      ## run loop at least two times
      current.diff <- tol + 1
      first.run <- FALSE
    } else {
      current.diff <- sum(abs(rhs - z))
    }
    z <- rhs
  }
  return(z[1, ])
}
