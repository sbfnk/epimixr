##' Calculates the age distribution of an epidemic
##'
##' calculates the age distribution in an epidemic setting using the iterative method
##' of: J Wallinga, P Teunis, M Kretschmar (2006)  Using Data on Social Contacts to Estimate Age-specific Transmission Parameters for Respiratory-spread Infectious Agents. Am J Epidemiol 164(10), 945-946.
##' @param mixing A mixing matrix or set of mixing matrices, as returned by \code{\link{contact_matrix}}
##' @param R0 basic reproduction number
##' @param immunity proportion immune before the epidemic
##' @param final.size.start starting value for inidence
##' @param tol tolerance for stopping the iteration
##' @export
##' @return A matrix of the final size(s) (proportion of susceptibles infected) in each age group (one row per matrix contained in \code{mixing})
##' @examples
##' mixing <- socialmixr::contact_matrix(survey=socialmixr::polymod, age.limits = c(0, 5, 10))
##' epidemic_age_dist(mixing, R0=5, immunity=0.50)
epidemic_age_dist <- function(mixing, R0,
                              immunity = 0,
                              final.size.start = 0.01,
                              tol = 1e-5)
{

  if (!("list" %in% class(mixing))) {
    mixing <- list(matrices=list(list(matrix=mixing)))
  } else if (!("matrices" %in% names(mixing)))
  {
    mixing$matrices <- list(list(matrix=mixing$matrix))
  }

  ret <- lapply(seq_along(mixing$matrices), function(x)
  {
    ## initialise variables
    contact.matrix <- mixing$matrices[[x]][["matrix"]]

    z <- rep(final.size.start, nrow(contact.matrix))
    last.z <- rep(0, nrow(contact.matrix))
    first.run <- TRUE

    ## calculate next-generation matrix
    ngm <- R0 * contact.matrix / eigen(contact.matrix)$values[1]

    ## distribute immunity across age groups
    if (length(immunity) == 1) {
      immunity <- rep(immunity, nrow(contact.matrix))
    }

    ## initialise variables
    ## correct for immunity
    ingm <- ngm %*% diag(1 - immunity)

    ## set to greater than tol for first time the loop is run
    current.diff <- tol + 1

    ## loop until difference between estimates is smaller than tolerance
    while (current.diff > tol) {
      rhs <- 1 - exp(- z %*% (ingm))
      if (first.run == TRUE) {
        ## run loop at least two times
        current.diff <- tol + 1
        first.run <- FALSE
      } else {
        current.diff <- sum(abs(rhs - z))
      }
      z <- rhs
    }
    colnames(z) <- colnames(contact.matrix)
    return(z)
  })

  return(as.matrix(do.call(rbind, ret)))
}
