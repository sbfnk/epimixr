##' Estimate reproduction number from contact information and immunity profile
##'
##' This takes a contact survey to derive a contact matrix and rescales contacts
##' to represent contacts with susceptibles. This is then combined with
##' information on the basic reproduction number R0 to calculate the effective
##' or net reproduction number.
##' @param mixing_matrix A mixing matrix, as returned by
##'   \code{socialmixr::contact_matrix}
##' @param immunity immunity profile; this should be given as a vector of the
##'   same length as the number of rows/columns of the mixing matrix; each
##'   element of the vector should contain a value <1 representing the
##'   proportion of the population immune in the corresponding age group; any
##'   element set to "herd" will be set to 1-1/R0
##' @param vector if TRUE, will return the eigenvector corresponding to the
##'   dominant eigenvector instead of adjusted immunity; this corresponds to the
##'   expected stable age distribution of infections in case of an outbreak
##' @return a list contain vectors of adjusted immunities
##' @importFrom socialmixr contact_matrix
##' @author Sebastian Funk
##' @export
##' @examples
##' library("socialmixr")
##' mixing <- contact_matrix(survey = polymod, age.limits
##'   = c(0, 5, 10))
##' adjust_immunity(mixing$matrix, immunity = c(0, 0.5, 0.8))
adjust_immunity <- function(mixing_matrix, immunity, vector = FALSE) {
  ## rescale by immunity
  mixing_immunised <- mixing_matrix %*% diag(1 - immunity)
  if (vector) {
    v <- eigen(mixing_immunised)$vectors[, 1]
    return(v / sum(v))
  } else {
    ## calculate adjImm = 1 - ratio of largest eigenvalues between immunised /
    ## completely susceptible
    return(1 - Re(eigen(mixing_immunised)$values[1]) /
      Re(eigen(mixing_matrix)$values[1]))
  }
}
