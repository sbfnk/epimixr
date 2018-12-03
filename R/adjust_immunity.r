##' Estimate reproduction number from contact information and immunity profile
##'
##' This takes a contact survey to derive a contact matrix and rescales contacts to represent contacts with susceptibles. This is then combined with information on the basic reproduction number R0 to calculate the effective or net reproduction number.
##' @param mixing A mixing matrix or set of mixing matrices, as returned by \code{\link{contact_matrix}}
##' @param demography A data frame with a \code{population} column, indicating the populations of the age groups given in \code{immunity}; if it is not given, it will be assumed that the mixing matrix is already normalised
##' @param immunity immunity profiles; this should be given as a vector of the same length as the number of rows/columns of the mixing matrix (or matrices); if given as a vector, each element of the vector should contain a value <1 representing the proportion of the population immune in the corresponding age group; any element set to "herd" will be set to 1-1/R0; if given as a data.frame, the number of rows need to be \code{n}, and every row will correspond to one of the samples
##' @return a list contain vectors of adjusted immunities
##' @importFrom socialmixr contact_matrix
##' @author Sebastian Funk
##' @export
adjust_immunity <- function(mixing, demography, immunity)
{
    ret <- c()

    if (!("matrices" %in% names(mixing)))
    {
        mixing$matrices <- list(mixing$matrix)
    }
    if (is.vector(immunity))
    {
        immunity <- matrix(rep(immunity, length(mixing$matrices)), byrow=TRUE,
                           nrow=length(mixing$matrices))
    }
    demography_given <- !missing(demography)
    ret <- lapply(seq_along(mixing$matrices), function(x)
    {
        if (demography_given) {
            ## rescale to a_{ij} n_j delta_i
            mixing_normalised <-
                t(t(mixing$matrices[[x]] * demography$population) /
                  demography$population)
        } else {
            mixing_normalised <- mixing$matrices[[x]]
        }
        ## rescale by immunity
        mixing_immunised <- mixing_normalised * (1 - immunity[x, ])
        ## calculate adjImm = 1 -
        return(1 - Re(eigen(mixing_immunised)$values[1])/
               Re(eigen(mixing_normalised)$values[1]))
    })

    return(ret)
}
