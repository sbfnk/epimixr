##' Estimate reproduction number from contact information and immunity profile
##'
##' This takes a contact survey to derive a contact matrix and rescales contacts to represent contacts with susceptibles. This is then combined with information on the basic reproduction number R0 to calculate the effective or net reproduction number.
##' @param mixing A mixing matrix or set of mixing matrices, as returned by \code{\link{contact_matrix}}
##' @param immunity immunity profiles; this should be given as a vector of the same length as the number of rows/columns of the mixing matrix (or matrices); if given as a vector, each element of the vector should contain a value <1 representing the proportion of the population immune in the corresponding age group; any element set to "herd" will be set to 1-1/R0; if given as a data.frame, the number of rows need to be \code{n}, and every row will correspond to one of the samples
##' @return a list contain vectors of adjusted immunities
##' @importFrom socialmixr contact_matrix
##' @author Sebastian Funk
##' @export
##' @examples
##' mixing <- socialmixr::contact_matrix(survey=socialmixr::polymod, age.limits = c(0, 5, 10))
##' adjust_immunity(mixing, immunity = c(0, 0.5, 0.8))
adjust_immunity <- function(mixing, immunity)
{
    ret <- c()

    if (!("list" %in% class(mixing))) {
        mixing <- list(matrices=list(list(matrix=mixing)))
    } else if (!("matrices" %in% names(mixing)))
    {
        mixing$matrices <- list(list(matrix=mixing$matrix))
    }
    if (is.vector(immunity))
    {
        immunity <- matrix(rep(immunity, length(mixing$matrices)), byrow=TRUE,
                           nrow=length(mixing$matrices))
    }
    ret <- lapply(seq_along(mixing$matrices), function(x)
    {
        contact.matrix <- mixing$matrices[[x]][["matrix"]]
        ## rescale by immunity
        mixing_immunised <- contact.matrix %*% diag(1 - immunity[x, ])
        ## calculate adjImm = 1 -
        return(1 - Re(eigen(mixing_immunised)$values[1])/
               Re(eigen(contact.matrix)$values[1]))
    })

    return(unlist(ret))
}
