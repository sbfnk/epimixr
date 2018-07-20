##' Estimate reproduction number from contact information and immunity profile
##'
##' This takes a contact survey to derive a contact matrix and rescales contacts to represent contacts with susceptibles. This is then combined with information on the basic reproduction number R0 to calculate the effective or net reproduction number.
##' @param immunity immunity profiles; this should be given as a vector of the same length as \code{age.limits} or a data.frame; if given as a vector, each element of the vector should contain a value <1 representing the proportion of the population immune in the corresponding age group; any element set to "herd" will be set to 1-1/R0; if given as a data.frame, the number of rows need to be \code{n}, and every row will correspond to one of the samples
##' @param model "age-specific" (default) or "homogeneous"; if set to "homogeneous", the contact survey is not used to generate a contact matrix, but only for demographic information; beyond demography, random mixing is assumed.
##' @param ... any further parameters for \code{\link{contact_matrix}}
##' @return a list contain vectors of adjusted immunities
##' @importFrom socialmixr contact_matrix
##' @author Sebastian Funk
##' @export
adjust_immunity <- function(survey, countries, immunity, n = 1, age.limits, model = c("age-specific", "homogeneous"), sample, ...)
{
    if (missing(survey))
    {
        stop("'survey' missing.")
    }

    if (missing(countries)) countries <- c()

    norm <- list()

    if (!missing(immunity) && is.vector(immunity)) {
      limits <- names(immunity)
      immunity <- data.frame(n=1:n, t(immunity))
      immunity$n <- NULL
      colnames(immunity) <- limits
    }

    if (missing(age.limits))
    {
        if (missing(immunity) || any(is.na(as.integer(colnames(immunity)))))
        {
            stop("Either 'age.limits' or 'immunity' must be given (as integers), otherwise I don't know what age groups to use.")
        } else age.limits <- as.integer(colnames(immunity))
    } else {
        age.limits <- as.integer(age.limits)
    }

    if (any(is.na(age.limits)) || any(diff(age.limits) <= 0))
    {
        stop("'age.limits' must be given (either directly or as names of the 'immunity' vector) as an increasing ",
             "integer vector of lower age limits.")
    }

    ret <- c()
    if (model == "age-specific")
    {
        mixing <-
            contact_matrix(survey, countries = countries, n = n, quiet = TRUE, age.limits = age.limits, ...)
        if (n == 1)
        {
            mixing$matrices <- list(mixing)
        }
        for (i in seq_along(mixing$matrices))
        {
          mixing_i <- mixing$matrices[[i]]
          ## rescale to a_{ij} n_j delta_i
          mixing_normalised <- t(t(mixing_i$matrix * mixing$demography$proportion) / mixing$demography$proportion)
          ## rescale by immunity
          mixing_immunised <- mixing_normalised * (1 - unlist(immunity[i, ]))
          ## calculate adjImm = 1 -
          ret <- c(ret, 1 - Re(eigen(mixing_immunised)$values[1])/Re(eigen(mixing_normalised)$values[1]))
        }
    } else if (model == "homogeneous")
    {
        ## create a mixing matrix to get the correct demography
        mixing <- contact_matrix(survey, split = TRUE, n = 1, quiet = TRUE, age.limits = age.limits, ...)
        ## calculate homogeneous probability v of being immune
        immunity_hom <-
            apply(immunity, 1, function(x) sum(mixing$demography$population * x)) /
            sum(mixing$demography$population)
        ## calculate R = (1-v) R_0
        ret <- c(ret, immunity_hom)
    }

    return(ret)
}

adjust_immunity_wrapper <- function(...)
{
    adjust_immunity(...)
}
