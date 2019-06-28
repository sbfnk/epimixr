##' Project immunity from a baseline via vaccination coverage rates
##'
##' @title Project immunity from a baseline
##' @param baseline.immunity baseline immunity, as a named vector; the names correspond to lower limits of the age groups, and the vector itself to the corresponding levels of immunity.
##' @param baseline.year year at which baseline immunity is taken (corresponding to a column in the \code{coverage} argument)
##' @param year year to project to
##' @param coverage coverage with multiple vaccine doses, given as a matrix in which each row is a dose and each (named) column a year
##' @param schedule the ages at which vaccines are given (in years).
##' @param maternal.immunity the proportion maternally immune.
##' @param efficacy vaccine efficacy.
##' @return a data frame of immunity levels by age group (as in \code{baseline.immunity}).
##' @author Sebastian Funk <sebastian.funk@lshtm.ac.uk>
##' @export
project_immunity <- function(baseline.immunity, baseline.year, year, coverage, schedule, maternal.immunity, efficacy)
{
  ## checks
  if (missing(baseline.immunity)) stop("baseline immunity must be provided")
  if (missing(baseline.year)) stop("baseline year must be provided")
  if (missing(year)) stop("'year' argument must be provided")
  if (!(year > baseline.year)) stop("'year' must be greater than 'baseline.year'")
  if (!missing(coverage) && dim(coverage)[1] != length(schedule))
  {
    stop("'coverage' must have a row for each element of 'schedule'")
  }

  ## convert baseline to annual immunity
  lower_age_limits <- as.integer(names(baseline.immunity))
  bdf <- data.frame(lower.age.limit=lower_age_limits,
                    immunity=baseline.immunity)
  df <- data.frame(lower.age.limit=
                     do.call(seq, as.list(range(lower_age_limits))))
  df[df$lower.age.limit %in% bdf$lower.age.limit, "immunity"] <- bdf$immunity
  ## fill
  df$immunity <- c(NA, na.omit(df$immunity))[cumsum(!is.na(df$immunity))+1]

  if (missing(coverage))
  {
    df <- df[!(df$lower.age.limit==0), ]
    df <- rbind(t(c(lower.age.limit=0, immunity=maternal.immunity)), df)
  } else
  {
    oldest <- df[nrow(df), ]
    min_age <- min(df$lower.age.limit)

    ## fill missing age groups with vaccination data
    while (min_age > schedule[1])
    {
      min_age <- min_age - 1
      df <-
        rbind(t(c(lower.age.limit=min_age,
                  immunity=coverage[1, as.character(baseline.year - min_age + 1)] *
                    efficacy)),
              df)
    }

    scaling_factor <-
      min(df[df$lower.age.limit==schedule[1], "immunity"] /
          coverage[1, as.character(baseline.year)], 1)

    if (dim(coverage)[2] > 1)
    {
      for (calc.year in seq(baseline.year + 1, year))
      {
        ## move all one age group up
        df$lower.age.limit <- df$lower.age.limit + 1
        ## implement vaccination schedule
        df <- df[df$lower.age.limit > schedule[1] + 1, ]
        first_years <-
          data.frame(lower.age.limit=seq(0, schedule[1] + 1),
                     immunity=c(rep(maternal.immunity, schedule[1]),
                                coverage[1, as.character(calc.year)] *
                                scaling_factor,
                                coverage[1, as.character(calc.year - 1)] * efficacy))
        df <- rbind(first_years, df)
        if (dim(coverage)[1] > 1)
        {
          for (j in seq(2, dim(coverage)[1]))
          {
            immunised <- 0
            for (k in seq(1, j-1)) {
              old_coverage <-
                coverage[k, as.character(calc.year - schedule[j] + schedule[k])]
              immunised <- immunised + (1 - immunised) * old_coverage * efficacy
            }
            df[df$lower.age.limit==schedule[j], "immunity"] <-
              min(1,
                  df[df$lower.age.limit==schedule[j], "immunity"] +
                  (1-immunised) * coverage[j, as.character(calc.year)] * efficacy)
          }
        }
      }
    }
    df <- df[df$lower.age.limit < oldest$lower.age.limit, ]
    df <- rbind(df, oldest)
    rownames(df) <- seq_len(nrow(df))
  }

  ## aggregate by age groups
  df$lower.age.limit <-
    socialmixr::reduce_agegroups(df$lower.age.limit, union(0, lower_age_limits))

  summarised <- by(df, list(df$lower.age.limit), function(x)
    c(
      lower.age.limit=unique(x$lower.age.limit),
      immunity=mean(x$immunity)
    ))
  df <- as.data.frame(do.call(rbind, summarised))
  ret <- df$immunity
  names(ret) <- df$lower.age.limit
  return(ret)
}
