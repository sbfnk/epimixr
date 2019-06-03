##' Project immunity from a baseline via vaccination coverage rates
##'
##' @title Project immunity from a baseline
##' @param baseline.immunity baseline immunity, as a named vector; the names correspond to lower limits of the age groups, and the vector itself to the corresponding levels of immunity.
##' @param coverage coverage with multiple vaccine doses, given as a matrix in which each row is a dose and each column a year; the first column is the coverage in the year the serological study was conducted
##' @param schedule the ages at which vaccines are given (in years).
##' @param maternal.immunity the proportion maternally immune.
##' @param efficacy vaccine efficacy.
##' @return a data frame of immunity levels by age group (as in \code{baseline.immunity}).
##' @author Sebastian Funk <sebastian.funk@lshtm.ac.uk>
project_immunity <- function(baseline.immunity, coverage, schedule, maternal.immunity, efficacy)
{
  ## checks
  if (missing(baseline.immunity)) stop("'coverage' must be provided")
  if (missing(coverage)) stop("'coverage' must be provided")
  if (length(unique(vapply(coverage, length, 0L))) > 1)
  {
    stop("All elements of 'coverage' must be the same length")
  }

  ## convert baseline to annual immunity
  lower_age_limits <- as.integer(names(baseline.immunity))
  bdf <- data.frame(lower.age.limit=lower_age_limits,
                    immunity=baseline.immunity)
  min_age <- min(bdf$lower.age.limit)
  max_age <- max(bdf$lower.age.limit)
  df <- data.frame(lower.age.limit=
                     do.call(seq, as.list(range(lower_age_limits))))
  df[df$lower.age.limit %in% bdf$lower.age.limit, "immunity"] <- bdf$immunity
  ## fill
  df$immunity <- c(NA, na.omit(df$immunity))[cumsum(!is.na(df$immunity))+1]

  scaling_factor <- min(df[df$lower.age.limit==1, "immunity"] / coverage[[1]][1], 1)
  oldest <- df[nrow(df), ]

  if (length(coverage[[1]]) > 1) {
    for (i in seq(2, length(coverage[[1]]))) {
      ## move all one age group up
      df$lower.age.limit <- df$lower.age.limit + 1
      ## implement vaccination schedule
      df <- df[df$lower.age.limit > 2, ]
      first_two_years <-
        data.frame(lower.age.limit=c(0, 1, 2),
                   immunity=c(maternal.immunity,
                              coverage[[1]][i-1] * scaling_factor * efficacy,
                              coverage[[1]][i] * efficacy))
      df <- rbind(first_two_years, df)
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
  df <- do.call(rbind, summarised)
  return(df)
}
