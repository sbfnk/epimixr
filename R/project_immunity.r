##' Project immunity from a baseline via vaccination coverage rates
##'
##' @title Project immunity from a baseline
##' @param baseline_immunity baseline immunity, as a named vector; the names
##'   correspond to lower limits of the age groups, and the vector itself to the
##'   corresponding levels of immunity.
##' @param baseline_year year at which baseline immunity is taken (corresponding
##'   to a column in the \code{coverage} argument)
##' @param year year to project to
##' @param coverage coverage with multiple vaccine doses, given as a matrix in
##'   which each row is a dose and each (named) column a year
##' @param schedule the ages at which vaccines are given (in years).
##' @param maternal_immunity the proportion maternally immune.
##' @param efficacy vaccine efficacy.
##' @return a data frame of immunity levels by age group (as in
##'   \code{baseline_immunity}).
##' @author Sebastian Funk <sebastian.funk@lshtm.ac.uk>
##' @importFrom stats na.omit
##' @importFrom socialmixr reduce_agegroups
##' @export
##' @examples
##' baseline_immunity <- c(`2` = 0.85, `5` = 0.9, `10` = 0.95)
##' coverage <- matrix(rep(0.9, 10), nrow = 2)
##' colnames(coverage) <- as.character(seq(2015, 2019))
##' project_immunity(
##'   baseline_immunity, 2018, 2019, coverage = coverage,
##'   schedule = c(1, 2), 0.5, 0.95
##' )
project_immunity <- function(baseline_immunity, baseline_year, year, coverage,
                             schedule, maternal_immunity, efficacy) {
  ## checks
  if (missing(baseline_immunity)) stop("baseline immunity must be provided")
  if (missing(baseline_year)) stop("baseline year must be provided")
  if (missing(year)) stop("'year' argument must be provided")
  if (!(year > baseline_year)) {
    stop("'year' must be greater than 'baseline_year'")
  }
  if (!missing(coverage)) {
    if (missing(schedule)) stop("'schedule' must be given if 'coverage' is")
    if (dim(coverage)[1] != length(schedule)) {
      stop("'coverage' must have a row for each element of 'schedule'")
    }
    if (missing(efficacy)) stop("'efficacy' must be provided if 'coverage' is")
  }
  if (missing(maternal_immunity)) stop("maternal immunity must be provided")

  ## convert baseline to annual immunity
  lower_age_limits <- as.integer(names(baseline_immunity))
  bdf <- data.frame(
    lower_age_limit = lower_age_limits,
    immunity = baseline_immunity
  )
  df <- data.frame(
    lower_age_limit =
      do.call(seq, as.list(range(lower_age_limits)))
  )
  df[df$lower_age_limit %in% bdf$lower_age_limit, "immunity"] <- bdf$immunity
  ## fill
  df$immunity <- c(NA, na.omit(df$immunity))[cumsum(!is.na(df$immunity)) + 1]

  if (missing(coverage)) {
    df <- df[!(df$lower_age_limit == 0), ]
    df <- rbind(t(c(lower_age_limit = 0, immunity = maternal_immunity)), df)
  } else {
    oldest <- df[nrow(df), ]
    min_age <- min(df$lower_age_limit)

    ## fill missing age groups with vaccination data
    while (min_age > schedule[1]) {
      min_age <- min_age - 1
      df <- rbind(t(c(
        lower_age_limit = min_age,
        immunity = unname(
          coverage[1, as.character(baseline_year - min_age + 1)]
        ) * efficacy
      )), df)
    }

    scaling_factor <- min(df[df$lower_age_limit == schedule[1], "immunity"] /
        coverage[1, as.character(baseline_year)], 1)

    if (dim(coverage)[2] > 1) {
      for (calc.year in seq(baseline_year + 1, year)) {
        ## move all one age group up
        df$lower_age_limit <- df$lower_age_limit + 1
        ## implement vaccination schedule
        df <- df[df$lower_age_limit > schedule[1] + 1, ]
        first_years <- data.frame(
          lower_age_limit = seq(0, schedule[1] + 1),
          immunity = c(
            rep(maternal_immunity, schedule[1]),
            coverage[1, as.character(calc.year)] * scaling_factor,
            coverage[1, as.character(calc.year - 1)] * efficacy
          )
        )
        df <- rbind(first_years, df)
        if (dim(coverage)[1] > 1) {
          for (j in seq(2, dim(coverage)[1])) {
            immunised <- 0
            for (k in seq(1, j - 1)) {
              old_coverage <- coverage[
                k, as.character(calc.year - schedule[j] + schedule[k])
              ]
              immunised <- immunised +
                (1 - immunised) * old_coverage * efficacy
            }
            df[df$lower_age_limit == schedule[j], "immunity"] <- min(
              1,
              df[df$lower_age_limit == schedule[j], "immunity"] +
              (1 - immunised) * coverage[j, as.character(calc.year)] *
              efficacy
            )
          }
        }
      }
    }
    df <- df[df$lower_age_limit < oldest$lower_age_limit, ]
    df <- rbind(df, oldest)
    rownames(df) <- seq_len(nrow(df))
  }

  ## aggregate by age groups
  df$lower_age_limit <- reduce_agegroups(
    df$lower_age_limit, union(0, lower_age_limits)
  )

  summarised <- by(df, list(df$lower_age_limit), function(x) {
    c(
      lower_age_limit = unique(x$lower_age_limit),
      immunity = mean(x$immunity)
    )
  })
  df <- as.data.frame(do.call(rbind, summarised))
  ret <- df$immunity
  names(ret) <- df$lower_age_limit
  return(ret)
}
