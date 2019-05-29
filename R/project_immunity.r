project_immunity <- function(baseline.immunity, coverage, schedule)
{
  ## checks
  if (missing(baseline.immunity)) stop("'coverage' must be provided")
  if (missing(coverage)) stop("'coverage' must be provided")
  if (length(unique(vapply(coverage, length, 0L))) > 1)
  {
    stop("All elements of 'coverage' must be the same length")
  }

  ## convert baseline to annual immunity
  bdf <- data.frame(lower.age.limit=as.integer(names(baseline.immunity)),
                    immunity=baseline.immunity)
  min_age <- min(bdf$lower.age.limit)
  max_age <- max(bdf$lower.age.limit)
  df <- data.frame(lower.age.limit=
                     do.call(seq, as.list(range(bdf$lower.age.limit))))
  df[df$lower.age.limit %in% bdf$lower.age.limit, "immunity"] <- bdf$immunity
  ## fill
  df$immunity <- c(NA, na.omit(df$immunity))[cumsum(!is.na(df$immunity))+1]

  for (i in seq_len(length(coverage[[1]]))) {
    ## move all one age group up
    df[seq(2, nrow(df) - 1), "immunity"] <- df[seq(1, nrow(df) - 2), "immunity"]
    ## implement vaccination schedule
    df[df$lower.age.limit==1, "immunity"] <- 0
  }
}
