% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/project_immunity.r
\name{project_immunity}
\alias{project_immunity}
\title{Project immunity from a baseline}
\usage{
project_immunity(baseline.immunity, baseline.year, year, coverage,
  schedule, maternal.immunity, efficacy)
}
\arguments{
\item{baseline.immunity}{baseline immunity, as a named vector; the names correspond to lower limits of the age groups, and the vector itself to the corresponding levels of immunity.}

\item{baseline.year}{year at which baseline immunity is taken (corresponding to a column in the \code{coverage} argument)}

\item{year}{year to project to}

\item{coverage}{coverage with multiple vaccine doses, given as a matrix in which each row is a dose and each (named) column a year}

\item{schedule}{the ages at which vaccines are given (in years).}

\item{maternal.immunity}{the proportion maternally immune.}

\item{efficacy}{vaccine efficacy.}
}
\value{
a data frame of immunity levels by age group (as in \code{baseline.immunity}).
}
\description{
Project immunity from a baseline via vaccination coverage rates
}
\author{
Sebastian Funk <sebastian.funk@lshtm.ac.uk>
}