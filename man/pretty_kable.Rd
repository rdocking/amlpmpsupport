% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility_functions.r
\name{pretty_kable}
\alias{pretty_kable}
\title{Set a function for pretty-printing kables}
\usage{
pretty_kable(df, ...)
}
\arguments{
\item{df}{Data frame containing data to format as a table}

\item{...}{Additional options to be passed to kable}
}
\value{
A kable table
}
\description{
Set a function for pretty-printing kables
}
\examples{
pretty_kable(data.frame(x = c(1000, 999), y = c(0.345:0.351)))
}
