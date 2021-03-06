% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility_functions.r
\name{nest_libraries_by_patient}
\alias{nest_libraries_by_patient}
\title{Add nested library information to a patient-level data frame}
\usage{
nest_libraries_by_patient(patients.df, libraries.df, key = "tfl_id")
}
\arguments{
\item{patients.df}{Data frame containing patient-level data.}

\item{libraries.df}{Data frame containing library-level data.}

\item{key}{Key to join the two data-frames by (default is 'tfl_id').}
}
\value{
patients.w_libs.df A data frame patients with nested library information
}
\description{
Add nested library information to a patient-level data frame
}
