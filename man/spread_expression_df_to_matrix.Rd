% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataframe_builders.r
\name{spread_expression_df_to_matrix}
\alias{spread_expression_df_to_matrix}
\title{Convert a 'tidy' expression DF to a matrix}
\usage{
spread_expression_df_to_matrix(df, key, value)
}
\arguments{
\item{df}{A tidy data frame.}

\item{key}{Column name for the column containing gene names.}

\item{value}{Column name for the column containing sample IDs.}
}
\value{
x A matrix containing the expression results
}
\description{
Assumes the only other columns other than 'key' and 'val' are expression measurements
}
