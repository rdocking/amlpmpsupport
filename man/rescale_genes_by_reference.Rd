% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cna_workflow.r
\name{rescale_genes_by_reference}
\alias{rescale_genes_by_reference}
\title{Rescale all genes by reference expression levels}
\usage{
rescale_genes_by_reference(expr.mat, reference_means)
}
\arguments{
\item{expr.mat}{a matrix of expression values}

\item{reference_means}{a named vector (With values for all genes, which is subset to the same genes retained in the matrix)}
}
\value{
expr.mat
}
\description{
Rescale all genes by reference expression levels
}
