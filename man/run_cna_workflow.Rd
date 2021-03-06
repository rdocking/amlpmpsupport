% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cna_workflow.r
\name{run_cna_workflow}
\alias{run_cna_workflow}
\title{Run the CNA workflow}
\usage{
run_cna_workflow(
  expr.mat,
  genes.bed,
  window_size = NA,
  mean_tpm_threshold = 0,
  cv_threshold = 0,
  reference_labels,
  scale_samples_by_zero = TRUE
)
}
\arguments{
\item{expr.mat}{A matrix of gene-expression values}

\item{genes.bed}{A bed-like file of gene coordinates}

\item{window_size}{Window size to use for rolling mean}

\item{mean_tpm_threshold}{TPM threshold to retain genes}

\item{cv_threshold}{CV threshold to retain genes}

\item{reference_labels}{Reference labels}

\item{scale_samples_by_zero}{Scale samples by zero?}
}
\value{
A re-scaled matrix
}
\description{
Run the CNA workflow
}
