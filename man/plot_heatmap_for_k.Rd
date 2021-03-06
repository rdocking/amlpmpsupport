% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting_functions.r
\name{plot_heatmap_for_k}
\alias{plot_heatmap_for_k}
\title{Plot a pheatmap for a single value of k}
\usage{
plot_heatmap_for_k(
  cluster_assignments.df,
  annotation_col.df,
  ann_colours.lst,
  cluster_input.mat,
  annotation_legend = FALSE
)
}
\arguments{
\item{cluster_assignments.df}{A data frame containing cluster assignments for samples}

\item{annotation_col.df}{A data frame containing sample annotations}

\item{ann_colours.lst}{A named list containing annotation colours}

\item{cluster_input.mat}{A matrix of expression values for the heatmap body.}

\item{annotation_legend}{Boolean - whether or not to include an annotation legend.}
}
\description{
Plot a pheatmap for a single value of k
}
