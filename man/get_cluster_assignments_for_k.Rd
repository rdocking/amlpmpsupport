% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clustering_helpers.r
\name{get_cluster_assignments_for_k}
\alias{get_cluster_assignments_for_k}
\title{Retrieve cluster assignments for a given value of k}
\usage{
get_cluster_assignments_for_k(k, cluster_obj)
}
\arguments{
\item{k}{Number of clusters}

\item{cluster_obj}{Clustering object returned from CCP clustering}
}
\value{
cluster_assignments.df A dataframe containing cluster assignments for each sample
}
\description{
Retrieve cluster assignments for a given value of k
}
