#' Retrieve cluster assignments for a given value of k
#'
#' @param k Number of clusters
#' @param cluster_obj Clustering object returned from CCP clustering
#'
#' @return cluster_assignments.df A dataframe containing cluster assignments for each sample
#' @export
get_cluster_assignments_for_k <- function(k, cluster_obj){
  # Extract the cluster membership from the clustering object
  cluster_for_k <- cluster_obj[[k]]
  cluster_assignments <-
    as.matrix(cluster_for_k$consensusClass[cluster_for_k$consensusTree[["order"]]])
  colnames(cluster_assignments) = "cluster"
  # Convert cluster assignments to a data frame
  cluster_assignments.df <-
    cluster_assignments %>%
    as.data.frame() %>%
    rownames_to_column('library_name')
  # Return the cluster assignments data frame
  return(cluster_assignments.df)
}
