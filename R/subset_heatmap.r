#' Generate a heatmap for a subset of genes
#'
#' @param gene_set A vector containing gene names to plot
#' @param annotation_col Data frame containing annotation columns for pheatmap
#' @param counts.mat A matrix of gene count data
#'
#' @return p A ggplot plot
#' @export
#'
subset_heatmap <- function(gene_set, annotation_col, counts.mat){
  # Find rows to keep
  rows_to_keep <- rownames(counts.mat) %in% gene_set
  # Subset the count matrix to the genes of interest
  counts.subset.mat <- counts.mat[rows_to_keep,]
  # log2 transform the data for plotting
  counts.subset.mat <- log2(counts.subset.mat + 1)
  # Generate the heatmap
  p <- pheatmap::pheatmap(counts.subset.mat,
                          annotation_col = annotation_col,
                          show_rownames = FALSE,
                          show_colnames = FALSE,
                          fontsize_col = 6)
  return(p)
}
