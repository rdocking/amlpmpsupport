#' Generate a prinicple component plot for a subset of genes, coloured by a given feature
#'
#' @param gene_set A vector containing gene names to subset
#' @param exp.design Experimental design data frame
#' @param feature_label Label in the data frame to facet by
#' @param counts.mat Matrix containing expression count data
#'
#' @return p A ggplot object containing the generated plot
#' @export
#'
subset_pc_plot <- function(gene_set, exp.design, feature_label, counts.mat){
  # Find rows to keep
  rows_to_keep <- rownames(counts.mat) %in% gene_set
  # Subset the count matrix to the genes of interest
  counts.subset.mat <- counts.mat[rows_to_keep,]
  # Calculate the principle components
  components <- prcomp(counts.subset.mat, center = F, scale = F)
  rotations <- as.data.frame(components$rotation)
  # Extract the first two components, and add labels from the design data frame
  pc1_pc2 <- subset(rotations, select = c("PC1", "PC2"))
  pc1_pc2$rna_seq_lib <- rownames(pc1_pc2)
  pc1_pc2_w_design <- dplyr::left_join(pc1_pc2, exp.design, by='rna_seq_lib')
  p <- ggplot2::ggplot(pc1_pc2_w_design, aes(x=PC1, y=PC2))
  p <- p + ggplot2::aes_string(colour=feature_label)
  # p <- p + geom_label_repel(box.padding = unit(0.5, "lines"), na.rm=TRUE)
  p <- p + ggplot2::geom_point(size = 5) + ggplot2::scale_colour_brewer(palette = "Set1")
  return(p)
}
