#' Generate a dot-plot of expression of selected genes, facetted by selected feature
#'
#' @param expression_df Data frame containing expression values
#' @param feature_label Label in the data frame to facet by
#'
#' @return p A ggplot object containing the generated plot
#' @export
#'
hit_plotter <- function(expression_df, feature_label){
  p <- ggplot2::ggplot(expression_df)
  p <- p + ggplot2::aes_string(x=feature_label, colour=feature_label)
  p <- p + ggplot2::aes(y=TPM)
  p <- p + ggplot2::geom_jitter(size=4) + ggplot2::scale_y_log10() + ggplot2::scale_colour_brewer(palette='Set1')
  p <- p + ggplot2::theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  p <- p + ggplot2::facet_wrap(~ gene) + ggplot2::theme(legend.position="none")
  return(p)
}

#' Generate a density-plot of expression of selected genes, facetted by selected feature
#'
#' @param expression_df Data frame containing expression values
#' @param feature_label Label in the data frame to facet by
#'
#' @return p A ggplot object containing the generated plot
#' @export
#'
density_plotter <- function(expression_df, feature_label){
  p <- ggplot2::ggplot(expression_df)
  p <- p + ggplot2::aes_string(fill=feature_label)
  p <- p + ggplot2::aes(x=TPM) + ggplot2::scale_x_log10()
  p <- p + ggplot2::geom_density(alpha=0.5) + ggplot2::scale_fill_brewer(palette='Set1')
  p <- p + ggplot2::facet_wrap(~ gene, scales = "free_y")
  return(p)
}

#' Generate a plot of enriched pathways
#'
#' @param gage_hits A data frame produced by the gage function
#' @param plot_limit The number of top pathways to plot
#'
#' @return p A ggplot object containing the generated plot
#' @export
#'
pathway_enrichment_plot <- function(gage_hits, plot_limit = 20){
  # Subset the data frame to the top plot_limit hits
  plot_subset <- head(gage_hits, plot_limit)
  # Arrange by q-value - reverse because of the coordinate flip below
  plot_subset$pathway <- factor(plot_subset$pathway,
                                levels = rev(plot_subset$pathway))
  p <- ggplot2::ggplot(plot_subset,
                       ggplot2::aes(x=pathway, y=q.val, size=set.size, colour=stat.mean))
  p <- p + ggplot2::geom_point() + ggplot2::scale_y_log10() + ggplot2::coord_flip()
  return(p)
}

#' Generate a volcano plot of expression hits
#'
#' @param voom_hits Data frame containing voom differential expression results
#' @param label_threshold Threshold for -log10(PVal) to add a label
#'
#' @return p A ggplot object containing the generated plot
#' @export
#'
volcano_plotter <- function(voom_hits, label_threshold = 30){
  # Use 'negative log10p' to get a Phred-like scale
  voom_hits$neg_log10_p <- log10(voom_hits$adj.P.Val) * -1
  # Add labels for top N hits
  voom_hits$label <- ifelse(voom_hits$neg_log10_p >= label_threshold,
                            voom_hits$gene,
                            NA)
  # Set up the plot
  p <- ggplot2::ggplot(voom_hits, aes(x=logFC, y=neg_log10_p,
                                      colour=neg_log10_p, label=label))
  p <- p + ggplot2::geom_point(alpha=0.6)
  p <- p + ggrepel::geom_label_repel(box.padding = unit(0.5, "lines"),
                                     na.rm=TRUE)
  return(p)
}

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

#' Generate a grid of plots with a shared legend
#'
#' Adapted from:
#'   https://github.com/tidyverse/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
#'   http://rpubs.com/sjackman/grid_arrange_shared_legend
#'   http://stackoverflow.com/questions/13649473/add-a-common-legend-for-combined-ggplots
#'
#' @param ... A set of plots to arrange with a shared legend
#'
#' @return A ggplot plot
#' @export
#'
grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  combined <- switch(position,
                     "bottom" = gridExtra::arrangeGrob(do.call(gridExtra::arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = grid::unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = gridExtra::arrangeGrob(do.call(gridExtra::arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = grid::unit.c(unit(1, "npc") - lwidth, lwidth)))
  grid::grid.newpage()
  grid::grid.draw(combined)
  # return gtable invisibly
  invisible(combined)
}

