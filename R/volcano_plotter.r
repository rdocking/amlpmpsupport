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
