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
