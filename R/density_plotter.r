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
