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
