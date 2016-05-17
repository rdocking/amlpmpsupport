#' Plot predictions from a randomForest model
#'
#' @param rf_predictions A vector of predictions from a randomForest model
#'
#' @return p A ggplot plot
#' @export
#'
plot_rf_predictions <- function(rf_predictions){
  # Convert to long for plotting
  rf_predictions.long <- tidyr::gather(rf_predictions,
                                       key=predicted_class, probability,
                                       -rna_seq_lib, -prediction, -prior_label)
  # Generate the plot
  p <- ggplot2::ggplot(rf_predictions.long,
                       aes(x=factor(predicted_class),
                           y=probability,
                           group=rna_seq_lib,
                           colour=prior_label,
                           label=rna_seq_lib))
  p <- p + ggplot2::geom_point() + ggplot2::geom_line()
  p <- p + ggplot2::scale_colour_brewer(palette='Set1')
  return(p)
}
