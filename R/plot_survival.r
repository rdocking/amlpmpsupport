#' Plot predictions from a randomForest model
#'
#' @param exp.design Experimental design data frame
#' @param label_col Column in the design data frame used to separate samples
#'
#' @return p A ggplot object containing the generated plot
#' @export
#'
plot_survival <- function(exp.design, label_col){
  # Generate a survival model using the survfit function
  surv <- survival::survfit(Surv(survival_time, survival_event) ~ label_col,
                            data = exp.design)
  p <- GGally::ggsurv(surv, CI = TRUE) + ggplot2::ggtitle('Kaplan-Meier Estimator')
  p
}
