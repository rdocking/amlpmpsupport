#' Fit a Cox proportional hazards model
#'
#' @param exp.design Experimental design data frame
#' @param label_col Column in the design data frame used to separate samples
#'
#' @return cx The cox proportional hazards model
#' @export
#'
fit_coxph <- function(exp.design, label_col){
  # Fit Cox proportional hazards model
  cx <- survival::coxph(survival::Surv(survival_time, survival_event) ~ label_col,
                        data=exp.design)
  return(cx)
}

#' Plot survival curves
#'
#' @param exp.design Experimental design data frame
#' @param label_col Column in the design data frame used to separate samples
#'
#' @return p A ggplot object containing the generated plot
#' @export
#'
plot_survival <- function(exp.design, label_col){
  # Generate a survival model using the survfit function
  surv <- survival::survfit(survival::Surv(survival_time, survival_event) ~ label_col,
                            data = exp.design)
  p <- GGally::ggsurv(surv, CI = TRUE) + ggplot2::ggtitle('Kaplan-Meier Estimator')
  p
}
