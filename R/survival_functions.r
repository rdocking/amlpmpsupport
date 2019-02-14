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

#' Modified survival plot
#'
#' @param surv.fit A survival object from survival::Surv
#' @param data The data frame used to generate the survival object
#' @param palette A named list containing plotting colours
#' @param title Plot title
#' @param ... Additional parameters to pass to ggsurvplot
#'
#' @return
#' @export
surv_panel_plot <- function(surv.fit, data, palette, title, pval.coord = c(600, 0.9), xlim = c(0, 800), ...){
  p <- survminer::ggsurvplot(surv.fit,
                  data = data,
                  # conf.int = TRUE,        # Add confidence interval
                  pval = TRUE,              # Add p-value
                  pval.coord = pval.coord,
                  pval.method = FALSE,
                  risk.table = TRUE,        # Add risk table
                  risk.table.col = "strata",  # Risk table color by groups)
                  risk.table.y.text = FALSE,
                  risk.table.height = 0.3,
                  risk.table.title = "",
                  ggtheme = theme_bw(),
                  ncensor.plot = FALSE,
                  legend = "none",
                  palette = palette,
                  legend.title = "",
                  xlab = "Time (Weeks)",
                  xlim = xlim,
                  break.x.by = 100,
                  fontsize = c(10),       # fontsize for number-at-risk table
                  font.legend = c(16),
                  font.main = c(20),
                  font.x = c(16),
                  font.y = c(16),
                  font.tickslab = c(16),
                  title = title,
                  surv.median.line = "none",
                  ...)

  # Hack into the object to tweak further
  p[["table"]][["labels"]][["x"]] <- ""
  p[["table"]][["labels"]][["y"]] <- ""
  p
}
