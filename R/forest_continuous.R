#' Forest Plot for Continuous Outcomes
#'
#' This function generates a contour-enhanced forest plot for continuous outcomes, such as mean differences or standardized mean differences.
#'
#' @param dat A data frame containing the study data.
#' @param measure Character. Effect size measure: "SMD" (standardized mean difference) or "MD" (mean difference).
#' @param method Character. Meta-analysis model fitting method (e.g., "REML", "DL").
#' @param xlab Character. Label for the x-axis.
#' @param title Character. Main title of the plot.
#' @param model Character. Type of meta-analysis model: "Random-effects" or "Fixed-effect".
#' @param estimator Character. Estimator for between-study variance (tau²) in random-effects model.
#' @param m_c_col Character. Name of the column for control group means.
#' @param sd_c_col Character. Name of the column for control group standard deviations.
#' @param n_c_col Character. Name of the column for control group sample sizes.
#' @param m_t_col Character. Name of the column for treatment group means.
#' @param sd_t_col Character. Name of the column for treatment group standard deviations.
#' @param n_t_col Character. Name of the column for treatment group sample sizes.
#' @param diamond.col Character. Color of the summary effect diamond.
#' @param study.col Character. Color of study labels.
#' @param CI.col Character. Color of the confidence interval lines.
#' @param Pred.Inter.col Character. Color for the prediction interval lines.
#' @param square.size Numeric. Size of the study effect squares.
#' @param contour_fill Character. Fill color for contour-enhanced regions.
#' @param text_size Numeric. Base size of text elements in the plot.
#' @param pred Logical. If TRUE, include prediction interval.
#' @param xpos Numeric. Horizontal position for table columns or annotations.
#' @param xlim Numeric vector of length 2. Limits for the x-axis.
#' @param hetero_x Logical. If TRUE, display heterogeneity statistics (I², Q, tau²) on the plot.
#'
#' @return A ggplot object representing the forest plot.
#' @export
#'
#' @examples forest.continuous(dat1)
forest.continuous <- function(dat,
                              measure = "SMD",
                              method = "REML",
                              xlab = "",
                              title = NULL,
                              model = "Random-effects",
                              estimator = "REML",
                              m_c_col = "mean_c",
                              sd_c_col = "sd_c",
                              n_c_col = "n_c",
                              m_t_col = "mean_t",
                              sd_t_col = "sd_t",
                              n_t_col = "n_t",
                              diamond.col = "red",
                              study.col = "blue",
                              CI.col = "blue",
                              Pred.Inter.col = "black",
                              square.size = 10,
                              contour_fill = c("gray95","gray80","gray60","gray40"),
                              text_size = 3.5,
                              pred = TRUE,
                              xpos = NULL,
                              xlim = NULL,
                              hetero_x = NULL) {

  # -------------------------
  # Compute effect size
  # -------------------------
  escalc_res <- metafor::escalc(measure=measure,
                                m1i=dat[[m_t_col]], sd1i=dat[[sd_t_col]], n1i=dat[[n_t_col]],
                                m2i=dat[[m_c_col]], sd2i=dat[[sd_c_col]], n2i=dat[[n_c_col]],
                                data=dat)

  # Random-effects meta-analysis
  meta_res <- metafor::rma(yi, vi, data=escalc_res, method=method)
  pi <- predict(meta_res, digits=3)

  # -------------------------
  # Prepare data for plotting
  # -------------------------
  dat_plot <- escalc_res %>%
    mutate(weight = 100 * (1/vi)/sum(1/vi),
           ci.lb = yi - 1.96*sqrt(vi),
           ci.ub = yi + 1.96*sqrt(vi),
           y = rev(seq_len(nrow(.))),
           Effect = sprintf("%.2f [%.2f - %.2f]", yi, ci.lb, ci.ub),
           WeightText = sprintf("%.1f%%", weight),
           MeanT = sprintf("%.2f [%.2f]", dat[[m_t_col]], dat[[sd_t_col]]),
           MeanC = sprintf("%.2f [%.2f]", dat[[m_c_col]], dat[[sd_c_col]])) %>%
    arrange(desc(yi))

  # -------------------------
  # Pooled effect diamond
  # -------------------------
  pooled <- data.frame(
    x = c(meta_res$ci.lb, meta_res$beta, meta_res$ci.ub, meta_res$beta, meta_res$ci.lb),
    y = c(0, 0.5, 0, -0.5, 0)
  )
  pooled_text <- sprintf("%.2f [%.2f - %.2f]", meta_res$beta, meta_res$ci.lb, meta_res$ci.ub)

  # -------------------------
  # Heterogeneity
  # -------------------------
  tau_ci <- confint(meta_res)$random
  tau2 <- meta_res$tau2
  I2 <- meta_res$I2
  Q <- meta_res$QE
  df <- meta_res$k - 1
  I2_lb <- max(0, (Q/df - 1)/Q * 100)
  I2_ub <- min(100, ((Q * qchisq(0.975, df)/df) - 1)/Q * 100)
  hetero_text <- paste0(model, " (", estimator, ")\n",
                        "Tau² = ", round(tau2,3), " [", round(tau_ci[1],3), " - ", round(tau_ci[2],3), "]\n",
                        "I² = ", round(I2,1), "% [", round(I2_lb,1), " - ", round(I2_ub,1), "%]\n",
                        "Q-test p = ", signif(meta_res$QEp,3))

  # -------------------------
  # Default xpos
  # -------------------------
  effect_min <- pi$pi.lb
  effect_max <- pi$pi.ub

  if(is.null(xpos)){
    xpos <- list(
      Study = effect_min - 4,
      MeanT = effect_min - 3,
      MeanC = effect_min - 1.5,
      Effect = effect_max + 0.3,
      Weight = effect_max + 3,
      PredInt = effect_max + 0.9
    )
  }

  # -------------------------
  # Heterogeneity x-position
  # -------------------------
  if(is.null(hetero_x)){
    hetero_x <- min(xpos$Study, effect_min) - 0.5
  }

  # -------------------------
  # Default xlim
  # -------------------------
  if(is.null(xlim)){
    xlim <- c(min(xpos$Study-1, hetero_x-0.5), xpos$Weight + 2)
  }

  # -------------------------
  # Contour shading
  # -------------------------
  ymax_plot <- max(dat_plot$y) + 2
  ymin_plot <- -2

  contour_left <- data.frame(
    xmin = c(effect_min-0.2, -0.8, -0.5, -0.2),
    xmax = c(-0.8, -0.5, -0.2, 0),
    ymin = ymin_plot,
    ymax = ymax_plot,
    fill = factor(c("Very Large","Large","Medium","Small"),
                  levels = c("Small","Medium","Large","Very Large"))
  )
  contour_right <- data.frame(
    xmin = c(0, 0.2, 0.5, 0.8),
    xmax = c(0.2, 0.5, 0.8, effect_max+0.1),
    ymin = ymin_plot,
    ymax = ymax_plot,
    fill = factor(c("Small","Medium","Large","Very Large"),
                  levels = c("Small","Medium","Large","Very Large"))
  )
  contour_all <- rbind(contour_left, contour_right)

  pred_y <- -0.6

  # -------------------------
  # Build plot
  # -------------------------
  p <- ggplot() +
    geom_rect(data = contour_all, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill), alpha = 0.7) +
    scale_fill_manual(values = contour_fill, name = "", guide = guide_legend(override.aes = list(alpha = 1))) +
    geom_vline(xintercept = 0, linetype = "solid", color = "black") +
    geom_errorbarh(data = dat_plot, aes(y = y, xmin = ci.lb, xmax = ci.ub), height = 0, color = CI.col) +
    geom_point(data = dat_plot, aes(y = y, x = yi, size = weight), shape = 15, color = study.col, show.legend = FALSE) +
    geom_polygon(data = pooled, aes(x = x, y = y), fill = diamond.col, alpha = 0.9)

  if(pred){
    p <- p +
      geom_errorbarh(data = data.frame(y=pred_y, xmin=pi$pi.lb, xmax=pi$pi.ub),
                     aes(y=y, xmin=xmin, xmax=xmax),
                     height = 0, size = 1.5, color = Pred.Inter.col) +
      annotate("text", x = xpos$PredInt, y = pred_y,
               label = paste0("[", round(pi$pi.lb,2), "-", round(pi$pi.ub,2), "]"),
               fontface="bold", hjust=0, size=text_size) +
      annotate("text", x = xpos$Study+.9, y = pred_y,
               label = "Prediction Interval", fontface="bold", hjust=0, size=text_size)
  }

  # -------------------------
  # Left & Right tables
  # -------------------------
  p <- p +
    annotate("text", x=xpos$Study, y=dat_plot$y, label=dat_plot$Study, hjust=0, size=text_size, family="mono") +
    annotate("text", x=xpos$MeanT, y=dat_plot$y, label=dat_plot$MeanT, hjust=0, size=text_size, family="mono") +
    annotate("text", x=xpos$MeanC, y=dat_plot$y, label=dat_plot$MeanC, hjust=0, size=text_size, family="mono") +
    annotate("text", x=xpos$Effect, y=dat_plot$y, label=dat_plot$Effect, hjust=0, size=text_size, family="mono") +
    annotate("text", x=xpos$Weight, y=dat_plot$y, label=dat_plot$WeightText, hjust=0, size=text_size, family="mono") +

    # Headers
    annotate("text", x=xpos$Study, y=max(dat_plot$y)+1, label="Study", fontface="bold", hjust=0, size=text_size+.3) +
    annotate("text", x=xpos$MeanT, y=max(dat_plot$y)+1, label="Treatment \nmean[sd]", fontface="bold", hjust=0, size=text_size) +
    annotate("text", x=xpos$MeanC, y=max(dat_plot$y)+1, label="Control \nmean[sd]", fontface="bold", hjust=0, size=text_size) +
    annotate("text", x=xpos$Effect, y=max(dat_plot$y)+1, label=paste0(measure," (95% CI)"), fontface="bold", hjust=0, size=text_size) +
    annotate("text", x=xpos$Weight, y=max(dat_plot$y)+1, label="Weight (%)", fontface="bold", hjust=0, size=text_size) +

    # Pooled effect
    annotate("text", x=xpos$Study+.9, y=0.3, label="Pooled effect", fontface="bold", hjust=0, size=text_size) +
    annotate("text", x=xpos$Effect, y=0.4, label=pooled_text, fontface="bold", hjust=0, size=text_size) +

    # Heterogeneity
    annotate("text", x=hetero_x, y=-1.1, label=hetero_text, hjust=0, fontface="italic", size=text_size) +

    scale_y_continuous(breaks = dat_plot$y, labels = rep("", length(dat_plot$y))) +
    labs(x=xlab, y=NULL, title=ifelse(is.null(title), paste0("Contour Enhanced Forest Plot (", measure, ")"), title)) +
    theme_minimal(base_size=14, base_family="mono") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.y = element_blank(),
          legend.position = "top") +
    coord_cartesian(xlim=xlim, ylim=c(pred_y-0.5, max(dat_plot$y)+1.5))

  return(p)
}

