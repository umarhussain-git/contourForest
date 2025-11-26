#' Create a Continuous Outcome Forest Plot with Subgroups
#'
#' This function generates a forest plot for continuous outcomes, including study-specific means,
#' standard deviations, sample sizes, effect sizes (SMD or MD), subgroup analyses, and overall pooled effects.
#'
#' @param dat A data frame containing the study data.
#' @param m_t_col Column name for the treatment group means (string).
#' @param sd_t_col Column name for the treatment group standard deviations (string).
#' @param n_t_col Column name for the treatment group sample sizes (string).
#' @param m_c_col Column name for the control group means (string).
#' @param sd_c_col Column name for the control group standard deviations (string).
#' @param n_c_col Column name for the control group sample sizes (string).
#' @param subgroup_col Column name indicating subgroup membership (string).
#' @param study_col Column name for study labels (string).
#' @param measure Effect measure: "SMD" (Standardized Mean Difference) or "MD" (Mean Difference) (string).
#' @param method Meta-analysis method for pooling: e.g., "REML", "DL" (string).
#' @param xlab Label for the x-axis (string).
#' @param xlim Numeric vector of length 2 giving x-axis limits, e.g., c(min, max) (default NULL for auto).
#' @param model Model type: "random" or "fixed" effects (string).
#' @param title Plot title (string).
#' @param diamond.col Color for the pooled effect diamond (string).
#' @param overall.col Color for overall pooled effect text/line (string).
#' @param study.col.col Color for study labels (string).
#' @param CI.col Color for confidence intervals (string).
#' @param Pred.Inter.col Color for prediction interval (string).
#' @param square.size Size of the study effect squares (numeric).
#' @param Pred.Int.size Size of the prediction interval diamond (numeric).
#' @param text_size Font size for annotations (numeric).
#' @param pred Logical; if TRUE, include prediction interval (default FALSE).
#' @param xpos List of x-axis positions for annotations (e.g., Study, Weight, Effect) (list).
#'
#' @return A ggplot object of the forest plot.
#' @export
#'
#' @examples
#' \dontrun{
#' df <- data.frame(
#'   Study = c("Study A", "Study B"),
#'   subgroup = c("Group 1", "Group 2"),
#'   mean_t = c(7.4, 5.9),
#'   sd_t = c(2.0, 1.9),
#'   n_t = c(50, 60),
#'   mean_c = c(5.2, 4.8),
#'   sd_c = c(1.8, 2.0),
#'   n_c = c(50, 60)
#' )
#' xpos <- list(Study=0, MeanT=1, MeanC=2, Effect=3, Weight=4)
#' forest.continuous.subgroup(df, m_t_col="mean_t", sd_t_col="sd_t", n_t_col="n_t",
#'                            m_c_col="mean_c", sd_c_col="sd_c", n_c_col="n_c",
#'                            subgroup_col="subgroup", study_col="Study",
#'                            measure="MD", method="REML", xlab="Mean Difference",
#'                            xlim=NULL, model="random", title="Forest Plot",
#'                            diamond.col="blue", overall.col="red", study.col.col="black",
#'                            CI.col="black", Pred.Inter.col="grey", square.size=3,
#'                            Pred.Int.size=4, text_size=3, pred=TRUE, xpos=xpos)
#' }
forest.continuous.subgroup <- function(dat,
                                       m_t_col = "mean_t",
                                       sd_t_col = "sd_t",
                                       n_t_col = "n_t",
                                       m_c_col = "mean_c",
                                       sd_c_col = "sd_c",
                                       n_c_col = "n_c",
                                       subgroup_col = "subgroup",
                                       study_col = "Study",
                                       measure = "SMD",
                                       method = "REML",
                                       xlab = NULL,
                                       xlim = NULL,
                                       model = "RE",  # "RE" or "FE"
                                       title = "Subgroup Forest Plot",
                                       diamond.col = "red",
                                       overall.col = "darkgreen",
                                       study.col.col = "blue",
                                       CI.col = "blue",
                                       Pred.Inter.col = "black",
                                       square.size = 5,
                                       Pred.Int.size = 1.5,
                                       text_size = 3.5,
                                       pred = TRUE,
                                       xpos = list(
                                         Study   = NULL,
                                         MeanT   = NULL,
                                         MeanC   = NULL,
                                         Effect  = NULL,
                                         Weight  = NULL,
                                         PredInt = NULL,
                                         Hetero  = NULL)) {

  method <- ifelse(toupper(model) == "RE", "REML", "FE")
  if(is.null(xlab)) xlab <- ifelse(measure=="SMD","Effect size (SMD)","Mean Difference (MD)")

  # ---- Effect sizes ----
  escalc_res <- escalc(
    measure = measure,
    m1i = dat[[m_t_col]],
    sd1i = dat[[sd_t_col]],
    n1i = dat[[n_t_col]],
    m2i = dat[[m_c_col]],
    sd2i = dat[[sd_c_col]],
    n2i = dat[[n_c_col]],
    data = dat
  )

  escalc_res$weight <- 100 * (1/escalc_res$vi)/sum(1/escalc_res$vi)
  escalc_res$ci.lb <- escalc_res$yi - 1.96*sqrt(escalc_res$vi)
  escalc_res$ci.ub <- escalc_res$yi + 1.96*sqrt(escalc_res$vi)

  # ---- y positions ----
  escalc_res <- escalc_res %>% arrange(.data[[subgroup_col]])
  escalc_res$y <- NA
  y_counter <- 0
  subgroup_gap <- 2
  for(g in unique(escalc_res[[subgroup_col]])){
    idx <- which(escalc_res[[subgroup_col]]==g)
    n <- length(idx)
    escalc_res$y[idx] <- seq(from=y_counter+n, to=y_counter+1, by=-1)
    y_counter <- y_counter + n + subgroup_gap
  }

  # ---- Overall meta-analysis ----
  overall <- rma(yi, vi, data=escalc_res, method=method)

  # ---- Prediction interval ----
  if(pred){
    pi <- predict(overall, digits=3)
    effect_min <- min(min(escalc_res$ci.lb), pi$pi.lb)
    effect_max <- max(max(escalc_res$ci.ub), pi$pi.ub)
  } else {
    effect_min <- min(escalc_res$ci.lb)
    effect_max <- max(escalc_res$ci.ub)
  }

  # ---- x positions for labels ----
  xpos <- list(
    Study   = effect_min - 6,
    MeanT   = effect_min - 4,
    MeanC   = effect_min - 2.5,
    Effect  = effect_max + 0.5,
    Weight  = effect_max + 3.5,
    PredInt = effect_max + 0.5,
    Hetero  = effect_min - 6
  )


  # Example adjustments
  xlim <- c(xpos$Study - 1, xpos$Weight + 2)

  # ---- Contour shading ----
  ymin_plot <- -3
  ymax_plot <- max(escalc_res$y)+2
  contour_fill <- c("Small"="gray95","Medium"="gray85","Large"="gray70","Very Large"="gray50")
  contour_all <- data.frame(
    xmin = c(effect_min-0.1, -0.8, -0.5, -0.2, 0, 0.2, 0.5, 0.8),
    xmax = c(-0.8, -0.5, -0.2, 0, 0.2, 0.5, 0.8, effect_max+0.1),
    ymin = ymin_plot,
    ymax = ymax_plot,
    fill = factor(c("Very Large","Large","Medium","Small","Small","Medium","Large","Very Large"),
                  levels=c("Small","Medium","Large","Very Large"))
  )

  # ---- Subgroup pooled effects & heterogeneity ----
  subgroup_effects <- lapply(unique(escalc_res[[subgroup_col]]), function(g){
    subset_g <- escalc_res %>% filter(.data[[subgroup_col]]==g)
    res <- rma(yi, vi, data=subset_g, method=method)
    tau2_ci <- confint(res, parm="tau2")$random
    I2 <- res$I2
    I2_lower <- max(0, I2 - 1.96*sqrt(res$vi/length(subset_g))) # crude approx
    I2_upper <- min(100, I2 + 1.96*sqrt(res$vi/length(subset_g)))
    data.frame(
      subgroup=g,
      b=res$b,
      ci.lb=res$ci.lb,
      ci.ub=res$ci.ub,
      y=min(subset_g$y)-0.7,
      tau2=res$tau2,
      tau2.lb=tau2_ci[1],
      tau2.ub=tau2_ci[2],
      I2=I2,
      I2.lb=I2_lower,
      I2.ub=I2_upper,
      pval=res$QEp
    )
  }) %>% bind_rows()

  # ---- Overall diamond ----
  y_overall <- -1.5
  overall_diamond <- data.frame(
    x=c(overall$ci.lb, overall$b, overall$ci.ub, overall$b, overall$ci.lb),
    y=c(y_overall, y_overall+0.4, y_overall, y_overall-0.4, y_overall)
  )

  # ---- Base plot ----
  p <- ggplot(escalc_res, aes(x=yi, y=y)) +
    geom_rect(data=contour_all, aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=fill),
              alpha=0.6, inherit.aes=FALSE) +
    scale_fill_manual(values=contour_fill, name=" ") +
    geom_segment(aes(x=0,xend=0,y=ymin_plot,yend=ymax_plot),
                 color="gray10", size=0.4, linetype="solid", inherit.aes=FALSE) +
    theme_minimal() +
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.title.x=element_text(face="bold"),
          axis.title.y=element_blank(),
          axis.text.x=element_text(),
          axis.text.y=element_blank(),
          legend.position="top") +
    xlab(xlab) +
    geom_errorbarh(aes(xmin=ci.lb,xmax=ci.ub), height=0, color=CI.col) +
    geom_point(aes(size=weight), shape=15, color=study.col.col) +
    scale_size_continuous(range=c(2,square.size), guide="none") +
    geom_polygon(data=overall_diamond, aes(x=x,y=y),
                 fill=overall.col, alpha=0.9, inherit.aes=FALSE) +
    # headers
    annotate("text", x=xpos$Study, y=ymax_plot, label="Study", fontface="bold", hjust=0, size=text_size) +
    annotate("text", x=xpos$MeanT, y=ymax_plot, label="Treatment\nMean[SD]", fontface="bold", hjust=0, size=text_size) +
    annotate("text", x=xpos$MeanC, y=ymax_plot, label="Control\nMean[SD]", fontface="bold", hjust=0, size=text_size) +
    annotate("text", x = xpos$Effect, y = ymax_plot,
             label = paste0(measure, " [95% CI]"),   fontface = "bold", hjust = 0, size = text_size) +

    annotate("text", x=xpos$Weight, y=ymax_plot, label="Weight (%)", fontface="bold", hjust=0, size=text_size)

  # ---- Study-level values ----
  p <- p +
    geom_text(aes(x=xpos$Study,label=.data[[study_col]]), hjust=0, size=text_size) +
    geom_text(aes(x=xpos$MeanT,label=sprintf("%.2f [%.2f]", .data[[m_t_col]], .data[[sd_t_col]])),
              hjust=0, size=text_size) +
    geom_text(aes(x=xpos$MeanC,label=sprintf("%.2f [%.2f]", .data[[m_c_col]], .data[[sd_c_col]])),
              hjust=0, size=text_size) +
    geom_text(aes(x=xpos$Effect,label=sprintf("%.2f [%.2f - %.2f]", yi, ci.lb, ci.ub)),
              hjust=0, size=text_size) +
    geom_text(aes(x=xpos$Weight,label=sprintf("%.1f%%", weight)), hjust=0, size=text_size)

  # ---- Subgroup diamonds + heterogeneity labels ----
  for(i in seq_len(nrow(subgroup_effects))){
    y_d <- subgroup_effects$y[i]
    poly_df <- data.frame(
      x = c(subgroup_effects$ci.lb[i], subgroup_effects$b[i],
            subgroup_effects$ci.ub[i], subgroup_effects$b[i],
            subgroup_effects$ci.lb[i]),
      y = c(y_d, y_d+0.4, y_d, y_d-0.4, y_d)
    )

    p <- p +
      geom_polygon(data=poly_df, aes(x=x, y=y),
                   fill=diamond.col, alpha=0.8, inherit.aes=FALSE) +

      # Left: subgroup label
      annotate("text", x = xpos$Study+ 1, y = y_d - 0.2,
               label = subgroup_effects$subgroup[i],
               hjust = 1, fontface = "bold", size = text_size) +

      # Right: numeric effect size without "MD"
      annotate("text", x = effect_max + 0.3, y = y_d,
               label = sprintf("%.2f [%.2f - %.2f]",
                               subgroup_effects$b[i],      # point estimate (MD)
                               subgroup_effects$ci.lb[i],  # lower 95% CI
                               subgroup_effects$ci.ub[i]), # upper 95% CI
               hjust = 0, fontface = "bold", size = text_size) +

      # Left: heterogeneity values
      annotate("text", x = effect_min - 1.5, y = y_d - 0.5,
               label = sprintf("Hetero:I²=%.1f%% [%.1f%%-%.1f%%]\nτ²=%.2f [%.2f-%.2f], p=%s",
                               subgroup_effects$I2[i],
                               subgroup_effects$I2.lb[i], subgroup_effects$I2.ub[i],
                               subgroup_effects$tau2[i],
                               subgroup_effects$tau2.lb[i], subgroup_effects$tau2.ub[i],
                               ifelse(subgroup_effects$pval[i] < 0.001, "<0.001",
                                      sprintf("%.3f", subgroup_effects$pval[i]))),
               hjust = 1, fontface = "italic", color = "black", size = text_size-.2)
  }

  # ---- Overall label ----
  p <- p +
    annotate("text",
             x = xpos$Study,  y = y_overall + 0.2, label = paste0("Pooled Effect (", model, ")"),
             hjust = 0,    size = text_size,       fontface = "bold")+
    annotate("text", x=xpos$Effect, y=y_overall+0.2,
             label=sprintf("%.2f [%.2f - %.2f]", overall$b, overall$ci.lb, overall$ci.ub),
             hjust=0, size=text_size, fontface="bold")+
    annotate("text",
             x = xpos$Study,  y = y_overall + 0.2, label = paste0("Pooled Effect (", model, ")"),
             hjust = 0,    size = text_size,       fontface = "bold")+
    annotate("text", x=xpos$Effect, y=y_overall+0.2,
             label=sprintf("%.2f [%.2f - %.2f]", overall$b, overall$ci.lb, overall$ci.ub),
             hjust=0, size=text_size, fontface="bold")+
    ## Label for Prediction Interval below the pooled effect
    annotate("text",
             x = xpos$Study,
             y = y_overall - 0.4,   # slightly below the pooled effect
             label = "Prediction Interval",
             hjust = 0, size = text_size, fontface = "bold", color = "black")

  # ---- Prediction interval ----
  if(pred){
    pi_y <- y_overall - 0.7
    p <- p +
      geom_errorbarh(aes(xmin=pi$pi.lb, xmax=pi$pi.ub, y=pi_y),
                     height=0, size=Pred.Int.size, color=Pred.Inter.col, inherit.aes=FALSE) +
      annotate("text", x=xpos$PredInt, y=pi_y,
               label=sprintf("[%.2f - %.2f]", pi$pi.lb, pi$pi.ub),
               hjust=0, fontface="bold", size=text_size)
  }

  if(!is.null(xlim)){
    p <- p + scale_x_continuous(limits=xlim)
  } else {
    p <- p + scale_x_continuous(limits=c(effect_min-0.5, effect_max+0.5))
  }
  # ---- Overall heterogeneity annotation ----
    p <- p +
    annotate("text",
             x = xpos$Study,   # left side of plot
             y = y_overall - 1.3,    # slightly below pooled effect/PI
             label = sprintf("Overall Hetero: I²=%.1f%%, τ²=%.2f,\n p=%s",
                             overall$I2,
                             overall$tau2,
                             ifelse(overall$QEp < 0.001, "<0.001", sprintf("%.3f", overall$QEp))),
             hjust = 0, fontface = "italic", size = text_size, color = "black")


  return(p)
}
