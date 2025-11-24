#' Contour-Enhanced Forest Plot for Binary Outcomes
#'
#' Creates a contour-enhanced forest plot for binary outcomes (e.g., risk ratios or odds ratios)
#' with subgroup-specific meta-analyses, overall effect, heterogeneity statistics, and optional
#' prediction intervals.
#'
#' @param dat A data frame containing the study-level data including events, sample sizes, and study labels.
#' @param subgroup Optional character vector specifying the subgroup for each study.
#' @param measure Character specifying the effect size measure (e.g., "RR", "OR").
#' @param method Method for meta-analysis (e.g., "REML", "DL").
#' @param xlab Label for the x-axis.
#' @param title Main title of the forest plot.
#' @param diamond.col Color for the overall effect diamond.
#' @param overall.col Color for the overall effect line.
#' @param study.col Color for study labels.
#' @param CI.col Color for confidence intervals.
#' @param Pred.Inter.col Color for prediction intervals (if displayed).
#' @param square.size Size of the squares representing individual study weights.
#' @param Pred.Int.size Line width of prediction interval.
#' @param xlim Numeric vector of length 2 specifying the x-axis limits.
#' @param text_size Font size for text annotations.
#' @param xpos Numeric vector specifying x positions for annotations.
#' @param study_x Numeric value for the x-position of study labels.
#' @param val_x Numeric value for the x-position of effect estimates.
#' @param pred Logical; if TRUE, prediction intervals are displayed.
#' @importFrom dplyr group_modify mutate arrange summarise left_join
#' @import dplyr
#'
#' @return A `ggplot` object representing the contour-enhanced forest plot.
#' @export
#'
#' @examples forest.binary.subgroup(dat)
forest.binary.subgroup <- function(dat,
                                   subgroup = NULL,
                                   measure = "RR",
                                   method = "REML",
                                   xlab = "Risk Ratio (RR)",
                                   title = "Subgroup Forest Plot",
                                   diamond.col="red",
                                   overall.col="darkgreen",
                                   study.col="blue",
                                   CI.col="blue",
                                   Pred.Inter.col="black",
                                   square.size=8,
                                   Pred.Int.size = 2,
                                   xlim = c(-2, 3.5),
                                   text_size = 3.5,
                                   xpos = list(
                                     EventsT = -1.2,
                                     EventsC = -0.4,
                                     Effect = 2.5,
                                     Weight = 3.2
                                   ),
                                   study_x = -1.5,
                                   val_x = 2.6,
                                   pred=TRUE) {

  # ---- Individual study effect sizes ----
  dat <- dat %>%
    mutate(
      yi = log((events_t / n_t)/(events_c / n_c)),
      sei = sqrt(1/events_t - 1/n_t + 1/events_c - 1/n_c),
      weight = 100 * (1/sei^2) / sum(1/sei^2),
      ci.lb = yi - 1.96*sei,
      ci.ub = yi + 1.96*sei
    )

  # ---- Subgroup meta-analysis ----
  sub_results <- dat %>%
    group_by(subgroup) %>%
    group_modify(~{
      m <- rma(ai=.x$events_t, n1i=.x$n_t,
               ci=.x$events_c, n2i=.x$n_c,
               data=.x, measure=measure, method=method)
      tau_ci <- confint(m, parm="tau2")$random
      I2.se <- sqrt(2/(m$k - 1)) * 100

      tibble(
        subgroup = unique(.x$subgroup),
        yi = round(m$b,2),
        ci.lb = round(m$ci.lb,2),
        ci.ub = round(m$ci.ub,2),
        tau2 = round(m$tau2,2),
        tau2.lb = round(tau_ci[1],2),
        tau2.ub = round(tau_ci[2],2),
        I2 = round(m$I2,2),
        I2.lb = round(max(0, m$I2 - 1.96*I2.se),2),
        I2.ub = round(min(100, m$I2 + 1.96*I2.se),2),
        pval = round(m$QEp,3)
      )
    })

  # ---- Overall random-effects model ----
  overall <- rma(ai=events_t, n1i=n_t, ci=events_c, n2i=n_c, data=dat, measure=measure, method=method)
  overall_tau_ci <- confint(overall, parm="tau2")$random
  I2.se_overall <- sqrt(2/(overall$k - 1)) * 100
  overall_I2_lb <- round(max(0, overall$I2 - 1.96*I2.se_overall),2)
  overall_I2_ub <- round(min(100, overall$I2 + 1.96*I2.se_overall),2)

  # ---- Prediction interval ----
  if(pred){
    pi <- predict(overall, digits=3, transf=exp)
  }

  # ---- y-axis positions ----
  dat <- dat %>% arrange(subgroup)
  dat$y <- NA
  y_counter <- 0
  subgroup_gap <- 2
  for (g in unique(dat$subgroup)) {
    idx <- which(dat$subgroup == g)
    n <- length(idx)
    dat$y[idx] <- seq(from = y_counter + n, to = y_counter + 1, by = -1)
    y_counter <- y_counter + n + subgroup_gap
  }

  # ---- Subgroup diamonds & heterogeneity labels ----
  diamonds <- list()
  hetero <- data.frame()
  for (g in unique(dat$subgroup)) {
    subset_g <- dat %>% filter(subgroup == g)
    res <- sub_results %>% filter(subgroup == g)
    y_diamond <- min(subset_g$y) - 0.7

    diamonds[[g]] <- data.frame(
      x = c(exp(res$ci.lb), exp(res$yi), exp(res$ci.ub), exp(res$yi), exp(res$ci.lb)),
      y = c(y_diamond, y_diamond + 0.4, y_diamond, y_diamond - 0.4, y_diamond)
    )

    hetero <- rbind(hetero, data.frame(
      subgroup = g,
      y = y_diamond - 0.8,
      label = paste0("I² = ", res$I2, "% [", res$I2.lb, "-", res$I2.ub,
                     "], τ² = ", res$tau2, " [", res$tau2.lb, "-", res$tau2.ub,
                     "], p = ", res$pval)
    ))
  }

  # ---- Overall diamond & labels ----
  y_overall <- -1.5
  overall_diamond <- data.frame(
    x = c(exp(overall$ci.lb), exp(round(overall$b,2)), exp(overall$ci.ub), exp(round(overall$b,2)), exp(overall$ci.lb)),
    y = c(y_overall, y_overall + 0.4, y_overall, y_overall - 0.4, y_overall)
  )
  overall_label <- paste0("I² = ", round(overall$I2,2), "% [", overall_I2_lb, "-", overall_I2_ub,
                          "], τ² = ", round(overall$tau2,2), " [", round(overall_tau_ci[1],2), "-", round(overall_tau_ci[2],2), "]")
  model_label <- paste0("Random effects (", method, ")")

  # ---- Contours ----
  ymax_plot <- max(dat$y) + 2
  ymin_plot <- -3
  contour_left <- data.frame(
    xmin = c(0, 1/2.0, 1/1.5, 1/1.2),
    xmax = c(1/2.0, 1/1.5, 1/1.2, 1),
    ymin = ymin_plot,
    ymax = ymax_plot,
    fill = factor(c("Very Large","Large","Medium","Small"), levels=c("Small","Medium","Large","Very Large"))
  )
  contour_right <- data.frame(
    xmin = c(1, 1.2, 1.5, 2.0),
    xmax = c(1.2, 1.5, 2.0, 2.5),
    ymin = ymin_plot,
    ymax = ymax_plot,
    fill = factor(c("Small","Medium","Large","Very Large"), levels=c("Small","Medium","Large","Very Large"))
  )
  contour_all <- rbind(contour_left, contour_right)

  # ---- Labels for subgroup headers ----
  group_labels <- dat %>% group_by(subgroup) %>% summarise(y_label = max(y) + 0.8)
  subgroup_labels <- sub_results %>%
    left_join(group_labels, by="subgroup") %>%
    mutate(label = subgroup)

  # ---- Plot ----
  p <- ggplot(dat, aes(x=exp(yi), y=y)) +

    # Contours
    geom_rect(data=contour_all, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=fill),
              alpha=0.9, inherit.aes=FALSE) +
    scale_fill_manual(values=c("Small"="gray95","Medium"="gray85","Large"="gray70",
                               "Very Large"="gray50"), name="") +

    # Vertical line at 1
    geom_segment(aes(x = 1, xend = 1, y = ymin_plot, yend = ymax_plot),
                 linetype = "solid", color = "black", size = 0.8, inherit.aes = FALSE) +

    # Study points & CIs
    geom_errorbarh(aes(xmin=exp(ci.lb), xmax=exp(ci.ub)), height=0, color=CI.col) +
    geom_point(aes(size=weight), shape=15, color=study.col) +
    scale_size_continuous(range=c(2, square.size), guide="none") +

    # Headers
    annotate("text", x=study_x, y=max(dat$y)+1.5, label="Study", hjust=0, size=text_size+.3, fontface="bold") +
    annotate("text", x=c(xpos$EventsT, xpos$EventsC, xpos$Effect, xpos$Weight),
             y=max(dat$y)+1.5,
             label=c("Events (T)","Events (C)", paste0(measure," (95% CI)"), "Weight (%)"),
             fontface="bold", hjust=0, size=text_size+.3) +

    # Study labels
    geom_text(aes(x=xpos$EventsT, label=paste0(events_t,"/",n_t)), hjust=0, size=text_size) +
    geom_text(aes(x=xpos$EventsC, label=paste0(events_c,"/",n_c)), hjust=0, size=text_size) +
    geom_text(aes(x=xpos$Effect, label=paste0(round(exp(yi),2)," [",round(exp(ci.lb),2),"-",round(exp(ci.ub),2),"]")),
              hjust=0, size=text_size) +
    geom_text(aes(x=xpos$Weight, label=paste0(round(weight,1),"%")), hjust=0, size=text_size) +
    annotate("text", x=study_x, y=dat$y, label=dat$Study, hjust=0, size=text_size) +

    # Subgroup diamonds
    lapply(diamonds, function(df) geom_polygon(data=df, aes(x=x, y=y), fill=diamond.col, alpha=0.8, inherit.aes=FALSE)) +

    # Subgroup labels & effect sizes
    geom_text(data=subgroup_labels,
              aes(x=study_x, y=y_label, label=paste0("Subgroup: ", label)),
              hjust=0, size=text_size+0.5, fontface="bold", inherit.aes=FALSE) +
    lapply(names(diamonds), function(g) {
      df <- diamonds[[g]]
      res <- sub_results %>% filter(subgroup==g)
      y_d <- df$y[1]
      geom_text(aes(x=val_x, y=y_d,
                    label=paste0(round(exp(res$yi),2), " [",
                                 round(exp(res$ci.lb),2), "-",
                                 round(exp(res$ci.ub),2), "]")),
                hjust=0, size=text_size+0.5, fontface="bold", inherit.aes=FALSE)
    }) +

    # Heterogeneity
    geom_text(data=hetero, aes(x=study_x, y=y, label=label),
              hjust=0, size=text_size-0.2, fontface="italic", inherit.aes=FALSE) +

    # Overall diamond + labels + model
    geom_polygon(data=overall_diamond, aes(x=x, y=y), fill=overall.col, alpha=0.9, inherit.aes=FALSE) +
    annotate("text", x=study_x, y=y_overall + 0.2, label="Overall", hjust=0, size=text_size+0.5, fontface="bold") +
    annotate("text", x=study_x, y=y_overall - 1.8, label=model_label,
             hjust=0, size=text_size, fontface="bold", color="black") +
    annotate("text", x=val_x, y=y_overall + 0.2,
             label=paste0(round(exp(overall$b),2), " [", round(exp(overall$ci.lb),2), "-", round(exp(overall$ci.ub),2), "]"),
             hjust=0, size=text_size+0.5, fontface="bold") +
    annotate("text", x=study_x, y=y_overall - 0.6, label=overall_label, hjust=0, size=text_size-0.2,
             fontface="italic") +
    # theme

    labs(x=xlab, y=NULL, title=title) +
    theme_minimal(base_size=14) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.y = element_blank(),
          legend.position = "top") +
    coord_cartesian(xlim = xlim, ylim=c(ymin_plot, max(dat$y)+3))

  # ---- Add prediction interval ----
  if(pred){
    p <- p +
      geom_errorbarh(aes(xmin=pi$pi.lb, xmax=pi$pi.ub, y=y_overall - 0.8),
                     height=0, size=Pred.Int.size, color=Pred.Inter.col, inherit.aes=FALSE) +
      annotate("text", x=xlim[1]+0.1, y=y_overall -1,
               label="Prediction interval", fontface="bold", hjust=0, size=text_size) +
      annotate("text", x=val_x, y=y_overall - 1,
               label=paste0("[", round(pi$pi.lb,2), "-", round(pi$pi.ub,2), "]"),
               fontface="bold", color="black", hjust=0, size=text_size)
  }

  return(p)
}
