## contourForest
*contourForest* is an R package for creating contour-enhanced forest plots for meta-analysis of both continuous and binary outcomes.
It integrates seamlessly with metafor and ggplot2, offering flexible customization for:
Significance contours
Effect size labels
Colors and themes
Prediction and confidence intervals

This package is designed to make forest plots more informative, interpretable, and visually appealing for publications and presentations.

## Installation

 Install 'remotes' if not already installed

#` install.packages("remotes")

remotes::install_github("umarhussain-git/contourForest")

### Usage Example
#` Load the package
library(contourForest)

#` --------------------------------


library(contourForest)

# Load included example dataset
data(dat)

### Preview data
head(dat)

### Generate a binary contour-enhanced forest plot
forest.binary(dat, measure = "OR")

### Generate a binary contour-enhanced forest plot with subgroup
forest.binary.subgroup(
  dat,
  subgroup = NULL,
  measure = "RR",
  method = "REML",
  xlab = "Risk Ratio (RR)",
  title = "Subgroup Forest Plot",
  diamond.col = "red",
  overall.col = "darkgreen",
  study.col = "blue",
  CI.col = "blue",
  Pred.Inter.col = "black",
  square.size = 8,
  Pred.Int.size = 2,
  xlim = c(2,4),
  text_size = 3.5,
  xpos = list(EventsT = -1.4, EventsC = -0.9, Effect = 2.5, Weight = 3.2),
  study_x = -1.8,
  val_x = 2.6,
  pred = TRUE
)

### Load included example dataset
data(dat1)

### Preview data
> head(dat1)
    Study mean_c sd_c n_c mean_t sd_t n_t
1 Study A    5.2  1.8  50    7.4  2.0  50
2 Study B    4.8  2.0  60    5.9  1.9  60
3 Study C    6.1  2.1  45    8.0  2.2  45
4 Study D    5.5  1.9  55    3.8  2.1  55
5 Study E    4.9  2.2  50    4.2  2.0  50
6 Study F    5.8  1.7  40    8.1  2.3  40

### Generate a continuous contour-enhanced forest plot
forest.continuous(dat1, measure = "SMD")


# --------------------------------

