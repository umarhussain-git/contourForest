# contourforest

**Contour-Enhanced Forest Plots for Meta-Analysis**

This R package provides functions to create **contour-enhanced forest plots** for meta-analysis of **continuous** and **binary** outcomes. It integrates with `metafor` and `ggplot2` and allows customization of colors, labels, prediction intervals, and subgroup analyses.

---

```
# contourforest

A detailed guide for using the `contourforest` package.  

## 📖 Read the Full Guide
For step-by-step instructions and examples, visit the RPubs guide:  
[https://rpubs.com/umarhussain/1382129](https://rpubs.com/umarhussain/1382129)

## 🚀 How to Run on Different Types of Datasets
The guide covers running `contourforest` on various dataset formats and scenarios. Follow the instructions in the RPubs tutorial to get started quickly.


```

## Installation

```r
# Install remotes if not already installed
# install.packages('remotes')

remotes::install_github('umarhussain-git/contourforest')
```

## Functions 


bcg	- BCG Vaccine Trials Dataset
dat1- 	Hypothetical dataset: GAD-7 Anxiety Scores by Socioeconomic Status
forest_bin()	- Contour-Enhanced Binary Outcome Forest Plot
forest_bin_subgroup()	- Subgroup Forest Plot for Binary Outcome Meta-analysis
forest_cont() - 	Contour-enhanced Forest Plot for Continuous Outcomes
forest_cont_subgroup() - 	Forest Plot for Subgroup Meta-Analysis
## Example Datasets
dat – Hypothetical dataset: Depression by Socioeconomic Status

dat1 – Hypothetical dataset: GAD-7 Anxiety Scores by Socioeconomic Status
library(contourForest)

# Example 1: Binary Outcomes
#` Load example dataset


data(bcg())

#' Generate a binary forest plot (Odds Ratio)

```r
forest_bin(
  dat = bcg(),
  measure = "OR",
  xlab = "Odds Ratio",
  title = "BCG Vaccine Meta-analysis",
  tlim = c(0, 2.3),
  contour_left_min = c(0,0.3,0.5,0.7),
  contour_left_max = c(0.3,0.5,0.7,1),
  contour_right_min = c(1,1.2,1.5,1.8),
  contour_right_max = c(1.2,1.5,1.8,2.5)
)
```

library(contourforest)

#' Load example dataset

data(dat1)

#' View first rows

head(dat1)
### Example 2: Continuous Outcomes
#' Generate a continuous forest plot (Standardized Mean Difference or Mean Difference)
```r
forest_cont(
  dat1,
  measure = "MD",
  xlab = "Mean Difference",
  study_x = -9,
  sort = "effect",
  hetero_x = -12,
  treatment_x = -7,
  control_x = -5,
  effect_x = 5.5,
  weight_x = 10,
  PredInt_x = 7
)

forest_cont(
  dat1,
  measure = "SMD",
  xlab = "Standardized Mean Difference",
  hetero_x = -9.9,
  study_x = -7,
  sort = "effect",
  treatment_x = -5,
  control_x = -3,
  effect_x = 2.5,
  weight_x = 4,
  PredInt_x = 4
)
forest_cont(
dat1,
study.col = "darkgreen",
CI.col = "black",
diamond.col = "red",
Pred.Inter.col = "black",
measure = "SMD",
sort = "effect",
xlab = "Standardized Mean Difference",
contour_fill = c("gray90","gray70","gray50", "gray30"),
hetero_x = -9.9,
study_x = -7,
square.size = 9,
treatment_x = -5,
control_x = -3.2,
text_size = 4,
effect_x = 2.5,
weight_x = 5.8,
PredInt_x = 3
)
```

### Example 3: Binary Outcomes with Subgroups
library(contourForest)

#, Generate a binary forest plot with subgroups
```r
forest_bin_subgroup(bcg(), measure = 'OR')
```

### Example 4: Continuous Outcomes with Subgroups
library(contourforest)

#` Generate a continuous forest plot with subgroups

```r
# Using the built-in dataset dat1
forest_cont_subgroup(dat1)
```



## Features

Contour-enhanced visualization to highlight statistical significance

Supports Odds Ratios (OR), Risk Ratios (RR), Mean Differences (MD), and Standardized Mean Differences (SMD)

Optional prediction intervals

Customizable colors, labels, and annotation positions

Subgroup analysis with pooled effects

forest_cont() – Forest plot for continuous outcomes

forest_cont_subgroup() – Forest plot for continuous outcomes with subgroups
# Maintainer 
## Umar Hussain

drumarhussain@gmail.com

Assistant professor, Orthodontics

Saidu college of Dentistry


