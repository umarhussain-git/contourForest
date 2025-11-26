# contourForest

**Contour-Enhanced Forest Plots for Meta-Analysis**

This R package provides functions to create **contour-enhanced forest plots** for meta-analysis of **continuous** and **binary** outcomes. It integrates with `metafor` and `ggplot2` and allows customization of colors, labels, prediction intervals, and subgroup analyses.

---

## Installation

```r
# Install remotes if not already installed
# install.packages('remotes')

remotes::install_github('umarhussain-git/contourForest')
```

## Functions 

forest.binary() – Contour-enhanced forest plot for binary outcomes (OR or RR)

forest.binary.subgroup() – Contour-enhanced forest plot for binary outcomes by subgroups

## Example Datasets
dat – Hypothetical dataset: Depression by Socioeconomic Status

dat1 – Hypothetical dataset: GAD-7 Anxiety Scores by Socioeconomic Status
library(contourForest)

# Example 1: Binary Outcomes
#` Load example dataset
data(dat)

#' View first rows
head(dat)

#' Generate a binary forest plot (Odds Ratio)
forest.binary(dat, measure = 'OR')


library(contourForest)

#' Load example dataset
data(dat1)

#' View first rows
head(dat1)
### Example 2: Continuous Outcomes
#' Generate a continuous forest plot (Standardized Mean Difference)
forest.continuous(dat1, measure = 'SMD')

### Example 3: Binary Outcomes with Subgroups
library(contourForest)

#, Generate a binary forest plot with subgroups
forest.binary.subgroup(dat, measure = 'OR', subgroup_col = 'group')

### Example 4: Continuous Outcomes with Subgroups
library(contourForest)

#` Generate a continuous forest plot with subgroups
forest.continuous.subgroup(dat1, m_t_col = "mean_t", sd_t_col = "sd_t",

                           n_t_col = "n_t", m_c_col = "mean_c", sd_c_col = "sd_c",
                           
                           n_c_col = "n_c", subgroup_col = "group",
                           
                           study_col = "Study", measure = "SMD")


## Features

Contour-enhanced visualization to highlight statistical significance

Supports Odds Ratios (OR), Risk Ratios (RR), Mean Differences (MD), and Standardized Mean Differences (SMD)

Optional prediction intervals

Customizable colors, labels, and annotation positions

Subgroup analysis with pooled effects
forest.continuous() – Forest plot for continuous outcomes
# Maintainer 
## Umar Hussain

drumarhussain@gmail.com
Saidu college of Dentistry

forest.continuous.subgroup() – Forest plot for continuous outcomes with subgroups
