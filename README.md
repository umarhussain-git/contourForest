# contourForest
contourForest is an R package for creating contour-enhanced forest plots for meta-analysis of both continuous and binary outcomes.
It integrates seamlessly with metafor and ggplot2, offering flexible customization for:
Significance contours
Effect size labels
Colors and themes
Prediction and confidence intervals

This package is designed to make forest plots more informative, interpretable, and visually appealing for publications and presentations.

# Installation
## Install 'remotes' if not already installed
## install.packages("remotes")

remotes::install_github("umarhussain-git/contourForest")

# Usage Example
## Load the package
library(contourForest)

# --------------------------------
library(contourForest)

# Load included example dataset
data(dat)

# Preview data
head(dat)

# Generate a binary contour-enhanced forest plot
forest.binary(dat, measure = "OR")

# Generate a binary contour-enhanced forest plot with subgroup
forest_subgroup_binary(dat)

# Load included example dataset
data(dat1)

# Preview data
head(dat1)

# Generate a continuous contour-enhanced forest plot
forest.continuous(dat1, measure = "SMD")

# --------------------------------

