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
# Example usage with metafor model
# --------------------------------

