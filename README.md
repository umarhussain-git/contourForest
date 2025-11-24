# contourForest

This is my contour-enhanced forest plot package. It provides functions to create **contour-enhanced forest plots** for meta-analysis of **continuous** and **binary** outcomes. The package works with `metafor` and `ggplot2` and allows customization of colors, labels, and prediction intervals.

---

## Installation

Install the development version directly from GitHub:

```r
# install devtools if not already installed
# install.packages("devtools")

devtools::install_github("umarhussain-git/contourForest")




```r
# Load package
library(contourForest)

# -------------------------
# Example 1: Binary outcomes
# -------------------------

# Load example dataset (included in package)
data(dat)

# View the first rows
head(dat)

# Generate binary forest plot
forest.binary(dat, measure = "OR")

# ----------------------------
# Example 2: Continuous outcomes
# ----------------------------

# Load example dataset (included in package)
data(dat1)

# View the first rows
head(dat1)

# Generate continuous forest plot
forest.continuous(dat1, measure = "SMD")



✅ This version is clean, fully reproducible, and ready for your GitHub README.md.  

If you want, I can **also add a Subgroup example** using `dat$subgroup` so users can see how to make stratified forest plots directly in the README. It only needs a few extra lines. Do you want me to do that?
