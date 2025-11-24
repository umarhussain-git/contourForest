# contourForest
This is my contour-enhanced forest plot package. It provides functions to create **contour-enhanced forest plots** for meta-analysis of **continuous** and **binary** outcomes. The package works with `metafor` and `ggplot2` and allows customization of colors, labels, and prediction intervals.",
"",
"---",
"",
"## Installation",
"",
"```r",
"# install remotes if not already installed",
"# install.packages('remotes')",
"",
"remotes::install_github('umarhussain-git/contourForest')",
"```",
"",
"---",
"",
# Load the package
library(contourForest)

# -------------------------
# Example 1: Binary outcomes with subgroups
# -------------------------

# Load example dataset included in package
data(dat)

# View the first rows
head(dat)

# Generate binary forest plot with subgroups
forest.binary.subgroup(dat, measure = "OR")

# -------------------------
# Example 2: Continuous outcomes
# -------------------------

# Load example dataset included in package
data(dat1)

# View the first rows
head(dat1)

# Generate continuous forest plot
forest.continuous(dat1, measure = "SMD")


