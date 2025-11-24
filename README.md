# Set working directory
setwd("E:/AJO-EJO/metamed/metamed1")

# Define README.md path
file_path <- "README.md"

# Write content to README.md
writeLines(c(
"# contourForest",
"",
"This is my contour-enhanced forest plot package. It provides functions to create **contour-enhanced forest plots** for meta-analysis of **continuous** and **binary** outcomes. The package works with `metafor` and `ggplot2` and allows customization of colors, labels, and prediction intervals.",
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
"## Usage",
"",
"```r",
"# Load package",
"library(contourForest)",
"",
"# -------------------------",
"# Example 1: Binary outcomes",
"# -------------------------",
"",
"# Load example dataset (included in package)",
"data(dat)",
"",
"# View the first rows",
"head(dat)",
"",
"# Generate binary forest plot",
"forest.binary(dat, measure = 'OR')",
"",
"# ----------------------------",
"# Example 2: Continuous outcomes",
"# ----------------------------",
"",
"# Load example dataset (included in package)",
"data(dat1)",
"",
"# View the first rows",
"head(dat1)",
"",
"# Generate continuous forest plot",
"forest.continuous(dat1, measure = 'SMD')",
"```"
), con = file_path)

# Check that the file was created
file.exists(file_path)
