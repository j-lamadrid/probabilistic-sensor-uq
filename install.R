# 1. Install BiocManager to handle Bioconductor packages
install.packages("BiocManager", repos = "https://cloud.r-project.org")

# 2. Install all standard CRAN packages
install.packages(c(
  "shiny",
  "plotly",
  "DT",
  "dplyr",
  "tidyr",
  "e1071",
  "progressr",
  "pracma",
  "viridis",
  "rstan",
  "ggplot2",
  "here",
  "ggrepel",
  "GGally",
  "DiagrammeR",
  "ggridges"
), repos = "https://cloud.r-project.org")

# 3. Install Mass Spectrometry Bioconductor packages
BiocManager::install(c(
  "MSnbase",
  "mzR",
  "BiocParallel"
), update = FALSE, ask = FALSE)