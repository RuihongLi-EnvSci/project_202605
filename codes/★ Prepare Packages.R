# This code requires R version 4.5.0 to run.Several packages and functions cannot run in the latest version.

options(timeout = 9999)

pre.packages <- c(
  "adespatial",
  "agricolae",
  "betapart",
  "bipartite",
  "car",
  "codyn",
  "dplyr",
  "ecotraj",
  "effectsize",
  "extraoperators",
  "GGally",
  "ggplot2",
  "ggrepel",
  "ggtern",
  "Hmisc",
  "igraph",
  "influence.ME",
  "JWileymisc",
  "labdsv",
  "lme4",
  "lmerTest",
  "multilevelTools",
  "MuMIn",
  "openxlsx",
  "piecewiseSEM",
  "purrr",
  "reshape2",
  "rstatix",
  "sf",
  "sjstats",
  "smacof",
  "stringr",
  "tibble",
  "tidyr",
  "tidyverse",
  "vegan"
)

installed <- installed.packages()[, "Package"]
to_install <- setdiff(pre.packages, installed)

if (length(to_install) > 0) {
  install.packages(to_install, dependencies = TRUE)
}

invisible(lapply(pre.packages, library, character.only = TRUE))

if (!requireNamespace("pairwiseAdonis", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
  }
  
  devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
}

library(pairwiseAdonis)
