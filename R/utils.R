# R/utils.R
# Shared helper functions and re-exports for the hypertension-gate-paper project.

# Ensure commonly-used packages are attached for legacy helper code that relies
# on them being on the search path.
suppressPackageStartupMessages({
  library(data.table)
})

# Load existing helper scripts so their functions are available when this file
# is sourced from the manuscript or analysis scripts.

if (file.exists("shared.helpers.R")) {
  source("shared.helpers.R")
}

if (file.exists("report.helpers.R")) {
  source("report.helpers.R")
}

