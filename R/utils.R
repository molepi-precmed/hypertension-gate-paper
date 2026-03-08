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

# Format the GWAS.hit / reported.genes column for tables:
#   - "NR" alone (or all-NR comma-list) -> "+"
#   - mixed gene symbols + NR -> remove NR entries, italicise remaining symbols
#   - gene symbols only -> italicise
format_gwas_hit <- function(x, latex = TRUE) {
  sapply(x, function(val) {
    if (is.na(val) || val == "") return(val)
    parts <- trimws(strsplit(as.character(val), ",")[[1]])
    non_nr <- parts[parts != "NR"]
    if (length(non_nr) == 0) {
      return("+")
    } else if (latex) {
      return(paste(paste0("\\textit{", non_nr, "}"), collapse = ", "))
    } else {
      return(paste(non_nr, collapse = ", "))
    }
  }, USE.NAMES = FALSE)
}

