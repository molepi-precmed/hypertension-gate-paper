# R/config.R
# Resolves DATA_DIR based on execution context:
#   1. Local development: sibling hypertension-gate-data/ directory exists
#   2. Remote/CI build:   data/ subdirectory populated by fetch_data.R
#
# When running on a host with system R libraries, packages are loaded from
# .libPaths() (e.g. /usr/local/lib/R/site-library, /usr/lib/R/site-library,
# /usr/lib/R/library). Set USE_SYSTEM_R_LIBS=1 to use system libs only (no installs);
# otherwise missing packages are installed from CRAN.

LOCAL_DATA_DIR  <- file.path("..", "hypertension-gate-data")

# System library paths to prefer when present (e.g. build host).
SYSTEM_LIB_PATHS <- c(
  "/usr/local/lib/R/site-library",
  "/usr/lib/R/site-library",
  "/usr/lib/R/library"
)
use_system_libs <- nzchar(Sys.getenv("USE_SYSTEM_R_LIBS", "")) ||
  any(dir.exists(SYSTEM_LIB_PATHS))
system_only <- nzchar(Sys.getenv("USE_SYSTEM_R_LIBS", ""))  # no installs, fail if missing
if (use_system_libs) {
  existing_system <- SYSTEM_LIB_PATHS[dir.exists(SYSTEM_LIB_PATHS)]
  if (length(existing_system) > 0L) {
    current <- .libPaths()
    .libPaths(c(existing_system, current[!current %in% existing_system]))
  }
}

#' Load a package. When USE_SYSTEM_R_LIBS=1: load only from .libPaths(), do not install.
#' Otherwise: install from CRAN if not found.
ensure_library <- function(pkg, repos = getOption("repos")[1L]) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (system_only) {
      stop(
        "Required package '", pkg, "' is not available in .libPaths(). ",
        "Install it in the system library or unset USE_SYSTEM_R_LIBS to allow installs.",
        call. = FALSE
      )
    }
    if (repos == "@CRAN@" || !length(repos) || !nzchar(repos))
      repos <- "https://cloud.r-project.org"
    install.packages(pkg, repos = repos, quiet = TRUE)
  }
  library(pkg, character.only = TRUE)
}

# Packages required by the pipeline and rendered manuscript (for ensure_library).
REQUIRED_PACKAGES <- c(
  "flextable", "data.table", "ggplot2", "writexl",
  "bigstatsr", "ggcorrplot", "foreach", "doParallel",
  "igraph", "ggraph", "graphlayouts", "kableExtra", "rmarkdown"
)
FETCHED_DATA_DIR <- file.path("data")

if (dir.exists(LOCAL_DATA_DIR)) {
  DATA_DIR <- LOCAL_DATA_DIR
  message("[config] Using LOCAL data directory: ", normalizePath(DATA_DIR))
} else if (length(list.files(FETCHED_DATA_DIR, all.files = FALSE)) > 0) {
  DATA_DIR <- FETCHED_DATA_DIR
  message("[config] Using FETCHED data directory: ", normalizePath(DATA_DIR))
} else {
  stop(
    "[config] No data found.\n",
    "  - For local dev: ensure '../hypertension-gate-data/' exists.\n",
    "  - For remote build: run `make data` to fetch data via piggyback."
  )
}

# Helper: build an absolute path to a data file
data_path <- function(...) file.path(DATA_DIR, ...)

