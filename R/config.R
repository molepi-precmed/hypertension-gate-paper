# R/config.R
# Resolves DATA_DIR based on execution context:
#   1. Local development: sibling hypertension-gate-data/ directory exists
#   2. Remote/CI build:   data/ subdirectory populated by fetch_data.R

LOCAL_DATA_DIR  <- file.path("..", "hypertension-gate-data")
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

