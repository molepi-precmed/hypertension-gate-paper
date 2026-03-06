# R/fetch_data.R
# Downloads project data from GitHub Releases using piggyback.
# Skipped automatically if the local sibling data directory is present.

source("R/config.R")  # will stop early if local data exists (already resolved)

# If we reach here and DATA_DIR is the local directory, nothing to do.
if (DATA_DIR == file.path("..", "hypertension-gate-data")) {
  message("[fetch_data] Local data directory found. No download needed.")
  quit(save = "no", status = 0)
}

library(piggyback)

REPO        <- "molepi-precmed/hypertension-gate-paper"
DATA_TAG    <- "data-release"
DEST_DIR    <- "data"

if (!dir.exists(DEST_DIR)) dir.create(DEST_DIR, recursive = TRUE)

message("[fetch_data] Downloading data files from GitHub Release '", DATA_TAG, "'...")

piggyback::pb_download(
  file      = NULL,          # NULL = download all files attached to the release
  dest      = DEST_DIR,
  repo      = REPO,
  tag       = DATA_TAG,
  overwrite = FALSE,         # use_timestamps: only download if remote is newer
  .token    = gh::gh_token() # reads GITHUB_PAT env var; empty string for public repos
)

message("[fetch_data] Download complete. Files in: ", normalizePath(DEST_DIR))

