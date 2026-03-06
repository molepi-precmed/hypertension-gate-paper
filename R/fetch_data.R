# R/fetch_data.R
# Downloads project data from GitHub Releases using piggyback.
# Skipped automatically if the local sibling data directory is present.
# Do not source config.R here: when neither sibling dir nor data/ exist,
# config would stop() before we can run the download.

LOCAL_DATA_DIR <- file.path("..", "hypertension-gate-data")
if (dir.exists(LOCAL_DATA_DIR)) {
  message("[fetch_data] Local data directory found. No download needed.")
  quit(save = "no", status = 0)
}

library(piggyback)

REPO        <- "molepi-precmed/hypertension-gate-paper"
DATA_TAG    <- "data-release"
DEST_DIR    <- "data"

if (!dir.exists(DEST_DIR)) dir.create(DEST_DIR, recursive = TRUE)

# Use token if set (avoids 403: unauthenticated API limit is 60 req/hr)
token <- Sys.getenv("GITHUB_PAT", "")
if (token == "") {
  message("[fetch_data] GITHUB_PAT not set. If you see 403 errors, set it for higher rate limits.")
}

message("[fetch_data] Downloading data files from GitHub Release '", DATA_TAG, "'...")

piggyback::pb_download(
  file      = NULL,
  dest      = DEST_DIR,
  repo      = REPO,
  tag       = DATA_TAG,
  overwrite = FALSE,
  .token    = if (nzchar(token)) token else NULL
)

message("[fetch_data] Download complete. Files in: ", normalizePath(DEST_DIR))

