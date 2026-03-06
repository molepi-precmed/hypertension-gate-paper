# R/fetch_data.R
# Downloads the project data archive from a GitHub Release and extracts it into data/.
# Skipped automatically if the local sibling data directory is present.

LOCAL_DATA_DIR <- file.path("..", "hypertension-gate-data")
if (dir.exists(LOCAL_DATA_DIR)) {
  message("[fetch_data] Local data directory found. No download needed.")
  quit(save = "no", status = 0)
}

library(piggyback)

REPO         <- "molepi-precmed/hypertension-gate-paper"
DATA_TAG     <- "data-release"
DEST_DIR     <- "data"
ARCHIVE_NAME <- "hypertension-gate-data.tar.gz"

if (!dir.exists(DEST_DIR)) dir.create(DEST_DIR, recursive = TRUE)

# Use token if set (avoids 403: unauthenticated API limit is 60 req/hr)
token <- Sys.getenv("GITHUB_PAT", "")
if (token == "") {
  message("[fetch_data] GITHUB_PAT not set. If you see 403 errors, set it for higher rate limits.")
}

message("[fetch_data] Downloading ", ARCHIVE_NAME, " from GitHub Release '", DATA_TAG, "'...")
piggyback::pb_download(
  file      = ARCHIVE_NAME,
  dest      = DEST_DIR,
  repo      = REPO,
  tag       = DATA_TAG,
  overwrite = TRUE,
  .token    = if (nzchar(token)) token else NULL
)

archive_path <- file.path(DEST_DIR, ARCHIVE_NAME)
if (!file.exists(archive_path)) {
  stop("[fetch_data] Download failed: ", archive_path, " not found. Check that the release has asset '", ARCHIVE_NAME, "'.")
}

message("[fetch_data] Extracting archive into ", normalizePath(DEST_DIR), " ...")
untar(archive_path, exdir = DEST_DIR, tar = "internal")
unlink(archive_path)

message("[fetch_data] Done. Data are in ", normalizePath(DEST_DIR), " (eqtl/, pqtl/, mrresults/, etc.).")
