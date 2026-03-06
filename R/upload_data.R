# R/upload_data.R
# Uploads the local data directory as a single compressed archive to a GitHub Release.
# Run manually: Rscript R/upload_data.R
# Requires GITHUB_PAT environment variable with repo write permissions.

library(piggyback)

REPO        <- "molepi-precmed/hypertension-gate-paper"
DATA_TAG    <- "data-release"
SOURCE_DIR  <- file.path("..", "hypertension-gate-data")
ARCHIVE_NAME <- "hypertension-gate-data.tar.gz"

if (!dir.exists(SOURCE_DIR)) {
  stop("Local data directory not found: ", SOURCE_DIR)
}

# Create the release tag if it doesn't exist yet
existing_tags <- tryCatch(
  piggyback::pb_releases(repo = REPO)$tag_name,
  error = function(e) character(0)
)

if (!DATA_TAG %in% existing_tags) {
  message("[upload_data] Creating release tag '", DATA_TAG, "'...")
  piggyback::pb_new_release(repo = REPO, tag = DATA_TAG)
}

# Create archive: contents of SOURCE_DIR (eqtl/, pqtl/, *.csv, etc.) at root of tarball
message("[upload_data] Creating archive of ", SOURCE_DIR, " ...")
tarball_path <- tempfile(fileext = ".tar.gz")
owd <- setwd(SOURCE_DIR)
on.exit(setwd(owd))
files <- list.files(recursive = TRUE, include.dirs = FALSE)
if (length(files) == 0L) stop("No files found in ", SOURCE_DIR)
tar(tarball_path, files = files, compression = "gzip", tar = "internal")
setwd(owd)
on.exit(NULL)

message("[upload_data] Uploading ", ARCHIVE_NAME, " to release '", DATA_TAG, "'...")
piggyback::pb_upload(
  file      = tarball_path,
  repo      = REPO,
  tag       = DATA_TAG,
  name      = ARCHIVE_NAME,
  overwrite = TRUE
)
unlink(tarball_path)

message("[upload_data] Done. Single asset '", ARCHIVE_NAME, "' uploaded to release '", DATA_TAG, "'.")
