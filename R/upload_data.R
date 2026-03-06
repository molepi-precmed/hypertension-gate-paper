# R/upload_data.R
# Uploads all files from the local data directory to a GitHub Release.
# Run manually: Rscript R/upload_data.R
# Requires GITHUB_PAT environment variable with repo write permissions.

library(piggyback)

REPO        <- "molepi-precmed/hypertension-gate-paper"
DATA_TAG    <- "data-release"
SOURCE_DIR  <- file.path("..", "hypertension-gate-data")

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

# Gather all files to upload (exclude the manifest, upload it last)
all_files     <- list.files(SOURCE_DIR, full.names = TRUE, recursive = FALSE)
manifest_file <- file.path(SOURCE_DIR, "data_manifest.csv")
data_files    <- setdiff(all_files, manifest_file)

message("[upload_data] Uploading ", length(data_files), " data file(s)...")

piggyback::pb_upload(
  file      = data_files,
  repo      = REPO,
  tag       = DATA_TAG,
  overwrite = "use_timestamps"
)

# Upload manifest last so its timestamp confirms all data is present
if (file.exists(manifest_file)) {
  piggyback::pb_upload(
    file      = manifest_file,
    repo      = REPO,
    tag       = DATA_TAG,
    overwrite = "use_timestamps"
  )
}

message("[upload_data] All files uploaded to release '", DATA_TAG, "'.")

