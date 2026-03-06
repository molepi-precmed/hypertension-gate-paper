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

# Gather all files to upload (recursive: include files inside eqtl/, pqtl/, mrresults/, etc.)
all_files     <- list.files(SOURCE_DIR, full.names = TRUE, recursive = TRUE)
all_files     <- all_files[!file.info(all_files)$isdir]  # only regular files, not directories
manifest_file <- file.path(SOURCE_DIR, "data_manifest.csv")
data_files    <- setdiff(all_files, manifest_file)

BATCH_SIZE <- 50L  # upload in batches to avoid exhausting R's connection limit (~128)

message("[upload_data] Uploading ", length(data_files), " data file(s) in batches of ", BATCH_SIZE, "...")

for (i in seq(1L, length(data_files), by = BATCH_SIZE)) {
  batch <- data_files[seq(i, min(i + BATCH_SIZE - 1L, length(data_files)))]
  piggyback::pb_upload(
    file      = batch,
    repo      = REPO,
    tag       = DATA_TAG,
    overwrite = "use_timestamps"
  )
  gc()  # release connections between batches
  message("  uploaded ", min(i + BATCH_SIZE - 1L, length(data_files)), " / ", length(data_files))
}

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

