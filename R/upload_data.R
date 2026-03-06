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

# Use relative path (from SOURCE_DIR) as release asset name so download preserves structure (data/eqtl/..., data/pqtl/..., etc.)
base_dir <- normalizePath(SOURCE_DIR, mustWork = TRUE, winslash = "/")
rel_name <- function(f) {
  full <- normalizePath(f, mustWork = TRUE, winslash = "/")
  sub(paste0("^", gsub("([.?*+^$[\\\\]()])", "\\\\\\1", base_dir), "/?"), "", full)
}
data_names <- vapply(data_files, rel_name, character(1L))

# pb_upload expects one file and one name per call (API expects scalar name)
message("[upload_data] Uploading ", length(data_files), " data file(s)...")

for (k in seq_along(data_files)) {
  piggyback::pb_upload(
    file      = data_files[k],
    repo      = REPO,
    tag       = DATA_TAG,
    name      = data_names[k],
    overwrite = "use_timestamps"
  )
  if (k %% 50L == 0L) gc()
  if (k %% 50L == 0L) message("  uploaded ", k, " / ", length(data_files))
}

# Upload manifest last so its timestamp confirms all data is present
if (file.exists(manifest_file)) {
  piggyback::pb_upload(
    file      = manifest_file,
    repo      = REPO,
    tag       = DATA_TAG,
    name      = "data_manifest.csv",
    overwrite = "use_timestamps"
  )
}

message("[upload_data] All files uploaded to release '", DATA_TAG, "'.")

