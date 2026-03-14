# R/fetch_data.R
# Downloads the project data archive from a GitHub Release and extracts it into data/.
# Skipped if the local sibling data directory exists, or if data/ already contains files.
# Pass --update to force a re-download when the remote asset is newer than local data.

args   <- commandArgs(trailingOnly = TRUE)
UPDATE <- "--update" %in% args

LOCAL_DATA_DIR <- file.path("..", "hypertension-gate-data")
if (dir.exists(LOCAL_DATA_DIR)) {
  message("[fetch_data] Local data directory found. No download needed.")
  quit(save = "no", status = 0)
}

DEST_DIR <- "data"
existing <- setdiff(list.files(DEST_DIR, all.files = FALSE), ".gitkeep")

if (!UPDATE && length(existing) > 0L) {
  message("[fetch_data] data/ already contains files. No download needed.")
  message("[fetch_data] Run with --update to re-download if the remote data have changed.")
  quit(save = "no", status = 0)
}

library(piggyback)

REPO         <- "molepi-precmed/hypertension-gate-paper"
DATA_TAG     <- "data-release"
ARCHIVE_NAME <- "hypertension-gate-data.tar.gz"

# Use token if set (avoids 403: unauthenticated API limit is 60 req/hr)
token <- Sys.getenv("GITHUB_PAT", "")
if (token == "") {
  message("[fetch_data] GITHUB_PAT not set. If you see 403 errors, set it for higher rate limits.")
}

SENTINEL <- file.path(DEST_DIR, ".last-updated")

# If --update and local files exist, compare remote timestamp to sentinel mtime.
if (UPDATE && length(existing) > 0L) {
  message("[fetch_data] --update: checking remote asset timestamp...")
  info <- tryCatch(
    piggyback::pb_list(repo = REPO, tag = DATA_TAG,
                       .token = if (nzchar(token)) token else NULL),
    error = function(e) {
      message("[fetch_data] Could not fetch release info: ", conditionMessage(e))
      NULL
    }
  )
  if (!is.null(info)) {
    asset_row <- info[info$file_name == ARCHIVE_NAME, ]
    if (nrow(asset_row) > 0L) {
      remote_ts <- asset_row$timestamp[1L]
      if (!inherits(remote_ts, "POSIXct"))
        remote_ts <- as.POSIXct(remote_ts, tz = "UTC")
      local_ts <- if (file.exists(SENTINEL)) file.info(SENTINEL)$mtime else as.POSIXct(NA)
      if (!is.na(remote_ts) && !is.na(local_ts) && remote_ts <= local_ts) {
        message("[fetch_data] Remote asset (", format(remote_ts, tz = "UTC"), " UTC) is not newer than local data (",
                format(local_ts, tz = "UTC"), " UTC). Nothing to do.")
        quit(save = "no", status = 0)
      }
      message("[fetch_data] Remote asset updated (", format(remote_ts, tz = "UTC"),
              " UTC). Re-downloading...")
    }
  }
}

if (!dir.exists(DEST_DIR)) dir.create(DEST_DIR, recursive = TRUE)

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

writeLines(format(Sys.time(), tz = "UTC"), SENTINEL)
message("[fetch_data] Done. Data are in ", normalizePath(DEST_DIR), " (eqtl/, pqtl/, mrresults/, etc.).")
