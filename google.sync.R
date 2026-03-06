library(trackdown)

# Helper function for consistent parameters
track_paper <- function(action = "download") {
  switch(action,
    upload = trackdown::upload_file(
      file = "hypertension_paper.Rmd",
      gfile = "hypertension_paper",
      gpath = "hypertension_paper",
      hide_code = TRUE
    ),
    update = trackdown::update_file(
      file = "hypertension_paper.Rmd",
      gfile = "hypertension_paper",
      gpath = "hypertension_paper",
      hide_code = TRUE
    ),
    download = trackdown::download_file(
      file = "hypertension_paper.Rmd",
      gfile = "hypertension_paper",
      gpath = "hypertension_paper"
    ),
    render = trackdown::render_file(
      file = "hypertension_paper.Rmd",
      gfile = "hypertension_paper",
      gpath = "hypertension_paper"
    )
  )
}

# Usage: track_paper("download"), track_paper("update"), etc.
