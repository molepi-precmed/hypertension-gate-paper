# R/utils.R
# Shared helper functions and re-exports for the hypertension-gate-paper project.

# Ensure commonly-used packages are attached for legacy helper code that relies
# on them being on the search path.
suppressPackageStartupMessages({
  library(data.table)
})

# Load existing helper scripts so their functions are available when this file
# is sourced from the manuscript or analysis scripts.

if (file.exists("shared.helpers.R")) {
  source("shared.helpers.R")
}

if (file.exists("report.helpers.R")) {
  source("report.helpers.R")
}

# Format the GWAS.hit / reported.genes column for tables:
#   - "NR" alone (or all-NR comma-list) -> "+"
#   - mixed gene symbols + NR -> remove NR entries, italicise remaining symbols
#   - gene symbols only -> italicise
#' Plot a correlation heatmap with rows/columns ordered by hierarchical
#' clustering on squared correlations.
#'
#' @param cor_mat Numeric matrix of correlations (symmetric; rownames = colnames = labels)
#' @param colors Length-3 character vector: low, mid, high fill colours
#' @return A ggplot object
correlation_heatmap <- function(cor_mat,
                                colors = c("#2166AC", "white", "#B2182B")) {
  # Align columns to rownames (stored order may differ)
  cor_mat <- cor_mat[, rownames(cor_mat)]

  d   <- as.dist(1 - cor_mat^2)
  hc  <- hclust(d, method = "complete")
  ord <- hc$order
  genes <- rownames(cor_mat)[ord]
  mat   <- cor_mat[genes, genes]

  df <- data.frame(
    x = factor(rep(genes, each = length(genes)), levels = genes),
    y = factor(rep(genes, length(genes)),        levels = rev(genes)),
    r = as.vector(mat)
  )

  ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, fill = r)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
    ggplot2::geom_text(
      data = df[abs(df$r) >= 0.5, ],
      ggplot2::aes(label = sprintf("%.1f", r)),
      colour = "white", size = 2.5, fontface = "bold"
    ) +
    ggplot2::scale_fill_gradient2(
      low = colors[1], mid = colors[2], high = colors[3],
      midpoint = 0, limits = c(-1, 1), name = "Correlation"
    ) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_text(angle = 45, hjust = 1, face = "italic", size = 11),
      axis.text.y  = ggplot2::element_text(face = "italic", size = 11),
      axis.title   = ggplot2::element_blank(),
      panel.grid   = ggplot2::element_blank(),
      legend.position = "right"
    ) +
    ggplot2::coord_fixed()
}

format_gwas_hit <- function(x, latex = TRUE) {
  sapply(x, function(val) {
    if (is.na(val) || val == "") return(val)
    parts <- trimws(strsplit(as.character(val), ",")[[1]])
    non_nr <- parts[parts != "NR"]
    if (length(non_nr) == 0) {
      return("+")
    } else if (latex) {
      return(paste(paste0("\\textit{", non_nr, "}"), collapse = ", "))
    } else {
      return(paste(non_nr, collapse = ", "))
    }
  }, USE.NAMES = FALSE)
}

