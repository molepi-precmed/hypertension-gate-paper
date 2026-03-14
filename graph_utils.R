## graph_utils.R
## Functions to plot bipartite graph of GWAS-hit trans-eQTL clumps and target genes

library(data.table)
library(igraph)
library(ggraph)
library(graphlayouts)
library(ggplot2)

#' Plot bipartite graph of GWAS-hit trans-eQTL clumps and their target genes
#'
#' Identifies clumps containing GWAS hits, filters target genes appearing in
#' at least min_gwas_clumps such clumps, and draws a bipartite graph.
#'
#' @param clumps_dt data.table with clump, chrom, start, stop, targetgenes, nearby_genes
#' @param gwas_dt data.table from fetch_gwas_associations (with chrom, position, reported_genes)
#' @param coregenes_dt data.table with gene_symbol, pvalue_trans (for node sizing)
#' @param min_gwas_clumps Integer. Minimum number of GWAS-hit clumps a target gene
#'   must be connected to (default 2)
#' @return A ggplot object
plot_clump_target_graph <- function(clumps_dt, gwas_dt, coregenes_dt,
                                    min_gwas_clumps = 2,
                                    label_size = 2.8,
                                    count_size = 2.2,
                                    legend_text_size = NULL,
                                    clump_width = 1.0,
                                    target_width = NULL,
                                    count_nudge_y = -0.8,
                                    row_gap = 1.0,
                                    clump_node_size = 3.5,
                                    clump_label_nudge_y = 0.1,
                                    target_label_nudge_y = -0.15) {

    ## --- 1. Identify GWAS-hit clumps ---
    ## For each clump, check if any GWAS SNP falls within its boundaries
    clumps_dt <- copy(clumps_dt)
    clumps_dt[, gwas_genes := vapply(seq_len(.N), function(i) {
        chr <- as.character(chrom[i])
        start_bp <- start[i] * 1e6
        stop_bp  <- stop[i] * 1e6
        gc <- sub("^chr", "", chr)
        hits <- gwas_dt[chrom == gc &
                        !is.na(position) &
                        position >= start_bp &
                        position <= stop_bp]
        if (nrow(hits) == 0) return(NA_character_)
        rg <- hits$reported_genes[!is.na(hits$reported_genes)]
        if (length(rg) == 0) return("+")
        genes <- unique(trimws(unlist(strsplit(rg, ",\\s*"))))
        genes <- genes[!genes %in% c("NR", "Intergenic", "intergenic", "")]
        genes <- genes[!grepl("^LOC|^LINC|^RN[U7]|^MIR\\d|^SNORD", genes)]
        if (length(genes) == 0) return("+")
        paste(genes, collapse = "/")
    }, character(1))]

    ## Exclude HLA region (25-34 Mb on chromosome 6)
    clumps_dt <- clumps_dt[!(as.character(chrom) == "6" & stop >= 25 & start <= 34)]

    gwas_clumps <- clumps_dt[!is.na(gwas_genes)]
    if (nrow(gwas_clumps) == 0) {
        warning("No GWAS-hit clumps found")
        return(ggplot() + theme_void())
    }

    ## --- 2. Build edge list: clump -> target gene ---
    edges <- gwas_clumps[, {
        targets <- trimws(unlist(strsplit(targetgenes, ",\\s*")))
        targets <- targets[targets != "" & targets != "."]
        .(target = targets)
    }, by = .(clump, gwas_genes)]

    ## Filter to targets in coregenes_dt (Table 1)
    edges <- edges[target %in% coregenes_dt$gene_symbol]

    ## --- 3. Filter targets by min_gwas_clumps ---
    target_counts <- edges[, .(n_clumps = uniqueN(clump)), by = target]
    keep_targets <- target_counts[n_clumps >= min_gwas_clumps, target]
    edges <- edges[target %in% keep_targets]

    if (nrow(edges) == 0) {
        warning("No target genes connected to >= ", min_gwas_clumps, " GWAS-hit clumps")
        return(ggplot() + theme_void())
    }

    ## --- 4. Create clump labels (chromosome + GWAS gene names) ---
    clump_info <- gwas_clumps[clump %in% edges$clump,
                               .(label = ifelse(gwas_genes == "+",
                                                paste0("<", chrom, ":", sprintf("%.1f", start), ">"),
                                                paste0(chrom, ": ", gwas_genes))),
                               by = clump]
    ## Deduplicate: one label per clump (take first if multiple)
    clump_info <- clump_info[, .(label = label[1]), by = clump]
    edges <- merge(edges, clump_info, by = "clump")

    ## --- 5. Build igraph ---
    ## Node names: clump labels (prefixed to avoid collision) and target genes
    edges[, clump_node := paste0("clump_", clump)]
    el <- edges[, .(clump_node, target)]
    g <- graph_from_data_frame(el, directed = FALSE)

    ## Set node attributes
    V(g)$type <- V(g)$name %in% edges$target  # TRUE = target gene, FALSE = clump
    V(g)$is_clump <- !V(g)$type

    ## Clump labels
    clump_labels <- clump_info$label
    names(clump_labels) <- paste0("clump_", clump_info$clump)
    V(g)$display_label <- ifelse(
        V(g)$is_clump,
        clump_labels[V(g)$name],
        V(g)$name
    )

    ## Target gene node size: proportional to -log10(pvalue_trans)
    ## Direction of effect: sign of standardized log OR
    pval_map <- coregenes_dt[, .(gene_symbol,
                                  neglog10p = -log10(pvalue_trans),
                                  logor = estimate_trans * sqrt(gw.variance))]
    ## For genes with multiple entries, take the best (max -log10p)
    pval_map <- pval_map[, .SD[which.max(neglog10p)], by = gene_symbol]
    V(g)$neglog10p <- ifelse(
        V(g)$is_clump, NA_real_,
        pval_map$neglog10p[match(V(g)$name, pval_map$gene_symbol)]
    )
    ## Direction: "positive" (blue) or "negative" (red)
    V(g)$direction <- ifelse(
        V(g)$is_clump, NA_character_,
        ifelse(pval_map$logor[match(V(g)$name, pval_map$gene_symbol)] >= 0,
               "Positive", "Negative")
    )

    ## Edge count per target gene
    target_counts_final <- edges[target %in% V(g)$name[!V(g)$is_clump],
                                  .(n_edges = uniqueN(clump)), by = target]
    V(g)$n_edges <- ifelse(
        V(g)$is_clump, NA_integer_,
        target_counts_final$n_edges[match(V(g)$name, target_counts_final$target)]
    )

    ## --- 6. Two-layer bipartite layout ---
    ## Use stress layout to get sensible x-coordinates, then fix y by layer
    set.seed(42)
    layout_xy <- layout_with_stress(g)

    ## Clumps at top (y = 1), targets at bottom (y = 0)
    is_clump <- V(g)$is_clump
    layout_xy[is_clump, 2]  <- row_gap
    layout_xy[!is_clump, 2] <- 0

    ## Spread x-coordinates evenly within each layer to avoid overlap
    clump_idx  <- which(is_clump)
    target_idx <- which(!is_clump)

    ## Order clumps by their stress-layout x to preserve neighbourhood structure
    clump_order  <- clump_idx[order(layout_xy[clump_idx, 1])]
    target_order <- target_idx[order(layout_xy[target_idx, 1])]

    ## Evenly space each layer; widen target spacing for label readability
    n_clumps  <- length(clump_order)
    n_targets <- length(target_order)
    layout_xy[clump_order, 1]  <- seq(0, clump_width, length.out = n_clumps)
    if (is.null(target_width))
        target_width <- max(1, n_targets - 1) * 0.05  # 0.05 units between targets
    target_centre <- clump_width / 2
    target_start  <- target_centre - target_width / 2
    layout_xy[target_order, 1] <- seq(target_start, target_start + target_width,
                                       length.out = n_targets)

    ## Create manual layout data frame for ggraph
    layout_df <- data.frame(
        x = layout_xy[, 1],
        y = layout_xy[, 2]
    )

    ## --- 7. Plot with ggraph ---
    p <- ggraph(g, layout = "manual", x = layout_df$x, y = layout_df$y) +
        geom_edge_link(colour = "grey60", width = 0.3, alpha = 0.4) +
        ## Clump nodes: black squares
        geom_node_point(
            aes(filter = is_clump),
            shape = 15, size = clump_node_size, colour = "black"
        ) +
        ## Target gene nodes: circles sized by -log10(p), coloured by direction
        geom_node_point(
            aes(filter = !is_clump, size = neglog10p, colour = direction),
            shape = 16
        ) +
        ## Clump labels (above nodes, rotated)
        geom_node_text(
            aes(filter = is_clump, label = display_label),
            size = 2.5, fontface = "bold", colour = "black",
            angle = 70, hjust = 0, nudge_y = clump_label_nudge_y
        ) +
        ## Target gene labels (below nodes, rotated)
        geom_node_text(
            aes(filter = !is_clump, label = display_label),
            size = label_size, fontface = "italic", colour = "grey20",
            angle = 70, hjust = 1, nudge_y = target_label_nudge_y
        ) +
        ## Edge count below target gene labels
        geom_node_text(
            aes(filter = !is_clump, label = n_edges),
            size = count_size, colour = "grey40",
            hjust = 0.5, nudge_y = count_nudge_y
        ) +
        scale_colour_manual(
            name = "Slope",
            values = c("Positive" = "steelblue", "Negative" = "firebrick"),
            na.translate = FALSE
        ) +
        scale_size_area(
            name = expression(-log[10](italic(P))),
            max_size = 8
        ) +
        scale_y_continuous(expand = expansion(mult = c(0.55, 0.45))) +
        scale_x_continuous(expand = expansion(mult = c(0.12, 0.12))) +
        guides(
            colour = guide_legend(title.position = "top"),
            size   = guide_legend(title.position = "top")
        ) +
        theme_void() +
        theme(legend.position = "bottom",
              legend.box = "horizontal",
              legend.margin = margin(t = -8, b = 2, l = 0, r = 0, "mm"),
              ## ordering is top, right, bottom, left
              plot.margin = margin(10, 1, 10, 1, "mm"),
              legend.text  = if (!is.null(legend_text_size))
                                 element_text(size = legend_text_size) else NULL,
              legend.title = if (!is.null(legend_text_size))
                                 element_text(size = legend_text_size) else NULL)

    return(p)
}
