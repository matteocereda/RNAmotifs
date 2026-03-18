#!/usr/bin/env Rscript
# RNAmotifs - Sequence conservation profile visualisation (PhyloP)
# Copyright (C) 2014-2026 Matteo Cereda
# SPDX-License-Identifier: GPL-2.0-or-later

.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))

args <- commandArgs(TRUE)
if (length(args) < 4)
    stop("Usage: Rscript conservation_profile.R <conservation_tsv> <output_pdf> <tetramer_order_file> <enriched_csv>")

cons_file      <- args[1]
output_pdf     <- args[2]
order_file     <- args[3]
enriched_csv   <- args[4]
in_exon        <- if (length(args) >= 5) as.numeric(args[5]) else 30
in_intron      <- if (length(args) >= 6) as.numeric(args[6]) else 300

script_dir <- tryCatch(dirname(sys.frame(1)$ofile),
    error = function(e) {
        args_all <- commandArgs(FALSE)
        f <- grep("--file=", args_all, value = TRUE)
        if (length(f) > 0) dirname(sub("--file=", "", f[1]))
        else getwd()
    })
source(file.path(script_dir, "config.R"))

d <- read.delim(cons_file, stringsAsFactors = FALSE)
if (nrow(d) == 0) { cat("No data.\n"); quit(save = "no", status = 0) }

# Read tetramer order with cluster assignments (written by selection.R)
order_df <- read.delim(order_file, stringsAsFactors = FALSE)
if ("cluster_id" %in% names(order_df)) {
    tet_order <- order_df$tetramer
    cluster_ids <- order_df$cluster_id
} else {
    tet_order <- readLines(order_file)
    tet_order <- tet_order[tet_order != ""]
    enr <- read.csv(enriched_csv, stringsAsFactors = FALSE)
    tc <- unique(enr[, c("tetramer", "cluster_id")])
    tc <- tc[!is.na(tc$cluster_id), ]
    cluster_ids <- tc$cluster_id[match(tet_order, tc$tetramer)]
    cluster_ids[is.na(cluster_ids)] <- 0
}
tet_order <- tet_order[tet_order %in% unique(d$tetramer)]
cluster_ids <- cluster_ids[seq_along(tet_order)]
if (length(tet_order) == 0) { cat("No matching tetramers.\n"); quit(save = "no", status = 0) }
cluster_colors <- rainbow(max(cluster_ids))[cluster_ids]

cat("Plotting conservation heatmap for", length(tet_order), "tetramers...\n")

# Splicing map geometry
exon <- in_exon; intron <- in_intron; gap <- 10; regions <- 2 * max(exon, intron) + gap
offset <- data.frame(
    reg   = c("1", "2", "3", "4"),
    start = c(0, exon + intron + gap,
              2 * exon + 2 * intron + 2 * gap,
              3 * exon + 3 * intron + 3 * gap),
    plus  = c(exon, exon + 2 * intron + gap,
              3 * exon + 2 * intron + 2 * gap,
              4 * intron + 3 * exon + 3 * gap),
    end   = c(exon + intron, 2 * exon + 2 * intron + gap,
              3 * exon + 3 * intron + 2 * gap,
              4 * exon + 4 * intron + 3 * gap))

pos_region <- c((regions - exon):(regions + intron),
                (2 * regions - intron):(2 * regions + exon),
                (3 * regions - exon):(3 * regions + intron),
                (4 * regions - intron):(4 * regions + exon))

d$reg <- floor((d$position + intron) / regions)
d$val <- d$position - (d$reg * regions)
d$plot <- NA
for (r in c("1", "2", "3", "4")) {
    sel <- which(d$reg == r & !is.na(d$reg))
    if (length(sel) > 0)
        d$plot[sel] <- d$val[sel] + offset$plus[offset$reg == r]
}
d <- d[d$position %in% pos_region & !is.na(d$plot), ]

d$cat_label <- factor(d$category, levels = c(1, -1),
                      labels = c("Enhanced", "Silenced"))
d <- d[d$category %in% c(1, -1), ]

# Factor tetramers in RNA map order (top-ranked at top, matching rev(ord) in selection.R)
d$tetramer <- factor(d$tetramer, levels = rev(tet_order))
d <- d[!is.na(d$tetramer), ]

# Gaps
gap_starts <- c(offset$end[1] + 1, offset$end[2] + 1, offset$end[3] + 1)
gap_ends   <- c(offset$start[2],   offset$start[3],   offset$start[4])
for (g in seq_along(gap_starts))
    d <- d[!(d$plot >= gap_starts[g] & d$plot <= gap_ends[g]), ]

exon_intron_boundaries <- c(offset$plus[1], offset$plus[2] + 1,
                            offset$plus[3], offset$plus[4] + 1)

# Data-driven colour range
valid <- d$phylop_score[d$n_exons > 0]
sc_min <- min(valid, na.rm = TRUE)
sc_max <- max(valid, na.rm = TRUE)
sc_mid <- (sc_min + sc_max) / 2
cat("  PhyloP range:", round(sc_min, 3), "to", round(sc_max, 3), "\n")

# Inferno-like palette (dark -> hot)
consPalette <- colorRampPalette(c("#000004", "#1b0c41", "#4a0c6b",
                                   "#781c6d", "#a52c60", "#cf4446",
                                   "#ed6925", "#fb9b06", "#f7d13d",
                                   "#fcffa4"))

n_tets <- length(tet_order)
n_levels <- length(levels(d$tetramer))

p <- levelplot(phylop_score ~ plot * tetramer | cat_label, data = d,
    col.regions = consPalette(100),
    at = seq(sc_min, sc_max, length.out = 101),
    aspect = "fill",
    layout = c(2, 1),
    xlim = c(-15, offset$end[4] + 5),
    between = list(x = 1),
    xlab = "", ylab = "",
    main = list("Sequence conservation (PhyloP)", cex = 0.9),
    colorkey = list(
        space = "left", width = 1, height = 0.5,
        labels = list(
            at = c(sc_min, sc_mid, sc_max),
            labels = c(sprintf("%.2f", sc_min),
                       sprintf("%.2f", sc_mid),
                       sprintf("%.2f", sc_max)),
            cex = 0.55
        )
    ),
    strip = strip.custom(bg = GREY_BG, par.strip.text = list(cex = 0.85)),
    scales = list(x = list(draw = FALSE),
                  y = list(cex = 0.55, tck = c(0, 0))),
    par.settings = list(
        axis.line = list(lwd = 0.5, col = GREY_LINE),
        layout.widths = list(key.ylab.padding = 1, strip.left = 0),
        strip.border = list(col = GREY_LINE, lwd = 0.5)),
    panel = function(x, y, z, subscripts, ...) {
        # Cluster annotation strip (reversed to match factor order)
        rev_colors <- rev(cluster_colors)
        for (i in seq_along(tet_order)) {
            panel.rect(-15, i - 0.5, -5, i + 0.5,
                       col = rev_colors[i], border = NA)
        }
        # Exon shading
        panel.rect(offset$start[1], 0.5, offset$start[1] + exon, n_levels + 0.5,
                   col = "gray95", border = NA)
        panel.rect(offset$end[2] - exon + 1, 0.5, offset$end[2] + 1, n_levels + 0.5,
                   col = "gray95", border = NA)
        panel.rect(offset$start[3], 0.5, offset$start[3] + exon, n_levels + 0.5,
                   col = "gray95", border = NA)
        panel.rect(offset$end[4] - exon + 1, 0.5, offset$end[4] + 1, n_levels + 0.5,
                   col = "gray95", border = NA)
        # Heatmap
        panel.levelplot(x, y, z, subscripts, ...)
        # White gaps
        for (g in seq_along(gap_starts))
            panel.rect(gap_starts[g], 0, gap_ends[g], n_levels + 1,
                       col = "white", border = NA)
        # Boundaries
        panel.abline(v = exon_intron_boundaries, lty = 1, col = "black", lwd = 0.6)
        panel.abline(v = c(offset$start, offset$end + 1),
                     lty = "dotted", col = GREY_LINE, lwd = 0.4)
    }
)

pdf(file = output_pdf, width = 12, height = 10)
print(p, panel.height = list(n_tets * 0.4, "cm"), panel.width = list(10, "cm"))
dev.off()
cat("Conservation profile PDF:", output_pdf, "\n")
