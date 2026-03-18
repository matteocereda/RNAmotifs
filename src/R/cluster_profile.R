#!/usr/bin/env Rscript
# RNAmotifs - Cluster-averaged profile curves
# Copyright (C) 2014-2026 Matteo Cereda
# SPDX-License-Identifier: GPL-2.0-or-later
#
# For each tetramer cluster identified in the RNA splicing map, average
# position-specific scores across tetramers, normalise 0-1, and plot
# as curves in the same coordinate system as the RNA splicing map.

.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))

args <- commandArgs(TRUE)
if (length(args) < 4)
    stop("Usage: Rscript cluster_profile.R <profile_tsv> <enriched_csv> <output_pdf> <score_col> [title]")

profile_file  <- args[1]    # conservation_profile.tsv or structure_profile.tsv
enriched_csv  <- args[2]    # MRMs_<name>_*.csv
output_pdf    <- args[3]
score_col     <- args[4]    # "phylop_score" or "ss_score"
plot_title    <- if (length(args) >= 5) args[5] else "Cluster profile"
in_exon       <- if (length(args) >= 6) as.numeric(args[6]) else 30
in_intron     <- if (length(args) >= 7) as.numeric(args[7]) else 300

script_dir <- tryCatch(dirname(sys.frame(1)$ofile),
    error = function(e) {
        args_all <- commandArgs(FALSE)
        f <- grep("--file=", args_all, value = TRUE)
        if (length(f) > 0) dirname(sub("--file=", "", f[1]))
        else getwd()
    })
source(file.path(script_dir, "config.R"))

# --- Read data ----------------------------------------------------------------

prof <- read.delim(profile_file, stringsAsFactors = FALSE)
if (nrow(prof) == 0) { cat("No profile data.\n"); quit(save = "no", status = 0) }

enr <- read.csv(enriched_csv, stringsAsFactors = FALSE)

# Get unique tetramers with cluster assignments
tet_clusters <- unique(enr[, c("tetramer", "cluster_id")])
tet_clusters <- tet_clusters[!is.na(tet_clusters$cluster_id), ]

# Merge cluster info into profile data
prof <- merge(prof, tet_clusters, by = "tetramer")

# Keep only enhanced and silenced
prof <- prof[prof$category %in% c(1, -1), ]
if (nrow(prof) == 0) { cat("No data after filtering.\n"); quit(save = "no", status = 0) }

prof$cat_label <- ifelse(prof$category == 1, "Enhanced", "Silenced")

cat("Clusters:", length(unique(prof$cluster_id)), "\n")
cat("Tetramers:", length(unique(prof$tetramer)), "\n")

# --- Map to plot coordinates --------------------------------------------------

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

prof$reg <- floor((prof$position + intron) / regions)
prof$val <- prof$position - (prof$reg * regions)
prof$plot <- NA
for (r in c("1", "2", "3", "4")) {
    sel <- which(prof$reg == r & !is.na(prof$reg))
    if (length(sel) > 0)
        prof$plot[sel] <- prof$val[sel] + offset$plus[offset$reg == r]
}
prof <- prof[prof$position %in% pos_region & !is.na(prof$plot), ]

# Remove gap positions
gap_starts <- c(offset$end[1] + 1, offset$end[2] + 1, offset$end[3] + 1)
gap_ends   <- c(offset$start[2],   offset$start[3],   offset$start[4])
for (g in seq_along(gap_starts))
    prof <- prof[!(prof$plot >= gap_starts[g] & prof$plot <= gap_ends[g]), ]

# --- Average by cluster, position, category (un-normalised) ------------------

prof$score <- prof[[score_col]]
prof <- prof[prof$n_exons > 0, ]

agg <- aggregate(score ~ cluster_id + position + plot + cat_label,
                 data = prof, FUN = mean, na.rm = TRUE)

clusters <- sort(unique(agg$cluster_id))

agg$cluster_label <- paste0("Cluster ", agg$cluster_id)
agg$cluster_label <- factor(agg$cluster_label,
    levels = paste0("Cluster ", rev(clusters)))

# Exon-intron boundaries
exon_intron <- c(offset$plus[1], offset$plus[2] + 1,
                 offset$plus[3], offset$plus[4] + 1)

# Cluster colors (matching selection.R rainbow assignment)
cluster_cols <- rainbow(max(clusters))
names(cluster_cols) <- as.character(seq_len(max(clusters)))

# --- Build wide table for ribbon fill ----------------------------------------

cat_cols <- c("Enhanced" = "#e05555", "Silenced" = "#2166ac")
ribbon_cols <- c("Enhanced" = adjustcolor("#e05555", alpha.f = 0.4),
                 "Silenced" = adjustcolor("#2166ac", alpha.f = 0.4))
n_clusters <- length(clusters)

# Pivot to wide
agg_enh <- agg[agg$cat_label == "Enhanced",
               c("cluster_id", "cluster_label", "position", "plot", "score")]
agg_sil <- agg[agg$cat_label == "Silenced",
               c("cluster_id", "position", "score")]
names(agg_enh)[5] <- "enh"
names(agg_sil)[3] <- "sil"

wide <- merge(agg_enh, agg_sil, by = c("cluster_id", "position"), all = TRUE)
wide$enh[is.na(wide$enh)] <- NA
wide$sil[is.na(wide$sil)] <- NA

# Assign region
wide$region <- NA
wide$region[wide$plot >= offset$start[1] & wide$plot <= offset$end[1]] <- 1
wide$region[wide$plot >= offset$start[2] & wide$plot <= offset$end[2]] <- 2
wide$region[wide$plot >= offset$start[3] & wide$plot <= offset$end[3]] <- 3
wide$region[wide$plot >= offset$start[4] & wide$plot <= offset$end[4]] <- 4
wide <- wide[!is.na(wide$region) & !is.na(wide$plot), ]
wide <- wide[order(wide$cluster_id, wide$region, wide$plot), ]

# Smooth per cluster per region
smooth_vec <- function(x, k = 15) {
    n <- length(x)
    if (n < k) return(x)
    s <- as.numeric(stats::filter(x, rep(1/k, k), sides = 2))
    s[is.na(s)] <- x[is.na(s)]
    s
}

wide <- do.call(rbind, lapply(split(wide, list(wide$cluster_id, wide$region)), function(d) {
    if (nrow(d) < 3) return(d)
    d <- d[order(d$plot), ]
    d$enh_sm <- smooth_vec(d$enh)
    d$sil_sm <- smooth_vec(d$sil)
    d
}))
rownames(wide) <- NULL

wide_by_cluster <- split(wide, wide$cluster_label)

# Assign region to agg
agg$region <- NA
agg$region[agg$plot >= offset$start[1] & agg$plot <= offset$end[1]] <- 1
agg$region[agg$plot >= offset$start[2] & agg$plot <= offset$end[2]] <- 2
agg$region[agg$plot >= offset$start[3] & agg$plot <= offset$end[3]] <- 3
agg$region[agg$plot >= offset$start[4] & agg$plot <= offset$end[4]] <- 4
agg <- agg[!is.na(agg$region), ]

# Global y-range across all clusters
y_range <- range(agg$score, na.rm = TRUE)
y_pad <- diff(y_range) * 0.05
ylim <- c(y_range[1] - y_pad, y_range[2] + y_pad)

# Y-axis tick marks
y_at <- pretty(y_range, n = 3)
y_at <- y_at[y_at >= ylim[1] & y_at <= ylim[2]]

# --- Plot ---------------------------------------------------------------------

p <- xyplot(score ~ plot | cluster_label, data = agg,
    groups = cat_label,
    type = "l", lwd = 1.2,
    layout = c(1, n_clusters),
    aspect = "fill",
    xlim = c(-40, offset$end[4] + 10),
    ylim = ylim,
    between = list(y = 0),
    xlab = "", ylab = "",
    main = list(plot_title, cex = 0.85),
    strip = FALSE,
    strip.left = strip.custom(horizontal = TRUE,
                              bg = GREY_BG,
                              par.strip.text = list(cex = 0.8)),
    strip.left.lines = 0.5,
    auto.key = list(space = "top", columns = 2, lines = TRUE, points = FALSE,
                    cex = 0.65),
    par.settings = list(
        superpose.line = list(col = c(cat_cols["Enhanced"], cat_cols["Silenced"]),
                              lwd = 1.2),
        axis.line = list(lwd = 0.5, col = GREY_LINE),
        layout.widths = list(strip.left = 4),
        strip.border = list(col = GREY_LINE, lwd = 0.5)
    ),
    scales = list(
        x = list(draw = FALSE),
        y = list(at = y_at, cex = 0.5, tck = c(0, 1), alternating = 2)
    ),
    panel = function(x, y, subscripts, groups, ...) {
        # Cluster color strip on the left
        cl_label <- levels(agg$cluster_label)[packet.number()]
        cl_id <- sub("Cluster ", "", cl_label)
        panel.rect(-40, ylim[1], -30, ylim[2],
                   col = cluster_cols[cl_id], border = NA)

        # Exon shading
        panel.rect(offset$start[1], ylim[1], offset$start[1] + exon, ylim[2],
                   col = "gray96", border = NA)
        panel.rect(offset$end[2] - exon + 1, ylim[1], offset$end[2] + 1, ylim[2],
                   col = "gray96", border = NA)
        panel.rect(offset$start[3], ylim[1], offset$start[3] + exon, ylim[2],
                   col = "gray96", border = NA)
        panel.rect(offset$end[4] - exon + 1, ylim[1], offset$end[4] + 1, ylim[2],
                   col = "gray96", border = NA)

        # Get this panel's wide data for ribbon
        cl <- levels(agg$cluster_label)[packet.number()]
        wd <- wide_by_cluster[[cl]]

        # Draw ribbon fill between smoothed curves per region
        if (!is.null(wd) && nrow(wd) > 0) {
            for (ri in 1:4) {
                rd <- wd[wd$region == ri, ]
                if (nrow(rd) < 2) next
                rd <- rd[order(rd$plot), ]
                xx <- rd$plot; ye <- rd$enh_sm; ys <- rd$sil_sm
                for (i in seq_len(nrow(rd) - 1)) {
                    if (is.na(ye[i]) || is.na(ye[i+1]) ||
                        is.na(ys[i]) || is.na(ys[i+1])) next
                    seg_x <- c(xx[i], xx[i+1], xx[i+1], xx[i])
                    seg_y <- c(ye[i], ye[i+1], ys[i+1], ys[i])
                    mid_e <- (ye[i] + ye[i+1]) / 2
                    mid_s <- (ys[i] + ys[i+1]) / 2
                    fill <- if (mid_e >= mid_s) ribbon_cols["Enhanced"]
                            else ribbon_cols["Silenced"]
                    panel.polygon(seg_x, seg_y, col = fill, border = NA)
                }
            }
        }

        # Draw smoothed lines per region
        if (!is.null(wd) && nrow(wd) > 0) {
            for (ri in 1:4) {
                rd <- wd[wd$region == ri, ]
                if (nrow(rd) < 2) next
                rd <- rd[order(rd$plot), ]
                ok <- !is.na(rd$enh_sm) & !is.na(rd$sil_sm)
                if (sum(ok) > 1) {
                    panel.lines(rd$plot[ok], rd$enh_sm[ok],
                                col = cat_cols["Enhanced"], lwd = 1.5)
                    panel.lines(rd$plot[ok], rd$sil_sm[ok],
                                col = cat_cols["Silenced"], lwd = 1.5)
                }
            }
        }

        # White gaps
        for (g in seq_along(gap_starts))
            panel.rect(gap_starts[g], ylim[1] - 1, gap_ends[g], ylim[2] + 1,
                       col = "white", border = NA)

        # Exon-intron boundaries
        panel.abline(v = exon_intron, lty = 1, col = "black", lwd = 0.5)
        panel.abline(v = c(offset$start, offset$end + 1),
                     lty = "dotted", col = GREY_LINE, lwd = 0.3)
    }
)

pdf(file = output_pdf, width = 10, height = 10)
print(p, panel.height = list(1.5, "cm"), panel.width = list(10, "cm"))
dev.off()

cat("Cluster profile PDF:", output_pdf, "\n")
