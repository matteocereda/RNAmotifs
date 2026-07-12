#!/usr/bin/env Rscript
# IR Splicing Map — uses rnaScoreMap style from config.R
# Regions: R1 (upstream exon) | R2 (5' retained intron) | R3 (3' retained intron) | R4 (downstream exon)
#
# Usage: Rscript ir_splicing_map.R <results_dir> <output_pdf> [in_exon] [in_intron]

.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))

args <- commandArgs(TRUE)
if (length(args) < 2) stop("Usage: ir_splicing_map.R <results_dir> <output_pdf> [in_exon] [in_intron]")

rdir      <- args[1]
outfile   <- args[2]
in_exon   <- if (length(args) >= 3) as.numeric(args[3]) else 30
in_intron <- if (length(args) >= 4) as.numeric(args[4]) else 300

# Source config.R to get myPalette, panel.rnaplot, GREY_BG, GREY_LINE etc
script_dir <- tryCatch(dirname(sys.frame(1)$ofile),
    error = function(e) {
        args_all <- commandArgs(FALSE)
        f <- grep("--file=", args_all, value = TRUE)
        if (length(f) > 0) dirname(sub("--file=", "", f[1]))
        else getwd()
    })
source(file.path(script_dir, "config.R"))

exon <- in_exon; intron <- in_intron; gap <- 10
regions <- 2 * max(exon, intron) + gap

# --- IR-specific rnaScoreMap ---
# Same as rnaScoreMap but with IR offset table and exon shading
rnaScoreMap_IR <- function(data, rcol, ylabels = c("0.5", "1"), main = "",
                           exon = 30, intron = 300, gap = 4, regions = 610) {

    # IR offset: contiguous layout
    # [R1:exon] gap [R2:intron] gap [R3:intron] gap [R4:exon]

    r1_w <- exon; r2_w <- intron; r3_w <- intron; r4_w <- exon

    p1 <- 0
    p2 <- r1_w + gap
    p3 <- r1_w + gap + r2_w + gap
    p4 <- r1_w + gap + r2_w + gap + r3_w + gap

    offset <- data.frame(
        reg   = c("1", "2", "3", "4"),
        start = c(p1, p2, p3, p4),
        plus  = c(p1, p2, p3, p4),
        end   = c(p1 + r1_w, p2 + r2_w, p3 + r3_w, p4 + r4_w))

    total_width <- 2*exon + 2*intron + 3*gap

    # Map abstract positions to plot coords
    # C++ maps IR positions to: R1=[regions-exon, regions], R2=[2*regions, 2*regions+intron],
    # R3=[3*regions-intron, 3*regions], R4=[4*regions, 4*regions+exon]
    r1_lo <- regions - exon;    r1_hi <- regions
    r2_lo <- 2*regions;         r2_hi <- 2*regions + intron
    r3_lo <- 3*regions - intron; r3_hi <- 3*regions
    r4_lo <- 4*regions;         r4_hi <- 4*regions + exon

    rownames(data) <- NULL
    pf <- cbind(data, plot = NA)
    pf <- as.data.frame(pf)

    # Map each abstract position l to plot coordinate
    sel <- pf$l >= r1_lo & pf$l <= r1_hi
    pf$plot[sel] <- (pf$l[sel] - r1_lo) / max(1, r1_hi - r1_lo) * (offset$end[1] - offset$start[1]) + offset$start[1]

    sel <- pf$l >= r2_lo & pf$l <= r2_hi
    pf$plot[sel] <- (pf$l[sel] - r2_lo) / max(1, r2_hi - r2_lo) * (offset$end[2] - offset$start[2]) + offset$start[2]

    sel <- pf$l >= r3_lo & pf$l <= r3_hi
    pf$plot[sel] <- (pf$l[sel] - r3_lo) / max(1, r3_hi - r3_lo) * (offset$end[3] - offset$start[3]) + offset$start[3]

    sel <- pf$l >= r4_lo & pf$l <= r4_hi
    pf$plot[sel] <- (pf$l[sel] - r4_lo) / max(1, r4_hi - r4_lo) * (offset$end[4] - offset$start[4]) + offset$start[4]

    pf_real <- pf[!is.na(pf$plot), ]
    if (nrow(pf_real) == 0) { cat("No data after position mapping\n"); return(NULL) }

    ylim <- c(0, max(pf_real$s, na.rm=TRUE))
    pf_real$rcol <- rcol[as.character(pf_real$g)]

    # Key boundary positions
    ss5 <- (offset$end[1] + offset$start[2]) / 2   # 5'SS: between R1 exon and R2 intron
    ss3 <- (offset$end[3] + offset$start[4]) / 2   # 3'SS: between R3 intron and R4 exon

    # Tick positions and labels for top axis
    tick_at <- c(ss5, ss3)
    tick_labels <- c("5'SS", "3'SS")

    xyplot(s ~ plot | g, data = pf_real,
           intensity = pf_real$z,
           rcol = pf_real$rcol,
           hrect = ylim[2],
           aspect = "fill", type = "l",
           layout = c(1, length(levels(data$g))),
           ylim = ylim,
           xlim = c(-30, total_width + 10),
           subscripts = TRUE,
           strip = FALSE,
           strip.left = strip.custom(horizontal = TRUE,
                                     var.name = unique(pf_real$g),
                                     bg = GREY_BG),
           strip.left.lines = 0.5,
           panel = function(x, y, hrect, rcol, intensity, groups,
                            subscripts, ...) {
               # Shade R1 (upstream exon) and R4 (downstream exon)
               panel.rect(offset$start[1], -1, offset$end[1], hrect*1.1,
                          col = "gray96", border = FALSE)
               panel.rect(offset$start[4], -1, offset$end[4], hrect*1.1,
                          col = "gray96", border = FALSE)

               # White gaps at splice sites and intron midpoint
               panel.rect(offset$end[1], -1, offset$start[2], hrect*1.1,
                          col = "white", border = FALSE)
               panel.rect(offset$end[2], -1, offset$start[3], hrect*1.1,
                          col = "white", border = FALSE)
               panel.rect(offset$end[3], -1, offset$start[4], hrect*1.1,
                          col = "white", border = FALSE)

               # Draw the colored bars (same as ES panel.rnaplot)
               panel.rnaplot(x, y, ry = hrect,
                             rc = rcol[subscripts],
                             intensity = intensity[subscripts], ...)

               # Boundary lines: dotted at all region edges
               panel.abline(v = c(offset$start, offset$end + 1),
                            lty = "dotted", col = "black", lwd = 0.4)
               # Splice site lines: solid at exon/intron junctions
               panel.abline(v = c(ss5, ss3), lty = 1, col = "black", lwd = 0.6)
           },
           scales = list(
               x = list(at = tick_at,
                        labels = tick_labels,
                        tck = c(0, 1),
                        cex = 0.45,
                        alternating = 2),
               y = list(tck = c(0, 1),
                        at = c(ylim[2] / 2, ylim[2]),
                        tick.number = 2,
                        labels = ylabels,
                        alternating = 2),
               cex = 0.5),
           ylab = "", xlab = "", main = main,
           par.strip.text = list(cex = 0.8),
           par.settings = list(
               axis.line     = list(lwd = 0.5, col = GREY_LINE),
               layout.widths = list(key.ylab.padding = 1, strip.left = 4),
               strip.border  = list(col = GREY_LINE, lwd = 0.5),
               plot.line     = list(col = GREY_LINE),
               add.line      = list(lwd = 0.5)),
           legend = list(left = list(
               fun = draw.colorkey,
               args = list(key = list(
                   col = myPalette(100), at = 0:100, tick.number = 3,
                   labels = list(labels = c("100% S", "50% E,S", "100% E"),
                                 at = c(0, 50, 100),
                                 cex = 0.5),
                   raster = TRUE, width = 1, height = 0.5,
                   space = "left"),
                   draw = FALSE))))
}

# --- Read top 10 tetramers ---
top10 <- readLines(file.path(rdir, "enriched_tetramers_top10.txt"))
top10 <- top10[top10 != ""]
cat("Tetramers:", paste(top10, collapse=", "), "\n")

# --- Read bootstrap data and compute Fisher scores per tetramer ---
bs_file <- list.files(rdir, pattern="^bootstrap_.*\\.tsv$", full.names=TRUE)[1]
if (is.na(bs_file)) stop("No bootstrap TSV found in ", rdir)

bs <- read.delim(bs_file, stringsAsFactors=FALSE)

# Position scoring: read e- and s- CSVs if they exist; otherwise compute from BED
e_csv <- list.files(rdir, pattern="^e-.*\\.csv$", full.names=TRUE)
s_csv <- list.files(rdir, pattern="^s-.*\\.csv$", full.names=TRUE)

if (length(e_csv) > 0 && length(s_csv) > 0) {
    cat("Reading score CSVs...\n")
    e_mat <- as.matrix(read.csv(e_csv[1], check.names=FALSE))
    s_mat <- as.matrix(read.csv(s_csv[1], check.names=FALSE))
    # Strip quotes from column names
    colnames(e_mat) <- gsub('"', '', colnames(e_mat))
    colnames(s_mat) <- gsub('"', '', colnames(s_mat))

    ES <- e_mat + s_mat
    rown <- 1:nrow(ES)

    # Build datar in the rnaScoreMap format
    ord <- top10[top10 %in% colnames(ES)]
    if (length(ord) == 0) {
        cat("No top10 tetramers found in score CSVs, trying BED fallback\n")
        ord <- NULL
    }
} else {
    ord <- NULL
}

if (is.null(ord) || length(ord) == 0) {
    # Fallback: compute from BED files
    cat("Computing scores from BED files...\n")
    all_datar <- NULL

    for (tet in top10) {
        bed_path <- file.path(rdir, "data", "counts", "nr", paste0(tet, ".bed"))
        if (!file.exists(bed_path))
            bed_path <- file.path(rdir, "data", "counts", "r", paste0(tet, ".bed"))
        if (!file.exists(bed_path)) next

        bed <- read.table(bed_path, sep="\t")
        colnames(bed) <- c("position", "cat", "count")

        # Aggregate by position
        agg_e <- aggregate(count ~ position, data=bed[bed$cat == 1,], FUN=sum)
        agg_s <- aggregate(count ~ position, data=bed[bed$cat == -1,], FUN=sum)
        all_pos <- sort(unique(bed$position))

        enh_v <- rep(0, length(all_pos)); names(enh_v) <- all_pos
        sil_v <- rep(0, length(all_pos)); names(sil_v) <- all_pos
        if (nrow(agg_e) > 0) enh_v[as.character(agg_e$position)] <- agg_e$count
        if (nrow(agg_s) > 0) sil_v[as.character(agg_s$position)] <- agg_s$count

        ES_v <- enh_v + sil_v
        z_v <- ifelse(ES_v > 0, enh_v / ES_v, 0.5)

        all_datar <- rbind(all_datar, data.frame(
            s = ES_v, e = enh_v, l = all_pos, g = tet, z = z_v,
            stringsAsFactors = FALSE
        ))
    }

    if (is.null(all_datar) || nrow(all_datar) == 0) {
        cat("No data to plot\n"); quit(save="no", status=1)
    }

    ms <- max(all_datar$s, na.rm=TRUE)
    all_datar$s <- all_datar$s / ms
    all_datar$z <- all_datar$z / ms * ms  # keep ratio
    # Recompute z as e/s ratio (0-1)
    all_datar$z <- ifelse(all_datar$e + (all_datar$s * ms - all_datar$e) > 0,
                          all_datar$e / (all_datar$s * ms), 0.5)
    all_datar$g <- factor(all_datar$g, levels = rev(top10))

    ord <- top10
    datar <- all_datar
} else {
    # Use score CSVs
    datar <- NULL
    for (i in ord) {
        datar <- rbind(datar, cbind(
            "s" = ES[,i], "e" = e_mat[,i], "l" = rown, "g" = rep(i, nrow(ES))
        ))
    }
    datar <- as.data.frame(datar, stringsAsFactors=FALSE)
    rownames(datar) <- NULL
    datar[, c("s","e","l")] <- apply(datar[, c("s","e","l")], 2, as.numeric)

    ms <- max(datar$s)
    datar$z <- datar$e / datar$s
    datar$g <- factor(datar$g, rev(ord))
    datar$s <- datar$s / ms
    datar$z <- datar$z / ms
}

# Cluster colors (rainbow for each tetramer)
rcol <- rainbow(length(ord))
names(rcol) <- ord

cat("Generating IR splicing map with", length(ord), "tetramers...\n")

prn <- rnaScoreMap_IR(datar, rcol,
                      ylabels = c("", as.character(floor(ms))),
                      main = "PTBP1 Intron Retention",
                      exon = in_exon, intron = in_intron,
                      gap = 4, regions = regions)

if (!is.null(prn)) {
    prn <- update(prn, lwd = 0.8)
    pdf(file = outfile, height = 10, width = 10)
    print(prn, panel.height = list(0.4, "cm"), panel.width = list(20, "cm"))
    dev.off()
    cat("Wrote:", outfile, "\n")
} else {
    cat("No plot generated\n")
}
