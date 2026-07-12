#!/usr/bin/env Rscript
# RNAmotifs - Tetramer selection and RNA splicing map visualisation
# Copyright (C) 2014-2026 Matteo Cereda
# SPDX-License-Identifier: GPL-2.0-or-later

.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))

args <- commandArgs(TRUE)
if (length(args) < 5)
    stop("Usage: Rscript selection.R <results_dir> <name> <n_bootstraps> <p_fisher> <p_empirical> [top_n] [counts_dir]")

wd      <- args[1]
name    <- args[2]
n_boot  <- as.numeric(args[3])
cFisher <- as.numeric(args[4])
cEmp    <- as.numeric(args[5])
top_n   <- if (length(args) >= 6) as.numeric(args[6]) else 0
counts_dir <- if (length(args) >= 7) args[7] else ""
in_exon    <- if (length(args) >= 8) as.numeric(args[8]) else 30
in_intron  <- if (length(args) >= 9) as.numeric(args[9]) else 300

pp <- "./"
# If counts_dir is provided, use it for filelists; otherwise use pp
cp <- if (nchar(counts_dir) > 0) paste0(counts_dir, "/") else pp

script_dir <- tryCatch(dirname(sys.frame(1)$ofile),
    error = function(e) {
        args_all <- commandArgs(FALSE)
        f <- grep("--file=", args_all, value = TRUE)
        if (length(f) > 0) dirname(sub("--file=", "", f[1]))
        else getwd()
    })
source(file.path(script_dir, "config.R"))

setwd(wd)

wr   <- c("r/", "nr/")
# Compute region centers: spacing must exceed 2*max(in_exon,in_intron)
# so adjacent regions don't overlap in position space
region_spacing <- 2 * max(in_exon, in_intron) + 10  # 10 = gap
regions <- cumsum(c(region_spacing, rep(region_spacing, 3)))
ROWN <- c((regions[1] - in_exon):(regions[1] + in_intron),
          (regions[2] - in_intron):(regions[2] + in_exon),
          (regions[3] - in_exon):(regions[3] + in_intron),
          (regions[4] - in_intron):(regions[4] + in_exon))
rown <- ROWN

cat("\n[selection] Loading bootstrap results...\n")
res <- read.delim(paste0(pp, "bootstrap_", n_boot, ".tsv"),
                  stringsAsFactors = FALSE)

# --- Select significantly enriched motifs -------------------------------------

# Adaptive per-region Fisher cutoff: the 1st percentile of each region's
# combined (enhanced + silenced) bootstrap Fisher p-values, capped at 0.05.
# This matches the RNAMaRs discovery selector (mars/selection_of_tetramers.R)
# so the plain CLI and the mars pipeline select identically. cFisher is retained
# only for the per-region where.sign annotation below, NOT for selection.
p_cutoff <- sapply(1:3, function(r) {
    q <- quantile(c(res[[paste0("r", r, "enh_pFis")]],
                    res[[paste0("r", r, "sil_pFis")]]), 0.01, names = FALSE)
    if (q < 0.05) q else 0.05
})
cat(sprintf("[selection] adaptive Fisher cutoffs (1st pct, cap 0.05): %s\n",
            paste(signif(p_cutoff, 3), collapse = ", ")))

## Methods: keep pFis <= min(1st-percentile, 0.05). Use <= (not <): '<' wrongly
## dropped the most significant motifs when the cutoff collapses to 0 at wide
## windows (>=1% of motifs saturate the Fisher test at pFis==0).
sig <- subset(res,
    (r1enh_pFis <= p_cutoff[1] & r1enh_pEmp <= cEmp) | (r1sil_pFis <= p_cutoff[1] & r1sil_pEmp <= cEmp) |
    (r2enh_pFis <= p_cutoff[2] & r2enh_pEmp <= cEmp) | (r2sil_pFis <= p_cutoff[2] & r2sil_pEmp <= cEmp) |
    (r3enh_pFis <= p_cutoff[3] & r3enh_pEmp <= cEmp) | (r3sil_pFis <= p_cutoff[3] & r3sil_pEmp <= cEmp))

if (nrow(sig) == 0) {
    cat("No significant tetramers found.\n")
    quit(save = "no", status = 0)
}

cat("Significantly enriched motifs:", nrow(sig), "\n")

sig[, 1] <- as.character(sig[, 1])

# Build summary table with enhanced/silenced rows
s <- matrix(0, nrow = 2 * nrow(sig), ncol = 11,
            dimnames = list(NULL,
                c("tetramer", "full_motifs", "exonType", "is.sign",
                  "where.sign", "r1_pf", "r1_pe", "r2_pf", "r2_pe",
                  "r3_pf", "r3_pe")))
s[, 1] <- rep(sig[, 1], 2)
s[, 3] <- c(rep("enh", nrow(sig)), rep("sil", nrow(sig)))

s[seq_len(nrow(sig)),
  c("r1_pf", "r1_pe", "r2_pf", "r2_pe", "r3_pf", "r3_pe")] <-
    as.matrix(sig[, c("r1enh_pFis", "r1enh_pEmp", "r2enh_pFis",
                       "r2enh_pEmp", "r3enh_pFis", "r3enh_pEmp")])

s[(nrow(sig) + 1):nrow(s),
  c("r1_pf", "r1_pe", "r2_pf", "r2_pe", "r3_pf", "r3_pe")] <-
    as.matrix(sig[, c("r1sil_pFis", "r1sil_pEmp", "r2sil_pFis",
                       "r2sil_pEmp", "r3sil_pFis", "r3sil_pEmp")])

s <- as.data.frame(s, stringsAsFactors = FALSE)
s[, 4:ncol(s)] <- lapply(s[, 4:ncol(s)], as.numeric)

x1 <- apply(s[, c("r1_pe", "r2_pe", "r3_pe")], 1, function(x) x <= cEmp)
x2 <- apply(s[, c("r1_pf", "r2_pf", "r3_pf")], 1, function(x) x <= cFisher)

s$where.sign <- apply(x1 + x2, 2, function(x) paste(which(x == 2), collapse = ","))
s$is.sign    <- s$where.sign != ""

tmp <- strsplit(s[, 1], "")
s$full_motifs <- sapply(tmp, function(x) paste(IUPAC_CODE[x], collapse = ""))

tets <- unique(sig[, 1])

# --- Exon counts --------------------------------------------------------------

flist_path <- paste0(cp, wr[1], "filelist_count.tsv")
fn <- read.delim(paste0(cp, wr[1],
    read.table(flist_path, stringsAsFactors = FALSE, header = FALSE)[1, 1]))
et <- table(fn$type)
CEone    <- et["1"]
CEminone <- et["-1"]
CEzero   <- et["0"]

# --- Positional Fisher tests --------------------------------------------------

cat("Calculating positional Fisher's test...\n")
p.ena <- p.sil <- NULL
for (f in 1:2) {
    ll <- getTables(cp, wr[f], tets)
    if (length(ll) > 0) {
        p.ena <- cbind(p.ena,
            lmb.cluster.fisher.test(ll[["enh"]], CEone, ll[["cont"]], CEzero))
        p.sil <- cbind(p.sil,
            lmb.cluster.fisher.test(ll[["sil"]], CEminone, ll[["cont"]], CEzero))
    }
}

# --- Fisher's method: combined score ------------------------------------------

cat("Calculating combined Fisher's method score...\n")
ES <- matrix(0, ncol = length(tets), nrow = nrow(p.ena),
             dimnames = list(rownames(p.ena), tets))
for (i in tets) ES[, i] <- (-2) * (log(p.ena[, i]) + log(p.sil[, i]))

# --- Sort and cluster tetramers -----------------------------------------------

cat("Sorting and clustering tetramers...\n")
auc  <- apply(ES, 2, sum)
tets <- names(sort(auc, decreasing = TRUE))

score.sort <- (-2) * log(p.ena) - (-2) * log(p.sil)

st <- if (length(tets) == 1) list(tets, 1) else SortingAndType(score.sort, tets)
ord  <- st[[1]]
rcol <- rainbow(max(st[[2]]))[st[[2]]]
names(rcol) <- ord

df_ord <- do.call(cbind, st)
s$cluster_id <- as.numeric(df_ord[match(s$tetramer, df_ord[, 1]), 2])

# Write tetramer order file (ord = splicing map display order, with cluster ids)
order_df <- data.frame(tetramer = ord, cluster_id = st[[2]], stringsAsFactors = FALSE)
write.table(order_df, file = paste0(pp, "tetramer_order.txt"),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
cat("Tetramer order (splicing map):", length(ord), "tetramers\n")

# Save enriched tetramers table
san <- function(x) gsub("\\.", "_", as.character(x))
out_csv <- paste0(pp, "MRMs_", name, "_emp-", san(cEmp),
                  "_fisher-", san(cFisher), "_nBoot-", n_boot, ".csv")
write.csv(s, file = out_csv, row.names = FALSE)
cat("Enriched tetramers table:", out_csv, "\n")

# --- Apply top-N filter if requested ------------------------------------------

if (top_n > 0 && top_n < length(ord)) {
    cat("Filtering to top", top_n, "ranked tetramers...\n")
    ord  <- ord[1:top_n]
    rcol <- rcol[1:top_n]
    names(rcol) <- ord
}

# --- RNA splicing maps --------------------------------------------------------

cat("Plotting RNA splicing maps...\n")

e_score <- (-2) * log(p.ena)
s_score <- (-2) * log(p.sil)
datar <- NULL
for (i in ord)
    datar <- rbind(datar,
        data.frame(s = ES[, i], e = e_score[, i], l = rown,
                   g = rep(i, nrow(ES)), stringsAsFactors = FALSE))

ms       <- max(datar$s)
datar$z  <- datar$e / datar$s
datar$g  <- factor(datar$g, rev(ord))
datar$s  <- datar$s / ms
datar$z  <- datar$z / ms

prn <- rnaScoreMap(datar, rcol, ylabels = c("", as.character(floor(ms))),
                   exon = in_exon, intron = in_intron,
                   regions = region_spacing)
prn <- update(prn, lwd = 0.8)

pdf_path <- paste0(pp, "MRMs_", name, "_emp-", san(cEmp),
                   "_fisher-", san(cFisher), "_nBoot-", n_boot, ".pdf")
pdf(file = pdf_path, height = 10)
print(prn, panel.height = list(0.4, "cm"), panel.width = list(10, "cm"))
dev.off()

rdata_path <- sub("\\.pdf$", ".Rdata", pdf_path)
df_enriched_tetramers <- s
save(prn, df_enriched_tetramers, file = rdata_path)

cat("RNA map PDF:", pdf_path, "\n")
cat("Done.\n")
