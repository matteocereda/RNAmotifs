# RNAmotifs - R configuration and utility functions
# Copyright (C) 2014-2026 Matteo Cereda
# SPDX-License-Identifier: GPL-2.0-or-later

# --- Package management ------------------------------------------------------

.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))

required_packages <- c("lattice", "latticeExtra", "bootstrap", "ggplot2",
                        "reshape2", "parallel")
for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE))
        install.packages(pkg, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE, quietly = TRUE)
}

options(warn = -1, stringsAsFactors = FALSE)

# --- Constants ----------------------------------------------------------------

IUPAC_CODE <- c("A", "C", "T", "G", "(A/G)", "(C/T)", "(G/C)", "(A/T)",
                "(G/T)", "(A/C)", "(C/G/T)", "(A/G/T)", "(A/C/T)", "(A/C/G)")
names(IUPAC_CODE) <- c("A", "C", "T", "G", "R", "Y", "S", "W",
                        "K", "M", "B", "D", "H", "V")

ROWN <- c(950:1200, 1800:2050, 2950:3200, 3800:4050)

myPalette <- colorRampPalette(c("blue2", "yellow", "red2"))
GREY_BG   <- "grey95"
GREY_LINE <- rgb(94, 90, 90, maxColorValue = 255)

# --- Progress bar lapply ------------------------------------------------------

lapply_pb <- function(X, FUN, ...) {
    pb <- txtProgressBar(min = 0, max = length(X), style = 3)
    env <- environment()
    counter <- 0
    wrapper <- function(...) {
        counter <<- counter + 1
        setTxtProgressBar(pb, counter)
        FUN(...)
    }
    res <- lapply(X, wrapper, ...)
    close(pb)
    res
}

# --- Fisher test per tetramer ------------------------------------------------

get_fisher <- function(dd, exon) {
    out_p <- numeric(0)
    sel_e <- dd$type == 1
    sel_c <- dd$type == 0
    sel_s <- dd$type == -1

    for (reg in c("hits_region1", "hits_region2", "hits_region3")) {
        sel_p <- dd[exon, reg] > 0
        sel_a <- dd[exon, reg] == 0

        cont_enh <- matrix(c(sum(sel_e & sel_p), sum(sel_c & sel_p),
                              sum(sel_e & sel_a), sum(sel_c & sel_a)), nrow = 2)
        cont_sil <- matrix(c(sum(sel_s & sel_p), sum(sel_c & sel_p),
                              sum(sel_s & sel_a), sum(sel_c & sel_a)), nrow = 2)

        out_p <- c(out_p,
                   fisher.test(cont_enh, alternative = "greater")$p.value,
                   fisher.test(cont_sil, alternative = "greater")$p.value)
    }
    out_p
}

fisher.region <- function(exon, myHitsData) {
    n <- length(myHitsData)
    p_mat <- matrix(ncol = 6, nrow = n)
    out <- lapply_pb(myHitsData, get_fisher, exon = exon)
    for (i in seq_len(n)) p_mat[i, ] <- out[[i]]
    p_adj <- apply(p_mat, 2, p.adjust, method = "BH")
    rownames(p_adj) <- names(myHitsData)
    colnames(p_adj) <- c("r1enh", "r1sil", "r2enh", "r2sil", "r3enh", "r3sil")
    p_adj
}

fisher.region.boot <- function(exon, myHitsData) {
    n <- length(myHitsData)
    p_mat <- matrix(ncol = 6, nrow = n)
    out <- lapply(myHitsData, get_fisher, exon = exon)
    for (i in seq_len(n)) p_mat[i, ] <- out[[i]]
    p_adj <- apply(p_mat, 2, p.adjust, method = "BH")
    rownames(p_adj) <- names(myHitsData)
    colnames(p_adj) <- c("r1enh", "r1sil", "r2enh", "r2sil", "r3enh", "r3sil")
    cat("=")
    p_adj
}

# --- Read tetramer positional data --------------------------------------------

getTables <- function(protein, place, tets) {
    flist_path <- paste0(protein, place, "filelist.txt")
    if (!file.exists(flist_path)) return(list())

    filelist <- read.delim(flist_path, header = FALSE, stringsAsFactors = FALSE)[, 1]
    overlap  <- pmatch(tets, filelist)

    if (length(overlap) == 1 && is.na(overlap)) return(list())

    filelist <- filelist[na.omit(overlap)]
    sel_tets <- tets[!is.na(overlap)]
    rown <- as.character(ROWN)
    tab  <- matrix(0, nrow = length(ROWN), ncol = length(sel_tets),
                   dimnames = list(rown, sel_tets))
    res  <- list(enh = tab, sil = tab, cont = tab)

    for (i in seq_along(filelist)) {
        fpath <- paste0(protein, place, filelist[i])
        df <- tryCatch(read.delim(fpath, fill = TRUE, skip = 1, header = FALSE),
                       error = function(e) NULL)
        if (is.null(df) || nrow(df) <= 1) next
        colnames(df) <- c("pos", "cat", "no")
        tet <- sel_tets[i]
        for (cat_val in c(1, -1, 0)) {
            key <- switch(as.character(cat_val),
                          "1" = "enh", "-1" = "sil", "0" = "cont")
            chosen <- subset(df, cat == cat_val)
            rownames(chosen) <- as.character(chosen$pos)
            res[[key]][, tet] <- chosen[rown, "no"]
        }
    }
    for (i in seq_along(res)) res[[i]][is.na(res[[i]])] <- 0
    res
}

# --- Positional Fisher test ---------------------------------------------------

lmb.cluster.fisher.test <- function(d, tot_d, contr, tot_contr) {
    pvd <- matrix(0, nrow = nrow(d), ncol = ncol(d),
                  dimnames = list(rownames(d), colnames(d)))
    for (j in seq_len(ncol(pvd)))
        for (i in seq_len(nrow(pvd))) {
            tab <- matrix(c(d[i, j], contr[i, j],
                            tot_d - d[i, j], tot_contr - contr[i, j]), nrow = 2)
            pvd[i, j] <- fisher.test(tab, alternative = "greater")$p.value
        }
    pvd
}

# --- Tetramer sorting and clustering ------------------------------------------

getTrims <- function(mot) {
    sp <- unlist(strsplit(mot, ""))
    degen <- c("Y", "R", "W", "S")
    if (!any(sp %in% degen)) {
        return(c(paste0(sp[1:3], collapse = ""),
                 paste0(sp[2:4], collapse = "")))
    }
    trimers <- character(0)
    expand <- function(letter) {
        switch(letter,
               Y = c("C", "T"), R = c("A", "G"),
               W = c("T", "A"), S = c("C", "G"),
               letter)
    }
    for (x in expand(sp[1])) trimers <- c(trimers, paste0(x, sp[2], sp[3]))
    for (x in expand(sp[4])) trimers <- c(trimers, paste0(sp[2], sp[3], x))
    trimers
}

idAlignMot <- function(m, mall) {
    mot  <- getTrims(m)
    lall <- lapply(mall, getTrims)
    hits <- sapply(lall, function(x) length(na.omit(match(x, mot))))
    lens <- sapply(lall, length)
    which(hits / lens >= 0.5)
}

SortingAndType <- function(score, nn) {
    ord  <- character(0)
    type <- integer(0)
    ii   <- 1
    while (length(nn) > 0) {
        if (length(nn) > 1) {
            m   <- nn[1]
            nn  <- nn[-1]
            ids <- idAlignMot(m, nn)

            if (length(ids) == 0) {
                add <- m
            } else if (length(ids) == 1) {
                add <- c(m, nn[ids])
                nn  <- nn[-ids]
            } else {
                cc <- cp <- numeric(length(ids))
                for (k in seq_along(ids)) {
                    ct <- cor.test(score[, m], score[, nn[ids[k]]],
                                   method = "pearson")
                    cc[k] <- ct$estimate
                    cp[k] <- ct$p.value
                }
                xs  <- order(cc, cp, decreasing = TRUE)
                add <- c(m, nn[ids[xs]])
                nn  <- nn[-ids]
            }
            ord  <- c(ord, add)
            type <- c(type, rep(ii, length(add)))
        } else {
            ord  <- c(ord, nn)
            type <- c(type, ii)
            nn   <- NULL
        }
        ii <- ii + 1
    }
    list(ord, type)
}

# --- RNA splicing map visualisation -------------------------------------------

panel.rnaplot <- function(x, y, ..., ry = 29, rc = "red", intensity = 1) {
    cols <- myPalette(101)[findInterval(
        intensity,
        seq(0, max(intensity, na.rm = TRUE),
            length.out = 101))]
    for (i in seq_along(y))
        panel.polygon(c(x[i], x[i]), c(0, y[i]), border = cols[i], cex = 0.5)
    panel.rect(-30, 0, -10, ry, col = rc)
}

rnaScoreMap <- function(data, rcol, ylabels = c("0.5", "1"), main = "",
                        exon = 50, intron = 200, gap = 10, regions = 1000) {
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

    rownames(data) <- NULL
    pf <- cbind(data, reg = NA, val = NA, plot = NA)
    pf <- as.data.frame(pf)
    pf$reg <- floor((pf$l + intron) / regions)
    pf$val <- pf$l - (pf$reg * regions)
    for (r in c("1", "2", "3", "4"))
        pf$plot[pf$reg == r] <- pf$val[pf$reg == r] +
            offset$plus[offset$reg == r]

    pf_real <- pf[pf$l %in% pos_region, ]
    ylim <- c(0, max(pf_real$s))
    pf_real$rcol <- rcol[as.character(pf_real$g)]

    xyplot(s ~ plot | g, data = pf_real,
           intensity = pf_real$z,
           rcol = pf_real$rcol,
           hrect = ylim[2],
           aspect = "fill", type = "l",
           layout = c(1, length(levels(data$g))),
           ylim = ylim,
           xlim = c(-30, offset$end[4] + 10),
           subscripts = TRUE,
           strip = FALSE,
           strip.left = strip.custom(horizontal = TRUE,
                                     var.name = unique(pf_real$g),
                                     bg = GREY_BG),
           strip.left.lines = 0.5,
           panel = function(x, y, hrect, rcol, intensity, groups,
                            subscripts, ...) {
               panel.xblocks(
                   x = offset$start[1]:offset$end[4],
                   c(rep("0", 50), rep(NA, 201), rep(NA, 9),
                     rep(NA, 201), rep("3", 50), rep(NA, 9),
                     rep("4", 50), rep(NA, 201), rep(NA, 9),
                     rep(NA, 201), rep("7", 50)),
                   col = "gray96", border = FALSE, lwd = 0.5)
               panel.xyplot(x, y, ...)
               panel.xblocks(
                   x = offset$start[1]:offset$end[4],
                   c(rep(NA, 251), rep("0", 9), rep(NA, 251),
                     rep("0", 9), rep(NA, 251), rep("0", 9),
                     rep(NA, 251)),
                   col = "white", border = FALSE)
               panel.abline(
                   v = c(offset$plus[1], offset$plus[2] + 1,
                         offset$plus[3], offset$plus[4] + 1,
                         offset$start[1], offset$start[2],
                         offset$start[3], offset$start[4],
                         offset$end[1] + 1, offset$end[2] + 1,
                         offset$end[3] + 1, offset$end[4] + 1),
                   lty = "dotted", col = "black")
               panel.rnaplot(x, y, ry = hrect,
                             rc = rcol[subscripts],
                             intensity = intensity[subscripts], ...)
           },
           scales = list(
               x = list(rot = 90, labels = NULL, tck = c(0, 1)),
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
