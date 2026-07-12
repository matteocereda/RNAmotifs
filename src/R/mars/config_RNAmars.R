# Libraries ====

suppressPackageStartupMessages(library(circlize))
suppressMessages({
  library(ComplexHeatmap)
  library(dplyr)
  library(data.table)
  library(stringr)
  library(rlang)
  library(snow)
  library(viridis)
  library(lsa)
  library(optparse)
  library(BAMMtools)
})

options(stringsAsFactors = FALSE)


# Themes, palettes, handy variables ====

colors_map = list('cell_line' = c('HepG2'='#ffc300','K562'= '#7ac74f'),
                  'typeAS' = c('sil'='dodgerblue','enh'='firebrick'))


theme_big <- function(base_size = 12, base_family = "sans"){
  theme_classic()+
    theme(
      legend.text = element_text( size = 15), plot.title = element_text(size = 18), axis.title = element_text(size=18),
      axis.text = element_text(size=12),
      strip.text = element_text(size = 12)
    )
}

chr_levels = paste0('chr',c(1:22,'X','Y'))


# Region labels by event type ====

get_region_labels <- function(event_type = "SE") {
  if (event_type == "RI") {
    list(R1 = "upstream exon", R2 = "5' retained intron",
         R3 = "3' retained intron", R4 = "downstream exon",
         center = "retained intron", flanks = "flanking exons",
         title_suffix = " (Intron Retention)")
  } else {
    list(R1 = "upstream intron", R2 = "exon body",
         R3 = "downstream intron", R4 = "distal flank",
         center = "cassette exon", flanks = "flanking introns",
         title_suffix = "")
  }
}


# Functions ====


plotROI_Heatmap <- function(mat, contour_mat, regulation = c("enh", "sil"), name, palette = c("gray20", "chartreuse1","chartreuse3"), clustering = F) {
  typeAS = ifelse(colnames(mat) %like%'enh','enh','sil')
  mat = mat[, typeAS == regulation, drop = FALSE]
  contour_mat = contour_mat[, typeAS == regulation, drop = FALSE]
  typeAS = typeAS[typeAS == regulation]
  rng = range(mat, na.rm = TRUE)
  col_M = colorRamp2(seq(rng[1], rng[2], length.out = 3),palette)
  Heatmap(mat,
          name = "BS",
          col = col_M,
          cluster_columns = F,
          cluster_rows = clustering,
          column_split = typeAS,
          column_labels = gsub("_(enh|sil)$", "", colnames(mat)),
          top_annotation = HeatmapAnnotation(typeAS = typeAS, show_annotation_name = F, col = list(typeAS = c('sil' = 'dodgerblue','enh'='firebrick'))),
          row_title=NULL,
          split               = factor(rownames(mat), levels = rownames(mat)),
          gap                 = unit(0,"mm"),
          row_names_side   = "left",
          height              = unit(7, "mm")*nrow(mat),
          width               = unit(7, "mm")*ncol(mat),
          cell_fun            = function(j, i, x, y, width, height, fill) {
            grid.rect(x      = x,
                      y      = y,
                      width  = width,
                      height = height,
                      gp     = gpar(col = "grey", fill = "white", lwd = 1))
            grid.rect(x      = x,
                      y      = y,
                      width  = width*mat[i,j]*0.9,
                      height = height*mat[i,j]*0.9,
                      gp = gpar(fill = col_M(mat[i,j]),  col = contour_mat[i,j], lwd=1.5 )
            )
          },
          column_title = name
  )
}


getdim = function(H) {

  ht_width  = sum(component_width(H)) + unit(4, "mm")
  ht_width  = convertWidth(ht_width, "inch", valueOnly = TRUE)

  ht_height = sum(component_height(H)) + unit(4, "mm")
  ht_height = convertHeight(ht_height, "inch", valueOnly = TRUE)

  return(c(width = ht_width, height = ht_height))
}



plot_association_heatmap <- function(MAT, PARAMS, PEAK, sub_res_df, typeAS, input_name, ranking_method = 'rowmean', sub_deseq, BIND_MAT = NULL, top_mrms = 0, top_rbps = 0) {
  # top_mrms / top_rbps: if > 0, keep only the top-N MRMs (columns, by enrichment
  # score) and top-N RBPs (rows, by association-score rank). 0 = show all.
  # - MAT is a list of association score matrices.
  #   MAT[[params]] is the matrix where in each row there is the protein and in each column the tetramer
  # - PEAK: normalized binding profile from trained RBPs
  # - PARAMS: dataframe where each row is a protein with the corresponding optimal parameters [,c('pr', 'params', 'hw','ew')]
  # - sub_res_df: dataframe with RNAmotifs enrichment pvalues. [,c('input_pr', 'tetramer', 'r1enh_pFis', 'r2enh_pFis', 'r3enh_pFis', 'r1sil_pFis', 'r2sil_pFis', 'r3sil_pFis', 'key' )]
  #         comes from the function importSummaryFile() --> "HepG2", "HepG2", rbps_h, subfolder)
  # - typeAS = "sil" or "enh"
  # - ranking_method = default "rowmean"
  # - sub_deseq = deseq results from same RNAseq data that produced the list of input alternative exons.
  #               It must be in the DESeq2 format and contain the column `deg` (TRUE/FALSE), otherwise these thresholds will be used: abs(log2FC) >= 0.1 and padj < 0.1

  all_tets = unique(sub_res_df$tetramer)
  length(all_tets)
  params_to_use = unique(PARAMS$params)


  final_mat = matrix(NA, nrow = length(rbps), ncol = length(all_tets), dimnames = list(rbps, all_tets))
  tet_score_in_run_mat = matrix(NA, nrow = length(rbps), ncol = length(all_tets), dimnames = list(rbps, all_tets))

  # Fill final_mat only if tetramer is enriched in optimal params for row protein
  for (p in params_to_use) {
    mat = MAT[[p]]
    if (!is.matrix(mat) ) {next}
    mat = mat[, colnames(mat) %in% all_tets, drop = F]
    tetramers = colnames(mat)
    proteins=subset(PARAMS, params == p)$pr
    proteins = proteins[proteins %in% rownames(mat)]   # guard: skip RBPs with no scores in this run (e.g. degenerate enh set of a pure repressor)
    if (length(proteins) == 0) next
    final_mat[proteins,tetramers] <- mat[proteins,,drop=FALSE]

    # tetramers enrichment in the specific run
    sub_res_df_RUN = subset(sub_res_df, key==p)
    if(nrow(sub_res_df_RUN)>=1){
      tet_score_in_run = tapply(1:nrow(sub_res_df_RUN), sub_res_df_RUN$tetramer, function(x){sum(-2*log(sub_res_df_RUN[x, paste0("r",c(1,2,3),typeAS,"_pFis")] ) ) })
      tet_score_in_run_mat[proteins, names(tet_score_in_run)] = tet_score_in_run
    }
  }

  # Column annotation: tetramer scores region-wise
  if(nrow(sub_res_df)>1){
    tetr_pval_r1     = tapply(1:nrow(sub_res_df),sub_res_df$tetramer, function(x){sum(-2*log(sub_res_df[x, paste0("r1",typeAS,"_pFis")] )) })
    tetr_pval_r2     = tapply(1:nrow(sub_res_df),sub_res_df$tetramer, function(x){sum(-2*log(sub_res_df[x, paste0("r2",typeAS,"_pFis")] )) })
    tetr_pval_r3     = tapply(1:nrow(sub_res_df),sub_res_df$tetramer, function(x){sum(-2*log(sub_res_df[x, paste0("r3",typeAS,"_pFis")] )) })
  } else {
    tetr_pval_r1     = sum(-2*log(sub_res_df[, paste0("r1",typeAS,"_pFis")]))
    tetr_pval_r2     = sum(-2*log(sub_res_df[, paste0("r2",typeAS,"_pFis")]))
    tetr_pval_r3     = sum(-2*log(sub_res_df[, paste0("r3",typeAS,"_pFis")]))
    names(tetr_pval_r1) = unique(sub_res_df$tetramer)
    names(tetr_pval_r2) = unique(sub_res_df$tetramer)
    names(tetr_pval_r3) = unique(sub_res_df$tetramer)
  }

  # Cap Inf values (from -2*log(0)) to avoid colorRamp2 failures
  cap_inf <- function(x) { x[is.infinite(x)] <- max(x[is.finite(x)], na.rm=TRUE) * 1.1; x }
  tetr_pval_r1 <- cap_inf(tetr_pval_r1)
  tetr_pval_r2 <- cap_inf(tetr_pval_r2)
  tetr_pval_r3 <- cap_inf(tetr_pval_r3)

  all_tet_scores = colSums(rbind(tetr_pval_r1, tetr_pval_r2, tetr_pval_r3))
  all_tet_scores = all_tet_scores[order(all_tet_scores, decreasing = T)]

  ## GENE EXPRESSION ANNOTATION
  padj_thres = 0.1
  fc_thresh = 0.1

  # is there a user-provided DESeq file?
  deseq = !(all(sub_deseq$log2FoldChange == 0, na.rm = TRUE) && all(sub_deseq$padj == 1 | is.na(sub_deseq$padj)))

  if (deseq && 'deg' %in% colnames(sub_deseq)){
    fc = sub_deseq[,c('gene_name', 'log2FoldChange','padj', 'deg')]
    rownames(fc) = fc$gene_name
  } else if (deseq && !('deg' %in% colnames(sub_deseq))){
    message('[*] DESeq2 table does not contain `deg` column, abs(log2FC) >= ', fc_thresh, ' and padj <= ', padj_thres,' will be annotated as deg')
    fc = sub_deseq[,c('gene_name', 'log2FoldChange','padj')]
    fc$deg = ifelse(abs(fc$log2FoldChange) >= fc_thresh & fc$padj <= padj_thres, TRUE, FALSE)
    rownames(fc) = fc$gene_name
  } else {
    message('[*] DESeq2 table not provided!')
    fc = sub_deseq[,c('gene_name', 'log2FoldChange','padj')]
    rownames(fc) = fc$gene_name
  }

  if (deseq){
  sign_fc = sign(fc$log2FoldChange)

  shape_fc = case_when(
    sign_fc == 1 & fc$deg == TRUE   ~ 24,
    sign_fc == 1 & fc$deg == FALSE  ~ 2,
    sign_fc == -1 & fc$deg == TRUE  ~ 25,
    sign_fc == -1 & fc$deg == FALSE ~ 6,
    sign_fc == 0                    ~ 1
    )

  shape_fc[is.na(shape_fc) & sign_fc==1]<- 2
  shape_fc[is.na(shape_fc) & sign_fc== -1]<- 6

  # Direction colour: up-regulated (log2FC > 0) = red, down-regulated (< 0) = blue,
  # unchanged (0) = grey. DEG-vs-not is still conveyed by the shape (filled 24/25
  # vs open 2/6 triangles); the colour now encodes the sign of the fold change.
  color_fc = case_when(
    sign_fc == 1  ~ 'red',
    sign_fc == -1 ~ 'blue',
    TRUE          ~ 'grey'
  )

  names(shape_fc) = fc$gene_name
  names(color_fc) = fc$gene_name

  padj_score = ifelse(fc$padj<=padj_thres & !is.na(fc$padj), -log10(fc$padj), 1)
  names(padj_score) = fc$gene_name
  padj_score[!is.finite(padj_score)] = max(padj_score[is.finite(padj_score)])+1
  }

  ## Color for the top annotation cells
  if (typeAS=='sil') {
    typeAS_2 = 'Silenced'
    tetr_regions_col = c("white","darkblue")
  } else {
    typeAS_2 = 'Enhanced'
    tetr_regions_col = c("white","darkred")
  }

  # ── Build the WHOLE heatmap first; top-N is a pure DISPLAY subset (a zoom). ──
  # Ranking, colour scale and bubble normalisation are all computed on the full
  # matrix, so a top-N view keeps the SAME row order, mean(AS), colours and bubble
  # sizes as the full figure — it does not recompute anything, only crops.

  # Column order: all MRMs by enrichment score (all_tet_scores already sorted desc)
  colorder_full = names(all_tet_scores)
  final_mat = final_mat[, colorder_full, drop = F]

  # Row ranking over the WHOLE matrix (all MRMs) — independent of top-N
  if (ranking_method =='rowmean') {
    row_rank = sort(rowMeans(final_mat, na.rm=T), decreasing = T)
  } else if (ranking_method =='rowmax') {
    rmax = apply(final_mat, 1, function(x){ v = x[!is.na(x)]; if (length(v)==0) NA_real_ else max(v) })
    row_rank = sort(rmax, decreasing = T, na.last = NA)
  } else {
    stop(sprintf("Unknown ranking_method '%s' (expected 'rowmean' or 'rowmax')", ranking_method))
  }
  row_rank = row_rank[row_rank != 0]
  row_rank = row_rank[!is.na(row_rank)]

  # Colour scale + bubble normalisation fixed on the WHOLE heatmap (zoom-stable)
  mat_min_full = min(final_mat, na.rm = TRUE)
  mat_max_full = max(final_mat, na.rm = TRUE)

  # DISPLAY SUBSET: top-N RBPs (rows) and top-N MRMs (columns); 0 = all
  roworder = names(row_rank)
  if (top_rbps > 0 && length(roworder) > top_rbps) roworder = roworder[seq_len(top_rbps)]
  row_rank = row_rank[roworder]
  colorder = colorder_full
  if (top_mrms > 0 && length(colorder) > top_mrms) colorder = colorder[seq_len(top_mrms)]
  final_mat = final_mat[roworder, colorder, drop = F]
  motif_score = row_rank

  all_tet_scores = as.vector(all_tet_scores[colnames(final_mat)])

  # To count as weight the tetramer score
  tet_score_mat = matrix(rep(all_tet_scores, nrow(final_mat)), byrow = T, nrow =  nrow(final_mat))
  colnames(tet_score_mat)<- all_tet_scores
  rownames(tet_score_mat)<- roworder

  color_barplot_tet_score = rep('darkcyan', length(all_tet_scores))

  # Maximum peak in region (subset to displayed rows)
  binding = PEAK[roworder,, drop=F]
  roi = t(apply(PEAK[roworder,, drop=F], 1, function(x) { ifelse(x==max(x),'black','white') }))

  # reorder fc
  fc  = fc[roworder, ]
  if (deseq){
    shape_fc = shape_fc[roworder]
    color_fc = color_fc[roworder]
    padj_score = padj_score[roworder]
  }

  left_annot <- plotROI_Heatmap(binding, roi, regulation = typeAS, name = NULL)

  # color scale for main matrix — fixed on the whole heatmap (zoom-stable)
  if (mat_max_full > mat_min_full) {
    cols = colorRamp2(seq(mat_min_full, mat_max_full, length.out = 101), viridis(101))
  } else { cols = viridis(1)}

  motif_score_mat = as.matrix(motif_score)
  rowmean_heatmap <- Heatmap(motif_score_mat,
                            col = cols,
                            cluster_rows = F, cluster_columns = F,
                            height              = unit(8, "mm")*nrow(motif_score_mat),
                            width               = unit(8, "mm")*ncol(motif_score_mat))

  anno_pval  = anno_barplot(all_tet_scores,
                            border=F,
                            bar_width = 0.9,
                            gp = gpar(fill = color_barplot_tet_score,col=NA))

  max_r_pval = max(c(tetr_pval_r1, tetr_pval_r2, tetr_pval_r3), na.rm = TRUE)

  pval_col_fun = colorRamp2(c(0, max_r_pval), tetr_regions_col)

  top_anno   = columnAnnotation("CCRE" = anno_pval,
                                annotation_name_side = "left", annotation_name_rot = 0,
                                'R1' = anno_simple(tetr_pval_r1[colnames(final_mat)], which = "row", gp = gpar(col="darkgrey"), col = pval_col_fun),
                                'R2' = anno_simple(tetr_pval_r2[colnames(final_mat)], which = "row", gp = gpar(col="darkgrey"), col = pval_col_fun),
                                'R3' = anno_simple(tetr_pval_r3[colnames(final_mat)], which = "row", gp = gpar(col="darkgrey"), col = pval_col_fun)
  )

  pval_legend = Legend(title = paste0("CRE"), col_fun = pval_col_fun, at = c(0, round(max_r_pval/2,1), round(max_r_pval,1)))

  if (deseq) {
  mat_annot_fc = matrix(abs(fc$log2FoldChange))
  rownames(mat_annot_fc) = rownames(fc)
  far_left_anno = rowAnnotation(`|log2FC|` = anno_points(mat_annot_fc,
                                           pch = shape_fc,
                                           gp = gpar(col = color_fc, fill = color_fc),
                                           size = unit(4, "mm"),
                                           width = unit(20, "mm"),
                                           border     = F))
  } else {
  far_left_anno = NULL
  }

  h <- Heatmap(final_mat,
               col                 = cols,
               name                = 'AS',
               cluster_rows        = F,
               cluster_columns     = F,
               height              = unit(8, "mm")*nrow(final_mat),
               width               = unit(8, "mm")*ncol(final_mat),
               top_annotation      = top_anno,
               split               = factor(roworder, levels = roworder),
               gap                 = unit(0,"mm"),
               row_title           = NULL,
               show_row_names      = F,
               right_annotation    = rowAnnotation("splicingMaps"       = anno_empty(which="row",width=unit(8,"cm")),
                                                   "empty6"             = anno_empty(which="row",width=unit(10,"mm"), border=F)),
               left_annotation     = rowAnnotation(`mean(AS)` = anno_barplot( row_rank,  border     = F )),
               column_title        = "",
               column_names_side   = "top",
               row_names_side      = "left",

               cell_fun            = function(j, i, x, y, width, height, fill) {
                 w                 = 2*final_mat[i,j]/mat_max_full

                 grid.rect(x      = x,
                           y      = y,
                           width  = width,
                           height = height,
                           gp     = gpar(col = "white", fill = "white", lwd = 1))
                 grid.segments(x0 = x,
                               y0 = y-height*0.5,
                               x1 = x,
                               y1 = y+height*0.5,
                               gp = gpar(col="ivory3"))

                 grid.segments(x0 = x-width*0.5,
                               y0 = y,
                               x1 = x+width*0.5,
                               y1 = y,
                               gp = gpar(col="ivory3"))
                 grid.circle(
                   x  = x,
                   y  = y,
                   r  = w*unit(1.5,"mm"),
                   gp = gpar(fill = if (is.function(cols)) cols(final_mat[i,j]) else cols, col = 'grey')
                 )
               }
  )

  row_names_annot <- rowAnnotation(rn = anno_text(rownames(fc)))
  if (!is.null(far_left_anno)) {
    hm = row_names_annot + far_left_anno + left_annot + h
  } else {
    hm = row_names_annot + left_annot + h
  }
  return(list(hm, final_mat, row_rank, all_tet_scores, pval_legend))
}



add_splicing_maps_to_final_ht <- function(sub_splicingMaps, mat, params, inExon = 30, inIntron = 300) {
  gap    = 15
  max_r1 = 0

  # find maximum height of the splicing maps
  for(sl in 1:nrow(mat)){
    row_rbp      = rownames(mat)[sl]
    tmp_df   = subset(sub_splicingMaps, hw==params$hw[match(row_rbp, params$pr)] & ew==params$ew[match(row_rbp, params$pr)])[,c("height_enh","height_sil","height","col")]
    if(nrow(tmp_df)==0){next}

    if(typeAS=="both"){
      max_r1 = c(max_r1,max(c(0,as.numeric(tmp_df$height)),na.rm=T))
    } else if(typeAS=="enh"){
      max_r1 = c(max_r1,max(c(0,as.numeric(tmp_df$height_enh)),na.rm=T))
    } else if(typeAS=="sil"){
      max_r1 = c(max_r1,max(c(0,as.numeric(tmp_df$height_sil)),na.rm=T))
    }
  }

  max_r1 = max(max_r1)

  # Reference splicing-map length (rows = 4 regions x positions). Used to draw an
  # empty exon skeleton for any row whose parameter combo has no map (e.g. a combo
  # with no enriched tetramers), so EVERY row shows a map rather than a blank box.
  ref_len = 0
  for(sl in 1:nrow(mat)){
    rr = rownames(mat)[sl]
    td = subset(sub_splicingMaps, hw==params$hw[match(rr, params$pr)] & ew==params$ew[match(rr, params$pr)])
    if(nrow(td) > 0){ ref_len = nrow(td); break }
  }

  for(sl in 1:nrow(mat)){
    row_rbp      = rownames(mat)[sl]
    tmp_df   = subset(sub_splicingMaps, hw==params$hw[match(row_rbp, params$pr)] & ew==params$ew[match(row_rbp, params$pr)])[,c("height_enh","height_sil","height","col")]

    # No map for this combo: inject a zero map of the reference length so the
    # empty-skeleton branch (max(r1)==0) draws the exon diagram instead of skipping.
    if(nrow(tmp_df)==0){
      if(ref_len==0){next}
      tmp_df = data.frame(height_enh=rep(0,ref_len), height_sil=rep(0,ref_len),
                          height=rep(0,ref_len), col=rep(0,ref_len))
    }

    if(typeAS=="both"){
      r1   = as.numeric(tmp_df$height)
      r2   = as.numeric(tmp_df$col)
    } else if(typeAS=="enh"){
      r1   = as.numeric(tmp_df$height_enh)
      r2   = rep(1,length(r1))
    } else if(typeAS=="sil"){
      r1   = as.numeric(tmp_df$height_sil)
      r2   = rep(0,length(r1))
    }

    IE = length(tmp_df$height)/4 -1
    if(max(r1)==0){

      decorate_annotation("splicingMaps", slice = sl,{

        pushViewport(viewport(xscale = c(0.5, IE + 3*(IE + gap) + 0.5), yscale = c(0, 50*1.1)))

        # rectangular stripes representing exons
        grid.rect(x             = inExon/2 + c(1, inIntron*2+inExon+1 + gap, IE*2+2 + 2*gap, inIntron + 3*(IE + gap)),
                  y             = 0,
                  width         = inExon,
                  height        = 50*1.05,
                  just          = "bottom",
                  gp            = gpar(fill = "snow2", col = NA),
                  default.units = "native")

        # dashed lines at the splice sites and at the gap -+ 300 nt from the exons (between flanking and alternative exons)
        grid.rect(x             = c(IE+1 + 0*gap + gap/2,
                                    IE*2+2 + 1*gap + gap/2,
                                    IE*3+3 + 2*gap + gap/2),
                  width         = gap,
                  y             = 0,
                  height        = 50*1.1,
                  just          = "bottom",
                  default.units = "native",
                  gp            = gpar(col = "black", lwd = 0.5, lty = 2, fill = NA))

        popViewport()
      })

      next
    }


    col_fun = colorRamp(c("blue4","dodgerblue4","blue","yellow","red","firebrick3","red4"))

    decorate_annotation("splicingMaps", slice = sl,{

      col_map = rgb(col_fun(r2),maxColorValue = 255)

      pushViewport(viewport(xscale = c(0.5, IE + 3*(IE + gap) + 0.5), yscale = c(0, max_r1*1.1)))

      # plot the frame with maximum height equal to max_r1
      grid.rect(x             = inExon/2 + c(1, inIntron*2 + inExon+1 + gap, IE*2+2 + 2*gap, inIntron + 3*(IE + gap)),
                y             = 0,
                width         = inExon,
                height        = max_r1*1.05,
                just          = "bottom",
                gp            = gpar(fill = "snow2", col = NA),
                default.units = "native")

      invisible(lapply(0:3,function(n){

        grid.rect(x             = (IE + gap)*n + 1:IE,
                  y             = 0,
                  width         = 1,
                  height        = r1[IE*n+(1:IE)],
                  just          = "bottom",
                  gp            = gpar(fill = col_map[IE*n+(1:IE)], col = col_map[IE*n+(1:IE)]),
                  default.units = "native")

        grid.rect(x             = c(IE+1 + 0*gap + gap/2,
                                    IE*2+2 + 1*gap + gap/2,
                                    IE*3+3 + 2*gap + gap/2),
                  width         = gap,
                  y             = 0,
                  height        = max_r1*1.1,
                  just          = "bottom",
                  default.units = "native",
                  gp            = gpar(col = "black", lwd = 0.5,lty=2,fill = NA))

      }))

      grid.yaxis(at =  max_r1*1.1, main = F, label = round(max_r1*1.1))

      popViewport()

    })
  }
}




adj_reg_mat = function(regions, enr_tetramers){

  if(ncol(regions) == 2){

    regions           = cbind(enr_tetramers,t(as.matrix(regions$x)))
    dimnames(regions) = list(enr_tetramers,c("X",paste0("R",1:3)))
    regions           = as.data.frame(regions)

    return(regions)

  } else {
    return(regions)
  }
}


get_exons_with_tetramer_region1 = function(tet, RESULTS, SPLICING_FILE){
  # this function keeps only input file elements that have a hit of tetramer 'tet' in the region 1
  # RESULTS = RNAmotifs result of specific protein
  # SPLICING_FILE = RNAmotifs input file of the same protein

  sp      = read.table(SPLICING_FILE, sep=";", h=F)
  f       = list.files(RESULTS, tet, full.names = T, recursive = T)
  rc      = read.delim(f[grep("region_count.tsv",f)])
  bed     = read.table(f[grep("bed",f)], sep="\t", skip=1)
  sel     = subset(rc, type!="0" & (hits_region1!=0))
  sp.sel  = subset(sp, V1%in%sel$myRID)
  sel     = cbind(sp.sel[,3:8],sel[match(sp.sel$V1,sel$myRID),])

  if(length(sel$V3)==0){
    return(data.frame(key=character()))
  }

  sel$key = with(sel, paste0(V3, ":", V6, "-", V7))

  return(sel)
}

get_exons_with_tetramer_region2 = function(tet, RESULTS, SPLICING_FILE){
  # this function keeps only input file elements that have a hit of tetramer 'tet' in the region 2
  # RESULTS = RNAmotifs result of specific protein
  # SPLICING_FILE = RNAmotifs input file of the same protein

  sp      = read.table(SPLICING_FILE, sep=";", h=F)
  f       = list.files(RESULTS, tet, full.names = T, recursive = T)
  rc      = read.delim(f[grep("region_count.tsv",f)])
  bed     = read.table(f[grep("bed",f)], sep="\t", skip=1)
  sel     = subset(rc, type!="0" & (hits_region2!=0))
  sp.sel  = subset(sp, V1%in%sel$myRID)
  sel     = cbind(sp.sel[,3:8],sel[match(sp.sel$V1,sel$myRID),])

  if(length(sel$V3)==0){
    return(data.frame(key=character()))
  }

  sel$key = with(sel, paste0(V3, ":", V6, "-", V7))

  return(sel)
}

get_exons_with_tetramer_region3 = function(tet, RESULTS,SPLICING_FILE){
  # this function keeps only input file elements that have a hit of tetramer 'tet' in the region 3
  # RESULTS = RNAmotifs result of specific protein
  # SPLICING_FILE = RNAmotifs input file of the same protein
  sp      = read.table(SPLICING_FILE, sep=";", h=F)
  f       = list.files(RESULTS, tet, full.names = T, recursive = T)
  rc      = read.delim(f[grep("region_count.tsv",f)])
  bed     = read.table(f[grep("bed",f)], sep="\t", skip=1)
  sel     = subset(rc, type!="0" & (hits_region3!=0))
  sp.sel  = subset(sp, V1%in%sel$myRID)
  sel     = cbind(sp.sel[,3:8],sel[match(sp.sel$V1,sel$myRID),])

  if(length(sel$V3)==0){
    return(data.frame(key = character()))
  }

  sel$key = with(sel, paste0(V3, ":", V6, "-", V7))

  return(sel)
}



lmb.cluster.fisher.test = function(d, tot.d, contr, tot.contr ){

  pvd           = matrix(0,nrow = nrow(d), ncol = ncol(d))
  colnames(pvd) = colnames(d)

  for(j in 1:ncol(pvd))
    for(i in 1:nrow(pvd)){
      test     = matrix( data=c( d[i,j], contr[i,j], tot.d - d[i,j] , tot.contr - contr[i,j] ), nrow=2 )
      pvd[i,j] = fisher.test(test, alternative = "greater")$p.value
    }

  pvd
}



SEeCLIPpeaks = function(datExpr,peaks,inIntron, inExon, n_cores){
  # datExpr: matrix with rows like 1;10046;chr6;+;43746655;43749692;43749824;43752277;1
  # peaks: matrix with rows like chr1	17496	17498

  m1 = matrix(rep(-inExon:inIntron,nrow(datExpr)),
              nrow  = nrow(datExpr),
              ncol  = length(-inExon:inIntron),
              byrow = T)

  m2 = matrix(rep(-inIntron:inExon,nrow(datExpr)),
              nrow  = nrow(datExpr),
              ncol  = length(-inIntron:inExon),
              byrow = T)

  pos = cbind(m1 + datExpr$V5,
              m2 + datExpr$V6,
              m1 + datExpr$V7,
              m2 + datExpr$V8)

  clus = snow::makeCluster(n_cores, type = "SOCK")
  snow::clusterExport(clus,list = c("datExpr","peaks","inIntron","inExon"))
  values_t = snow::parSapply(clus, 1:nrow(datExpr), function(x){

    m1 = matrix(rep(-inExon:inIntron,nrow(datExpr)),
                nrow  = nrow(datExpr),
                ncol  = length(-inExon:inIntron),
                byrow = T)

    m2 = matrix(rep(-inIntron:inExon,nrow(datExpr)),
                nrow  = nrow(datExpr),
                ncol  = length(-inIntron:inExon),
                byrow = T)

    pos = cbind(m1 + datExpr$V5,
                m2 + datExpr$V6,
                m1 + datExpr$V7,
                m2 + datExpr$V8)

    tmp0 = subset(peaks, V1==datExpr$V3[x])
    tmp  = subset(tmp0, V2%in%pos[x,] | V3%in%pos[x,])

    if(nrow(tmp)==0){
      return(pos[x,]*0)
    }

    tmp         = do.call(rbind, lapply(1:nrow(tmp),function(y){cbind(seq(as.numeric(tmp$V2[y]),as.numeric(tmp$V3[y])),1)}))
    a           = tmp[match(pos[x,],tmp[,1]),2]
    a[is.na(a)] = 0

    if(datExpr$V4[x]=="-"){a = a[length(a):1]}

    return(a)
  })

  stopCluster(clus)

  res = list(pos = pos, values = t(values_t))

  return(res)
}
