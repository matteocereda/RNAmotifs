## ___________________________
## Script name: generate_heatmap.R
## Purpose of script: Uses the Signal Recovery Rate (SRR) and Cosine Similarity (CS) score matrices derived from the previous step (compute_association_scores) to aggregate them and produce the final heatmap.
## ___________________________

library(optparse)

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="RNAmotifs input txt file", metavar="character"),
  make_option(c("-n", "--name"), type="character", default=NULL,
              help="RNAmotifs and RNAMaRs run name", metavar="character"),
  make_option(c("-d", "--RNAmotifs_res"), type="character", default=NULL,
              help="RNAmotifs results folder", metavar="character"),
  make_option(c("-r", "--repository"), type="character", default=NULL,
              help="RNAmars repository", metavar="character"),
  make_option(c("-c", "--cell_line"), type="character", default=NULL,
              help="Cell line to compare input exons to (HepG2 or K562)", metavar="character"),
  make_option(c("-e", "--deseq_file"), type="character", default=NULL,
              help="DESEQ2 expression file", metavar="character"),
  make_option(c("-p", "--cores"), type="character", default=NULL,
              help="Number of cpus", metavar="numerics"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Output directory", metavar="character"),
  make_option(c("-s", "--scripts_dir"), type="character", default=NULL,
              help="Path to mars R scripts directory", metavar="character"),
  make_option(c("--in_exon"), type="integer", default=30,
              help="Region extent into exons in bp (default: 30)", metavar="integer"),
  make_option(c("--in_intron"), type="integer", default=300,
              help="Region extent into introns in bp (default: 300)", metavar="integer"),
  make_option(c("--sort_by"), type="character", default="mean",
              help="Sort RBPs in the heatmap by association score: 'mean' (default) or 'max'", metavar="character"),
  make_option(c("--top_mrms"), type="integer", default=0,
              help="Show only the top-N enriched MRMs (columns, by enrichment score). 0 = all (default)", metavar="integer"),
  make_option(c("--top_rbps"), type="integer", default=0,
              help="Show only the top-N RBPs (rows, by association-score rank). 0 = all (default)", metavar="integer")
);


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

READY_RNAMOTIFS_INPUT = opt$input
name                  = opt$name
DIR_RESULTS_RNAMOTIFS = opt$RNAmotifs_res
repo                  = opt$repository
cell_line_target      = opt$cell_line
deseq_file            = opt$deseq_file
cpus                  = as.numeric(opt$cores)
DIR_OUTPUT            = opt$output
scripts_dir           = opt$scripts_dir
inExon_param          = opt$in_exon
inIntron_param        = opt$in_intron
ranking_method        = if (tolower(opt$sort_by) == "max") "rowmax" else "rowmean"
top_mrms_param        = opt$top_mrms
top_rbps_param        = opt$top_rbps


# deseq_file = '' # either '' or the tsv with these columns:
## baseMean	log2FoldChange	lfcSE	stat	pvalue	padj

if (is.na(deseq_file) || deseq_file == "NA") {
  deseq_file <- ""
}

# Set up =========================================================
source(file.path(scripts_dir, 'config_RNAmars.R'))


# TRAINING RBPS DEFINITION ================
if (cell_line_target == "HepG2"){
  rbps = c("HNRNPC", "HNRNPK", "HNRNPU", "NCBP2", "PRPF8", "PTBP1", "QKI", "RBFOX2", "RBM22", "SF3A3", "SF3B4", "SRSF1", "U2AF1", "U2AF2", "UCHL5")
} else if (cell_line_target == "K562"){
  rbps = c("AGGF1", "EFTUD2", "FXR1", "HNRNPU", "PRPF8", "PTBP1", "PUS1", "RBM15", "SF3B4", "SRSF1", "TARDBP", "U2AF1", "U2AF2")
}


## Optimal set of params
OPT    = read.csv(file.path(repo, "Tables/RNAmotifs_optimal_parameters.csv"), stringsAsFactors = FALSE)
PARAMS = subset(OPT, true_rbp %in% rbps)
colnames(PARAMS)[colnames(PARAMS)=='true_rbp']<- 'pr'
sel_par_sets <- unique(PARAMS[, c("hw", "ew")])
sel_par_sets$hw <- as.numeric(sel_par_sets$hw)
sel_par_sets$ew <- as.numeric(sel_par_sets$ew)

## Peak - try TSV first, fall back to RDS
peak_tsv <- file.path(repo, "Rdata", paste0(cell_line_target, "_binding_profile_PEAK_normalized.tsv"))
peak_rds <- file.path(repo, "Rdata", paste0(cell_line_target, "_binding_profile_PEAK_normalized.rds"))
if (file.exists(peak_tsv)) {
  PEAK = as.matrix(read.delim(peak_tsv, row.names=1, check.names=FALSE))
} else {
  PEAK = readRDS(peak_rds)
}


# Splicing maps  ================

message(noquote("\n[*] Creating splicing maps using RBP-specific optimized parameters..."))
all_splicingMaps = do.call(rbind,lapply(unique(OPT$params), function(y){
  params = unlist(str_split(y, '_'))[c(2,4)]

  path_folder   = paste0(DIR_RESULTS_RNAMOTIFS, name, '_',y)

  if(length(grep("^All_group_enh-",list.files(path=path_folder),value=T))==0){return()}

  # Read the first matching group-map CSV; return NULL if absent.
  read_group <- function(prefix){
    fs <- grep(paste0("^",prefix), list.files(path=path_folder), value=T)
    if(length(fs)==0) return(NULL)
    read.csv(paste0(path_folder,"/",fs[1]))[,1]
  }
  e_spl_map  = read_group("All_group_enh-")
  ES_spl_map = read_group("All_group_both-")
  all_sil    = read_group("All_group_sil-")

  # Direction-specific height: prefer the direction-CLUSTER map (Sil_group_sil- /
  # Enh_group_enh-), but many subjects have tetramers clustering in only ONE
  # direction, leaving the other empty. When the cluster-specific map is absent or
  # all-zero, fall back to the ALL-tetramer map over that exon set (All_group_sil- /
  # All_group_enh-) so every row shows its RNA splicing map rather than a blank line.
  sil_cluster = read_group("Sil_group_sil-")
  if(is.null(sil_cluster) || sum(sil_cluster, na.rm=TRUE)==0){
    SIL_spl_map = if(!is.null(all_sil)) all_sil else ES_spl_map*0
  } else {
    SIL_spl_map = sil_cluster
  }

  enh_cluster = read_group("Enh_group_enh-")
  if(is.null(enh_cluster) || sum(enh_cluster, na.rm=TRUE)==0){
    ENH_spl_map = if(!is.null(e_spl_map)) e_spl_map else ES_spl_map*0
  } else {
    ENH_spl_map = enh_cluster
  }

  tmp_col                 = e_spl_map/ES_spl_map
  tmp_col[is.na(tmp_col)] = 0

  return(suppressWarnings(cbind(pr         = name,
                                hw         = params[1],
                                ew         = params[2],
                                rbp        = name,
                                pos        = 1:length(SIL_spl_map),
                                height_enh = ENH_spl_map,
                                height_sil = SIL_spl_map,
                                height     = ES_spl_map,
                                col        = tmp_col)))
}))


all_splicingMaps            = as.data.frame(all_splicingMaps)
all_splicingMaps$pos        = as.numeric(all_splicingMaps$pos)
all_splicingMaps$height_enh = as.numeric(all_splicingMaps$height_enh)
all_splicingMaps$height_sil = as.numeric(all_splicingMaps$height_sil)

message(noquote("[*] RNAmotifs splicing maps created and ready for final heatmap."))

# DEseq2 ===============================
message(noquote("\n[*] Importing DESeq2 differential genes..."))
if (file.exists(deseq_file) & deseq_file %like% '.tsv') {
  deseq     = read.delim(deseq_file)
  sub_deseq = subset(deseq, gene_name %in% rbps)
} else if (file.exists(deseq_file) & deseq_file %like% '.rds') {
  deseq     = as.data.frame(readRDS(deseq_file))
  sub_deseq = subset(deseq, gene_name %in% rbps)
} else {
  # baseMean	log2FoldChange	lfcSE	stat	pvalue	padj
  sub_deseq = data.frame(log2FoldChange = rep(0, length(rbps)), padj =  rep(1, length(rbps)), gene_name =rbps)
}



# Import scores ===================
# Read TSV score files from C++ output: SCORE{1,2}_{enh,sil}_hw_X_ew_Y_{name}.tsv
# Each file is a matrix with rows=RBPs, cols=tetramers
read_score_tsvs <- function(dir, prefix, params_list, name) {
  result <- list()
  for (i in 1:nrow(params_list)) {
    key <- paste0("hw_", params_list$hw[i], "_ew_", params_list$ew[i])
    f <- file.path(dir, paste0(prefix, "_", key, "_", name, ".tsv"))
    if (file.exists(f)) {
      result[[key]] <- as.matrix(read.delim(f, row.names=1, check.names=FALSE))
    }
  }
  result
}

DIR_DIAG = file.path(DIR_OUTPUT, "diagnostics")
srr_sil         = read_score_tsvs(DIR_DIAG, "SCORE1_sil", sel_par_sets, name)
profile_sim_sil = read_score_tsvs(DIR_DIAG, "SCORE2_sil", sel_par_sets, name)
srr_enh         = read_score_tsvs(DIR_DIAG, "SCORE1_enh", sel_par_sets, name)
profile_sim_enh = read_score_tsvs(DIR_DIAG, "SCORE2_enh", sel_par_sets, name)

ENH = lapply(names(profile_sim_enh), function(params) {
  srr_mat = srr_enh[[params]]
  prof_mat = profile_sim_enh[[params]]
  res = prof_mat*srr_mat
  return(res)
})
names(ENH) = names(profile_sim_enh)
ENH = list(ENH)
names(ENH) = name


SIL = lapply(names(profile_sim_sil), function(params) {
  srr_mat = srr_sil[[params]]
  prof_mat = profile_sim_sil[[params]]
  res = prof_mat*srr_mat
  return(res)
})
names(SIL) = names(profile_sim_sil)
SIL = list(SIL)
names(SIL) = name



info = do.call(rbind, lapply(unique(sel_par_sets$hw), function(hw){
  do.call(rbind, lapply(unique(sel_par_sets$ew), function(ew){
    params       = paste0('hw_',hw,'_ew_',ew)

    summary_file = file.path(DIR_DIAG, paste0('summary_results_', name,'_',params,'.tsv'))
    if(!file.exists(summary_file)) {return()}
    else {
      res = read.delim(summary_file, check.names=FALSE)
    }
    res$key = params
    return(res)
  }))
}))


for ( typeAS in c('enh','sil')) {
  if (typeAS =='sil') {MAT = SIL[[name]]} else {MAT = ENH[[name]]}

  sub_res_df = subset(info, key %in% unique(PARAMS$params))
  if (typeAS =='sil') {
    sub_res_df = sub_res_df[which(sapply(1:nrow(sub_res_df), function(x){"Silenced" %in% sub_res_df[x,] | "Both" %in% sub_res_df[x,]})),]
  } else {
    sub_res_df = sub_res_df[which(sapply(1:nrow(sub_res_df), function(x){"Enhanced" %in% sub_res_df[x,] | "Both" %in% sub_res_df[x,]})),]
  }

  if (dim(sub_res_df)[1]!=0) {
    res = plot_association_heatmap(MAT            = MAT,
                                   PARAMS         = PARAMS,
                                   PEAK           = PEAK,
                                   sub_res_df     = sub_res_df,
                                   typeAS         = typeAS,
                                   input_name     = name,
                                   ranking_method = ranking_method,
                                   sub_deseq      = sub_deseq,
                                   top_mrms       = top_mrms_param,
                                   top_rbps       = top_rbps_param
    )
    p           = res[[1]]
    mat         = res[[2]]
    ranking     = res[[3]]
    tet_score   = res[[4]]
    pval_legend = res[[5]]

    message(noquote(paste0("\n[*] Saving heatmap for ", ifelse(typeAS == 'sil', 'silenced', 'enhanced'), " exons as rds file...")))
    saveRDS(mat, file.path(DIR_OUTPUT, paste0('final_mat_',typeAS,'.rds')))
    saveRDS(tet_score, file.path(DIR_OUTPUT, paste0('final_mat_tet_score_',typeAS,'.rds')))
    saveRDS(ranking, file.path(DIR_OUTPUT, paste0('final_mat_prot_score_',typeAS,'.rds')))

    if ( all(dim(mat) ==c(1,1))) {
      ht_dim = c(10,10)
    } else {
      pdf(NULL)
      ht_dim = getdim(draw(p))
      dev.off()
    }

    # Save plot
    message(noquote(paste0("[*] Printing heatmap for ", ifelse(typeAS == 'sil', "silenced", "enhanced"), " exons...")))
    figure_final_path = file.path(DIR_OUTPUT, paste0(cell_line_target,'_', name, '_', typeAS,'.pdf'))
    pdf(figure_final_path, width = ht_dim[1]+5, height = ht_dim[2])
    draw(p, annotation_legend_list = list(pval_legend))
    add_splicing_maps_to_final_ht(sub_splicingMaps = all_splicingMaps, mat = mat, params = PARAMS, inExon = inExon_param, inIntron = inIntron_param)
    dev.off()
  } else {next}
}
