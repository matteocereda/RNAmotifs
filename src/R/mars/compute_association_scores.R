## ___________________________
## Script name: compute_association_scores.R
## Purpose of script: Given the RNAmotifs output, compute the association score of enriched tetramers and trained RBPs, compared to either HepG2 or K562 cell line
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
  make_option(c("-p", "--cores"), type="character", default=NULL,
              help="Number of cpus", metavar="numerics"),
  make_option(c("-e", "--eclip_path"), type="character", default=NULL,
              help="eCLIP path", metavar="numerics"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="RNAMaRs directory", metavar="character"),
  make_option(c("-s", "--scripts_dir"), type="character", default=NULL,
              help="Path to mars R scripts directory", metavar="character"),
  make_option(c("--in_exon"), type="integer", default=30,
              help="Region extent into exons in bp (default: 30)", metavar="integer"),
  make_option(c("--in_intron"), type="integer", default=300,
              help="Region extent into introns in bp (default: 300)", metavar="integer"),
  make_option(c("-b", "--bootstraps"), type="integer", default=10000,
              help="Number of bootstrap iterations (default: 10000)", metavar="integer")
);



opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


READY_RNAMOTIFS_INPUT = opt$input
name                  = opt$name
DIR_RESULTS_RNAMOTIFS = opt$RNAmotifs_res
repo                  = opt$repository
cell_line             = opt$cell_line
cpus                  = as.numeric(opt$cores)
ECLIP_PEAKS_PATH      = opt$eclip_path
DIR_OUTPUT            = opt$output
scripts_dir           = opt$scripts_dir
bootstrapN            = opt$bootstraps

source(file.path(scripts_dir, 'config_RNAmars.R'))

cat("====================================================",
    "\nRun name              = ", name,
    "\nInput exons file      = ", READY_RNAMOTIFS_INPUT,
    "\nRNAmotifs results     = ", DIR_RESULTS_RNAMOTIFS,
    "\nRNAMaRs folder        = ", repo,
    "\nReference cell line   = ", cell_line,
    "\nCPUs                  = ", cpus,
    "\nRNAMaRs results       = ", DIR_OUTPUT,
    "\neCLIP peaks path      = ", ECLIP_PEAKS_PATH,
    "\n=================================================="
    )

# TRAINING RBPS DEFINITION ================
if (cell_line == "HepG2"){
  rbps = c("HNRNPC", "HNRNPK", "HNRNPU", "NCBP2", "PRPF8", "PTBP1", "QKI", "RBFOX2", "RBM22", "SF3A3", "SF3B4", "SRSF1", "U2AF1", "U2AF2", "UCHL5")
} else if (cell_line == "K562"){
  rbps = c("AGGF1", "EFTUD2", "FXR1", "HNRNPU", "PRPF8", "PTBP1", "PUS1", "RBM15", "SF3B4", "SRSF1", "TARDBP", "U2AF1", "U2AF2")
}

# HEATMAP PARAMETERS =====================
OPT <- read.csv(file.path(repo, "Tables/RNAmotifs_optimal_parameters.csv"), stringsAsFactors = FALSE)
opt_params = OPT[OPT$true_rbp %in% rbps, , drop = F]
colnames(opt_params)[colnames(opt_params)=='true_rbp']<- 'pr'
sel_par_sets <- unique(opt_params[, c("hw", "ew")])
sel_par_sets$hw <- as.numeric(sel_par_sets$hw)
sel_par_sets$ew <- as.numeric(sel_par_sets$ew)

# RNAmaps PARAMETERS =====================
inIntron              = opt$in_intron
inExon                = opt$in_exon
wd                    = 15
r                     = round(wd/2)
inIntron_adj          = inIntron + r
inExon_adj            = inExon + r
l                     = length(-inExon:inIntron)
pos_adj               = c(r + (1:l),
                          max(r + (1:l)) + 2*r+(1:l),
                          max(r + (1:l)  + 2*r+(1:l)) + 2*r+(1:l),
                          max(r + (1:l)  + 2*r+(1:l)  + 2*r+(1:l)) + 2*r+(1:l))



# 1. Recovery rate score ===========
AUC_vec  = readRDS(paste0(repo,'/Rdata/', cell_line, "_AUC_based_metric_downsampling.rds"))
AUC_sil  = vector("list", nrow(sel_par_sets))
AUC_enh  = vector("list", nrow(sel_par_sets))
names(AUC_sil) = paste0('hw_',sel_par_sets$hw, '_ew_', sel_par_sets$ew)
names(AUC_enh) = paste0('hw_',sel_par_sets$hw, '_ew_', sel_par_sets$ew)

message(noquote("\n\n================================="))
message(noquote(paste0("[*] STEP 1: Extracting AUC-based Signal Recovery rate for ", cell_line, " cell line")))
message(noquote("================================="))
for ( y in 1:nrow(sel_par_sets)) {
  params = sel_par_sets[y,]
  hw = params$hw
  ew = params$ew
  index = paste0("hw_", hw, "_ew_", ew)
  message(noquote(paste0("\n[*] Analyzing optimal parameters combination hw: ", hw, " and ew: ", ew)))

  folder_name   = paste0(name, "_hw_", hw, "_ew_", ew)
  path_folder   = file.path(DIR_RESULTS_RNAMOTIFS, folder_name)

  file_enr_tet  = file.path(path_folder, "enriched_tetramers.txt")
  if(!file.exists(file_enr_tet)){next()}

  enr_tetramers = read.table(file_enr_tet)$V1

  regions       = read.csv(file.path(path_folder,"regions.csv"))
  regions       = adj_reg_mat(regions, enr_tetramers)

  enh_tets = regions$X[apply(regions[2:4], 1, function(row) {'Enhanced'%in% row | 'Both' %in% row})]
  sil_tets = regions$X[apply(regions[2:4], 1, function(row) {'Silenced'%in% row | 'Both' %in% row})]

  rbps_opt_params = opt_params$pr[opt_params$hw == hw & opt_params$ew == ew]

  # AUC_weight matrix
  AUC_enh[[index]] = matrix( rep( AUC_vec[paste0(rbps_opt_params,'_enh')], length(enh_tets)), nrow = length(rbps_opt_params), dimnames = list(rbps_opt_params, enh_tets))
  AUC_sil[[index]] = matrix( rep( AUC_vec[paste0(rbps_opt_params,'_sil')], length(sil_tets)), nrow = length(rbps_opt_params), dimnames = list(rbps_opt_params, sil_tets))
}

message(noquote("\n[*] Saving signal recovery rate files (enhanced and silenced separately) ..."))
saveRDS(AUC_enh, paste0(DIR_OUTPUT, "SCORE1_enh_signal_recovery_rate_", name, ".rds"))
saveRDS(AUC_sil, paste0(DIR_OUTPUT, "SCORE1_sil_signal_recovery_rate_", name, ".rds"))


# 2. Binding scores from RNAmotifs results ===================
SIM_SIL = vector("list", nrow(sel_par_sets))
SIM_ENH = vector("list", nrow(sel_par_sets))
names(SIM_SIL) = paste0('hw_',sel_par_sets$hw, '_ew_', sel_par_sets$ew)
names(SIM_ENH) = paste0('hw_',sel_par_sets$hw, '_ew_', sel_par_sets$ew)
message(noquote("\n\n================================="))
message(noquote("[*] STEP 2: Computing cosine similarity between RBP binding profile and tetramer enrichment score"))
message(noquote("================================="))
res = data.frame()
for (y in 1:nrow(sel_par_sets)) {
  params       = sel_par_sets[y, ]
  hw = params$hw
  ew = params$ew
  index =  paste0('hw_',hw, '_ew_', ew)
  message(noquote(paste0("\n[*] Analyzing optimal parameters combination hw: ", hw, " and ew: ", ew)))

  summary_file = paste0(DIR_OUTPUT, "summary_results_", name, "_hw_", hw, "_ew_", ew, ".rds")

  folder_name   = paste0(name, "_hw_", hw, "_ew_", ew)
  path_folder   = file.path(DIR_RESULTS_RNAMOTIFS, folder_name)
  file_enr_tet  = file.path(path_folder, "enriched_tetramers.txt")

  if(!file.exists(file_enr_tet)){next()}
  load(file.path(path_folder, paste0("bootstrap_", bootstrapN, ".Rdata")))

  enr_tetramers = read.table(file_enr_tet)$V1

  regions       = read.csv(paste0(path_folder,"/regions.csv"))
  regions       = adj_reg_mat(regions, enr_tetramers)

  ## Save RNAmotifs enrichments and pvalues in a convenient way
  res = cbind(data.frame(pr        = name,
                         hw        = params[1],
                         ew        = params[2],
                         tetramer  = enr_tetramers,
                         row.names = NULL),
              regions[match(enr_tetramers, regions$X), c("R1", "R2", "R3")],
              outRes[enr_tetramers, paste0(rep(c("r1", "r2", "r3"), 2), c(rep("enh", 3), rep("sil", 3)), "_pFis")])
  saveRDS(res, summary_file)


  enh_tets = regions$X[apply(regions[2:4], 1, function(row) {'Enhanced'%in% row | 'Both' %in% row})]
  sil_tets = regions$X[apply(regions[2:4], 1, function(row) {'Silenced'%in% row | 'Both' %in% row})]

  # extracting only last run in case of multiple runs with the same name
  ff_silenced = list.files(path_folder,"^s-")[length(list.files(path_folder,"^s-"))]
  ff_enhanced = list.files(path_folder,"^e-")[length(list.files(path_folder,"^e-"))]
  RNAmaps_silenced       = read.csv(paste0(path_folder,"/",ff_silenced))
  RNAmaps_enhanced       = read.csv(paste0(path_folder,"/",ff_enhanced))

  # Input RNAmotifs (exons)
  datExpr    <<- read.table(READY_RNAMOTIFS_INPUT, header = F, sep = ";")

  # Extract binding profile of all the trained proteins RBPs, restricting eCLIP to input exons containing enriched tetramer
  rbps_opt_params = opt_params$pr[opt_params$hw == hw & opt_params$ew == ew]
  message(noquote(paste0("[*] RBPs ", paste(rbps_opt_params, collapse = ", ")," have the combination hw=", hw," ew=", ew, " as optimal parameters")))

  sim_sil_tmp = matrix(NA, nrow = length(rbps_opt_params), ncol = length(sil_tets), dimnames = list(rbps_opt_params, sil_tets))
  sim_enh_tmp = matrix(NA, nrow = length(rbps_opt_params), ncol = length(enh_tets), dimnames = list(rbps_opt_params, enh_tets))

  for(pr in rbps_opt_params) {

    message(noquote(paste0("\n[*] ", pr)))
    message(noquote(paste0("[*] Reading ", cell_line, " eCLIP peaks and overlapping with exons containing enriched tetramers...")))
    peaks      <<- read.delim(paste0(ECLIP_PEAKS_PATH, grep("ordered_merged.bed$",list.files(pattern = pr, path = ECLIP_PEAKS_PATH),value = T)[1]),header=F)
    Peaks       =  SEeCLIPpeaks(datExpr, peaks, inIntron_adj, inExon_adj, cpus)


    ## ________________________________________________________________
    ## > Similarity in ENHANCED profiles restricted to tetramer tetr -----
    if (length(enh_tets)>=1) {
      message(noquote("[*] Computing cosine similarity in enhanced profiles..."))
      similarity_enh = sapply(enh_tets, function(tetr) {
        ind         = which(subset(regions, X==tetr)[,2:4]!=0)
        # index of input exons that have enrichment in R1 or R2 or R3
        ind_datExpr = sort(unique(unlist(lapply(ind,function(w){
          get(paste0("get_exons_with_tetramer_region",w))(tetr, path_folder, READY_RNAMOTIFS_INPUT)$myRID
        } ))))

        # Extracting enhanced and constitutive exons
        CEs_tetr    = Peaks$values[intersect(ind_datExpr, which(datExpr$V9=='1')), ]
        CEs_const   = Peaks$values[which(datExpr$V9 == 0),]

        if(is.vector(CEs_tetr)){
          tmp = as.matrix(CEs_tetr > 0)
        } else {
          tmp = as.matrix(colSums(CEs_tetr > 0))
        }

        bind_profile_tetr = -2*log(lmb.cluster.fisher.test(tmp,                   # alternative exons with tetramer regions[z,] overlapping eCLIP at single position
                                                           length(ind_datExpr),   # total exons with tetramer regions[z,]
                                                           as.matrix(colSums(CEs_const > 0)),
                                                           sum(datExpr$V9 == 0)))

        x1 = bind_profile_tetr[pos_adj]
        x2 = RNAmaps_enhanced[,tetr]
        res = cosine(x1, x2)
        return(res)
      })
      sim_enh_tmp[pr, names(similarity_enh)] = unname(similarity_enh)
    } else {sim_enh_tmp = matrix(NA, nrow = length(rbps_opt_params), ncol = 0, dimnames = list(rbps_opt_params))}

    ## __________________________________________________________________
    ## > Similarity in SILENCED profiles restricted to tetramer tetr -----
    if (length(sil_tets)>=1) {
      message(noquote("[*] Computing cosine similarity in silenced profiles..."))
      similarity_sil = sapply(sil_tets, function(tetr) {
        ind         = which(subset(regions, X==tetr)[,2:4]!=0)
        # index of input exons that have enrichment in R1 or R2 or R3
        ind_datExpr = sort(unique(unlist(lapply(ind,function(w){
          get(paste0("get_exons_with_tetramer_region",w))(tetr, path_folder, READY_RNAMOTIFS_INPUT)$myRID
        } )
        )))

        # Extracting silenced and constitutive exons
        CEs_tetr    = Peaks$values[intersect(ind_datExpr, which(datExpr$V9=='-1')), ]
        CEs_const   = Peaks$values[which(datExpr$V9 == 0),]

        if(is.vector(CEs_tetr)){
          tmp = as.matrix(CEs_tetr > 0)
        } else {
          tmp = as.matrix(colSums(CEs_tetr > 0))
        }

        bind_profile_tetr = -2*log(lmb.cluster.fisher.test(tmp,                   # alternative exons with tetramer regions[z,] overlapping eCLIP at single position
                                                           length(ind_datExpr),   # total exons with tetramer regions[z,]
                                                           as.matrix(colSums(CEs_const > 0)),
                                                           sum(datExpr$V9 == 0)))

        x1 = bind_profile_tetr[pos_adj]
        x2 = RNAmaps_silenced[,tetr]
        res = cosine(x1, x2)
        return(res)
      })
      sim_sil_tmp[pr, names(similarity_sil)] = unname(similarity_sil)
    } else {sim_sil_tmp = matrix(NA, nrow = length(rbps_opt_params), ncol = 0, dimnames = list(rbps_opt_params))}


  }
  SIM_SIL[[index]] = sim_sil_tmp
  SIM_ENH[[index]] = sim_enh_tmp
  message(noquote("\n ==============================="))
}

message(noquote("\n\n[*] Saving cosine similarities (enhanced and silenced separately) ..."))
saveRDS(SIM_SIL, paste0(DIR_OUTPUT, "SCORE2_sil_profile_similarity_", name, ".rds"))
saveRDS(SIM_ENH, paste0(DIR_OUTPUT, "SCORE2_enh_profile_similarity_", name, ".rds"))
