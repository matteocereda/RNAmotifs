suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))

args = commandArgs(TRUE)

pr_folder = paste0("/",args[2])
res_dir = args[1]
bootstrapN = if (length(args) >= 3) as.numeric(args[3]) else 10000
pEmp = if (length(args) >= 4) as.numeric(args[4]) else 0.0005

load(paste0(res_dir,pr_folder,"/bootstrap_",bootstrapN,".Rdata"))
enr_tets <- as.character(read.table(paste0(res_dir,pr_folder,"/enriched_tetramers.txt"))$V1)

p_enh = outRes %>% select(c("r1enh_pFis","r1enh_pEmp","r2enh_pFis","r2enh_pEmp","r3enh_pFis","r3enh_pEmp"))
p_sil = outRes %>% select(c("r1sil_pFis","r1sil_pEmp","r2sil_pFis","r2sil_pEmp","r3sil_pFis","r3sil_pEmp"))

p_cutoff =  p_enh %>% plyr::rename(c("r1enh_pFis" = "r1sil_pFis",
                               "r1enh_pEmp" = "r1sil_pEmp",
                               "r2enh_pFis" = "r2sil_pFis",
                               "r2enh_pEmp" = "r2sil_pEmp",
                               "r3enh_pFis" = "r3sil_pFis",
                               "r3enh_pEmp" = "r3sil_pEmp")) %>%
  bind_rows(p_sil) %>%
  select(c("r1sil_pFis",
           "r2sil_pFis",
           "r3sil_pFis")) %>%
  plyr::rename(c("r1sil_pFis" = "r1_pFis_0.01",
           "r2sil_pFis" = "r2_pFis_0.01",
           "r3sil_pFis" = "r3_pFis_0.01")) %>%
  apply(2,quantile,0.01)

## Region-specific Fisher cutoff = min(1st-percentile, 0.05) (matches selection_of_tetramers.R).
p_cutoff = p_cutoff * (p_cutoff < 0.05) + 0.05 * (p_cutoff > 0.05)

res = outRes

res$is.sign = (

               (res$r1enh_pFis <= p_cutoff[1]  & res$r1enh_pEmp <= pEmp) |
               (res$r2enh_pFis <= p_cutoff[2]  & res$r2enh_pEmp <= pEmp) |
               (res$r3enh_pFis <= p_cutoff[3]  & res$r3enh_pEmp <= pEmp) |

               (res$r1sil_pFis <= p_cutoff[1]  & res$r1sil_pEmp <= pEmp) |
               (res$r2sil_pFis <= p_cutoff[2]  & res$r2sil_pEmp <= pEmp) |
               (res$r3sil_pFis <= p_cutoff[3]  & res$r3sil_pEmp <= pEmp)

)

res = res %>% subset(is.sign)

res$exonType = NA

enh_ind = which(  (res$r1enh_pFis <= p_cutoff[1]  & res$r1enh_pEmp <= pEmp) |
                  (res$r2enh_pFis <= p_cutoff[2]  & res$r2enh_pEmp <= pEmp) |
                  (res$r3enh_pFis <= p_cutoff[3]  & res$r3enh_pEmp <= pEmp))

sil_ind = which(  (res$r1sil_pFis <= p_cutoff[1]  & res$r1sil_pEmp <= pEmp) |
                  (res$r2sil_pFis <= p_cutoff[2]  & res$r2sil_pEmp <= pEmp) |
                  (res$r3sil_pFis <= p_cutoff[3]  & res$r3sil_pEmp <= pEmp) )

res$exonType[enh_ind] = "enh"
res$exonType[sil_ind] = "sil"
if(length(intersect(enh_ind,sil_ind))>0){
  res$exonType[intersect(enh_ind,sil_ind)] = "both"
}

res$R1 = "0"
res$R1[    which( res$r1enh_pFis <= p_cutoff[1] & res$r1enh_pEmp <= pEmp)    ] = "Enhanced"
res$R1[    which( res$r1sil_pFis <= p_cutoff[1] & res$r1sil_pEmp <= pEmp)    ] = "Silenced"
both_ind = which((res$r1enh_pFis <= p_cutoff[1] & res$r1enh_pEmp <= pEmp) & (res$r1sil_pFis <= p_cutoff[1]  & res$r1sil_pEmp <= pEmp))
if(length(both_ind)>0){
  res$R1[both_ind] = "Both"
}

res$R2 = "0"
res$R2[    which( res$r2enh_pFis <= p_cutoff[2] & res$r2enh_pEmp <= pEmp)    ] = "Enhanced"
res$R2[    which( res$r2sil_pFis <= p_cutoff[2] & res$r2sil_pEmp <= pEmp)    ] = "Silenced"
both_ind = which((res$r2enh_pFis <= p_cutoff[2] & res$r2enh_pEmp <= pEmp) & (res$r2sil_pFis <= p_cutoff[2]  & res$r2sil_pEmp <= pEmp))
if(length(both_ind)>0){
  res$R2[both_ind] = "Both"
}

res$R3 = "0"
res$R3[    which( res$r3enh_pFis <= p_cutoff[3] & res$r3enh_pEmp <= pEmp)    ] = "Enhanced"
res$R3[    which( res$r3sil_pFis <= p_cutoff[3] & res$r3sil_pEmp <= pEmp)    ] = "Silenced"
both_ind = which((res$r3enh_pFis <= p_cutoff[3] & res$r3enh_pEmp <= pEmp) & (res$r3sil_pFis <= p_cutoff[3]  & res$r3sil_pEmp <= pEmp))
if(length(both_ind)>0){
  res$R3[both_ind] = "Both"
}

if (dim(res)[1]>1) {
  x = as.matrix(res[,c("R1","R2","R3")])
  x = x[enr_tets,]
} else {
  x = res[enr_tets,c("R1","R2","R3")]
}

write.csv(x,file = paste0(res_dir,pr_folder,"/regions.csv"))
suppressWarnings({pdf(file=paste0(res_dir,pr_folder,"/regions.pdf"))
plot(x, col = c("purple","blue","red"), breaks = c("Both", "Silenced","Enhanced"), axis.row=NULL, axis.col=NULL, ann=F, asp=T,yaxt="n")
dev.off()
})
