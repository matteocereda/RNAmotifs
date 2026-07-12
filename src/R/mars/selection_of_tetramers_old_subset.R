#######################################################################################
## Selection of tetramers and RNA maps by Matteo Cereda
#######################################################################################

args = commandArgs(TRUE)

wd = args[1]
pr = args[2]
pp = paste0(pr,"/")

## Number of permutations
#-------------------------

bootstrapN = as.numeric(args[3])

## THRESHOLDS
#--------------

cFisher = as.numeric(args[4])
cEmp    = as.numeric(args[5])
inExon    = as.numeric(args[6])
inIntron  = as.numeric(args[7])

## CONFIGURATION
# ------------------
source("conf/config_RNAmotifs.R")
# Set your working directory here
setwd(wd)

# Folders with results of tetramer analysis
wr = c("r/","nr/")
ff = c("fisher-redundant-RSYW-percent-","fisher-not-redundant-percent-")


region_spacing = 2 * max(inExon, inIntron) + 10
regions = cumsum(c(region_spacing, rep(region_spacing, 3)))

rown = c((regions[1]-inExon):(regions[1]+inIntron),
         (regions[2]-inIntron):(regions[2]+inExon),
         (regions[3]-inExon):(regions[3]+inIntron),
         (regions[4]-inIntron):(regions[4]+inExon))

#---------------
# RUN SELECTION
#---------------
suppressWarnings(library(lattice))
suppressWarnings(library(latticeExtra))
suppressWarnings(library(bootstrap))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))

print.logo()

cat("\nSelecting significantly enriched motifs...\n")

load(paste0(pp,"bootstrap_",bootstrapN,".Rdata"))

res     = outRes

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

p_cutoff = p_cutoff * (p_cutoff < 0.05) + 0.05 * (p_cutoff > 0.05)

sig = subset(res,
             (r1enh_pFis <= p_cutoff[1] & r1enh_pEmp <= cEmp) | (r1sil_pFis <= p_cutoff[1] & r1sil_pEmp <= cEmp) |
               (r2enh_pFis <= p_cutoff[2] & r2enh_pEmp <= cEmp) | (r2sil_pFis <= p_cutoff[2] & r2sil_pEmp <= cEmp) |
               (r3enh_pFis <= p_cutoff[3] & r3enh_pEmp <= cEmp) | (r3sil_pFis <= p_cutoff[3] & r3sil_pEmp <= cEmp)
)

regs = read.csv(paste0(pp,"regions.csv"))
if(ncol(regs) == 2){

  regs = cbind(sig$tetramer,t(as.matrix(regs$x)))
  dimnames(regs) = list(sig$tetramer,c("X",paste0("R",1:3)))
  regs = as.data.frame(regs)
}

tetramers = unique(sig$tetramer)
groups = list(both=regs[,1],
              sil=regs[which(sapply(1:nrow(regs),function(x){"Silenced"%in% regs[x,]|"Both"%in% regs[x,]})),1],
              enh=regs[which(sapply(1:nrow(regs),function(x){"Enhanced"%in% regs[x,]|"Both"%in% regs[x,]})),1])

groups = groups[sapply(groups,function(x){!is.null(x)})]
sig2 = sig

for(gr_ind in names(groups)){
  sig = subset(sig2,tetramer %in% groups[[gr_ind]])


  if(nrow(sig)>0){
    cat("\nSignificantly enriched motifs:\t",nrow(sig),"\n")

    sig[,1]   = as.character(sig[,1])

    tets = unique(sig[,1])

    fileCounts = paste(pp,wr[1],"filelist_count.tsv",sep="",coll="")

    fn = read.delim( paste( pp,wr[1], read.table(fileCounts,stringsAsFactors=F,header=F)[1,1], sep="",coll=""))
    et = table(fn$type)
    CEone    = et["1"]
    CEminone = et["-1"]
    CEzero   = et["0"]
    # calculate fisher at single positions
    cat("\nCalculating fisher's test at single positions...\n")

    id_redundant = grep("W|S|Y|R", tets)

    ll = list(enh=rep(0,length(rown)),sil = rep(0,length(rown)), cont = rep(0,length(rown)))

    for(i in tets){
      for(j in 1:2){
        ll_tmp = getTables(pp,wr[j],i, inExon, inIntron)
        if(length(ll_tmp)>0){
          ll$enh  = ll$enh  + ll_tmp$enh
          ll$sil  = ll$sil  + ll_tmp$sil
          ll$cont = ll$cont + ll_tmp$cont
        }

      }
    }

    p.ena = p.sil = cbind()
    p.ena = cbind(p.ena, lmb.cluster.fisher.test(ll[["enh"]],max(c(CEone,max(ll[["enh"]]))),ll[["cont"]],max(c(CEzero,max(ll[["cont"]])))))
    p.sil = cbind(p.sil,lmb.cluster.fisher.test(ll[["sil"]],max(c(CEminone,max(ll[["sil"]]))),ll[["cont"]],max(c(CEzero,max(ll[["cont"]])))))

    tets = paste0("Group_",gr_ind); colnames(p.ena) = paste0("Group_",gr_ind); colnames(p.sil) = paste0("Group_",gr_ind)

    #------------------------------------------------------
    # FISHER S METHOD : SUM LOG2 PVALUE ENHANCED, SILENCED
    #------------------------------------------------------
    cat("\nCalculating fisher's method...\n")


    ES = matrix(0,nc=length(tets),nr=nrow(p.ena), dimnames=list(rownames(p.ena),tets))

    for( i in tets )	 ES[,i] = (-2)*( log(p.ena[,i])+log( p.sil[,i] ) )

    #----------------------------------
    # sort and cluster tetramers on ES
    #-----------------------------------
    cat("\nSorting and clustering tetramers on ES...\n")

    auc  = apply(ES,2,sum)
    tets = names(sort(auc,dec=T))

    score.sort =  (-2)*log(p.ena) - (-2)*log(p.sil)

    st   = SortingAndType(score.sort,tets)
    if(length(tets)==1) st = list(tets,1)
    ord  = st[[1]]
    rcol = rainbow(max( st[[2]] ))[ st[[2]] ]
    names(rcol) = ord

    #--------------------------------
    # plot RNA MAPS
    #--------------------------------
    cat("Plotting RNA maps...\n")

    suppressWarnings(library(reshape2))

    e = (-2)*log(p.ena)
    s = (-2)*log(p.sil)
    datar=cbind()

    write.csv(e,
              file=paste0(pp,ifelse(gr_ind=="both","All",ifelse(gr_ind=="sil","Sil","Enh")),"_group_enh-",format(Sys.time(),"%Y%m%d_%H_%M"),"-tets-",pr,"-emp-",cEmp,"-fisher-",cFisher,"-nBoot-",bootstrapN,".csv"),
              row.names=F)
    write.csv(s,
              file=paste0(pp,ifelse(gr_ind=="both","All",ifelse(gr_ind=="sil","Sil","Enh")),"_group_sil-",format(Sys.time(),"%Y%m%d_%H_%M"),"-tets-",pr,"-emp-",cEmp,"-fisher-",cFisher,"-nBoot-",bootstrapN,".csv"),
              row.names=F)
    write.csv(e+s,
              file=paste0(pp,ifelse(gr_ind=="both","All",ifelse(gr_ind=="sil","Sil","Enh")),"_group_both-",format(Sys.time(),"%Y%m%d_%H_%M"),"-tets-",pr,"-emp-",cEmp,"-fisher-",cFisher,"-nBoot-",bootstrapN,".csv"),
              row.names=F)

    for( i in ord )	datar = rbind(	datar,  cbind(  "s" = ES[,i], "e" = e[,i], "l" = rown, "g" = rep(i,nrow(ES))))

    datar = as.data.frame(datar,stringsAsFactors=F);rownames(datar)=NULL

    datar[,c("s","e","l")] = apply(datar[,c("s","e","l")],2,as.numeric)

    ms      = max(datar$s)
    datar$z = datar$e/datar$s
    datar$g = factor(datar$g,rev(ord))
    datar$s = datar$s/ms
    datar$z = datar$z/ms

    prn = rnaScoreMap(datar,rcol,ylabels=c("",as.character(floor(ms))), exon = inExon, intron = inIntron, regions = region_spacing)
    prn = update(prn,lwd=0.8)

    suppressWarnings({pdf(file=paste(pp,ifelse(gr_ind=="both","All",ifelse(gr_ind=="sil","Sil","Enh")),"_",length(unique(sig$tetramer)),"_",format(Sys.time(),"%Y%m%d_%H_%M"),"-tets-",pr,"-emp-",cEmp,"-fisher-",cFisher,"-nBoot-",bootstrapN,".pdf",sep="",collapse=""),
        height=10, width = 10)
    print(prn,panel.height=list(0.4,"cm"),panel.width=list(20,"cm"))
    dev.off()
   })

  }else{
    print("NO significant tetramers")
  }
}
