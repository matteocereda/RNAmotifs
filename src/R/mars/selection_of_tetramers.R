#######################################################################################
## Selection of tetramers and RNA maps by Matteo Cereda
#######################################################################################

args = commandArgs(TRUE)

if(length(args)<5) stop("Error: you must specify your working directory, your protein folder, number of boostraps,pFisher cutoff and pEmpirical cutoff")


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

p_enh = outRes %>% dplyr::select(c("r1enh_pFis","r1enh_pEmp","r2enh_pFis","r2enh_pEmp","r3enh_pFis","r3enh_pEmp"))
p_sil = outRes %>% dplyr::select(c("r1sil_pFis","r1sil_pEmp","r2sil_pFis","r2sil_pEmp","r3sil_pFis","r3sil_pEmp"))

p_cutoff =  p_enh %>% plyr::rename(c("r1enh_pFis" = "r1sil_pFis",
                                     "r1enh_pEmp" = "r1sil_pEmp",
                                     "r2enh_pFis" = "r2sil_pFis",
                                     "r2enh_pEmp" = "r2sil_pEmp",
                                     "r3enh_pFis" = "r3sil_pFis",
                                     "r3enh_pEmp" = "r3sil_pEmp")) %>%
  bind_rows(p_sil) %>%
  dplyr::select(c("r1sil_pFis",
                  "r2sil_pFis",
                  "r3sil_pFis")) %>%
  plyr::rename(c("r1sil_pFis" = "r1_pFis_0.01",
                 "r2sil_pFis" = "r2_pFis_0.01",
                 "r3sil_pFis" = "r3_pFis_0.01")) %>%
  apply(2,quantile,0.01)

## Region-specific Fisher cutoff = min(1st-percentile of the pFis distribution, 0.05):
##   keep pFis <= 0.05, OR <= the 1st percentile when that 1st percentile is <= 0.05.
p_cutoff = p_cutoff * (p_cutoff < 0.05) + 0.05 * (p_cutoff > 0.05)
print(p_cutoff)
## Use <= (not <): keeps the most significant motifs when the cutoff collapses to 0
## (wide windows where >=1% of motifs saturate the Fisher test at pFis==0).
sig = subset(res,
               (r1enh_pFis <= p_cutoff[1] & r1enh_pEmp <= cEmp) | (r1sil_pFis <= p_cutoff[1] & r1sil_pEmp <= cEmp) |
               (r2enh_pFis <= p_cutoff[2] & r2enh_pEmp <= cEmp) | (r2sil_pFis <= p_cutoff[2] & r2sil_pEmp <= cEmp) |
               (r3enh_pFis <= p_cutoff[3] & r3enh_pEmp <= cEmp) | (r3sil_pFis <= p_cutoff[3] & r3sil_pEmp <= cEmp)
)


if(nrow(sig)>0){
	cat("\nSignificantly enriched motifs:\t",nrow(sig),"\n")

	sig[,1]   = as.character(sig[,1])

	s = matrix(0,nr=2*nrow(sig),nc=10,dimnames=list(NULL,c("tetramer","exonType","is.sign","where.sign","r1_pf","r1_pe","r2_pf","r2_pe","r3_pf","r3_pe")))
	s[,1] = rep(sig[,1],2)
	s[,2] = c(rep("enh",nrow(sig)),rep("sil",nrow(sig)))
	s[1:nrow(sig),c("r1_pf","r1_pe","r2_pf","r2_pe","r3_pf","r3_pe")] = as.matrix(sig[,c("r1enh_pFis","r1enh_pEmp","r2enh_pFis","r2enh_pEmp","r3enh_pFis","r3enh_pEmp")])
	s[(nrow(sig)+1):nrow(s),c("r1_pf","r1_pe","r2_pf","r2_pe","r3_pf","r3_pe")] = as.matrix(sig[,c("r1sil_pFis","r1sil_pEmp","r2sil_pFis","r2sil_pEmp","r3sil_pFis","r3sil_pEmp")])
	s = as.data.frame(s,stringsAsFactors=F)
	s[,3:ncol(s)]=apply(s[,3:ncol(s)],2,as.numeric)
	x1=	apply(s[,c("r1_pe","r2_pe","r3_pe")],1,function(x) x<=cEmp)
	x2=	apply(s[,c("r1_pf","r2_pf","r3_pf")],1,function(x) x<=cFisher)


	s$where.sign = apply(x1+x2,2,function(x) paste(which(x==2),collapse=","))
	s$is.sign=s$where.sign!=""

	df_enriched_tetramers=s

	tets = unique(sig[,1])

	fileCounts = paste(pp,wr[1],"filelist_count.tsv",sep="",coll="")

	fn = read.delim( paste( pp,wr[1], read.table(fileCounts,stringsAsFactors=F,header=F)[1,1], sep="",coll=""))
	et = table(fn$type)

	CEone    = et["1"]
	CEminone = et["-1"]
	CEzero   = et["0"]

	# calculate fisher at single positions
	cat("\nCalculating fisher's test at single positions...\n")

	p.ena = p.sil = cbind()
	for( f in 1:2){
		ll = getTables(pp,wr[f],tets, inExon, inIntron)
		if(length(ll)>0){
			p.ena = cbind(p.ena,lmb.cluster.fisher.test(ll[["enh"]],CEone,ll[["cont"]],CEzero))
			p.sil = cbind(p.sil,lmb.cluster.fisher.test(ll[["sil"]],CEminone,ll[["cont"]],CEzero))
		}
	}


	#------------------------------------------------------
	# FISHER S METHOD : SUM LOG2 PVALUE ENHANCED, SILENCED
	#------------------------------------------------------
	cat("\nCalulating fisher's method...\n")


	ES = matrix(0,nc=length(tets),nr=nrow(p.ena), dimnames=list(rownames(p.ena),tets))

	for( i in tets )	 ES[,i] = (-2)*( log(p.ena[,i])+log( p.sil[,i] ) )

	e = (-2)*log(p.ena)
	s = (-2)*log(p.sil)


	#----------------------------------
	# sort and cluster tetramers on ES
	#-----------------------------------
	cat("\nSorting and clustering tetramers on ES...\n")

	auc  = apply(ES,2,sum)
	tets = names(sort(auc,dec=T))

	score.sort = e-s

	###
 write.csv(ES,
           file=paste0(pp,"ES-",format(Sys.time(),"%Y%m%d_%H_%M"),"-tets-",pr,"-emp-",cEmp,"-fisher-",cFisher,"-nBoot-",bootstrapN,".csv"),
           row.names=F)
 write.csv(e,
           file=paste0(pp,"e-",format(Sys.time(),"%Y%m%d_%H_%M"),"-tets-",pr,"-emp-",cEmp,"-fisher-",cFisher,"-nBoot-",bootstrapN,".csv"),
           row.names=F)
 write.csv(s,
           file=paste0(pp,"s-",format(Sys.time(),"%Y%m%d_%H_%M"),"-tets-",pr,"-emp-",cEmp,"-fisher-",cFisher,"-nBoot-",bootstrapN,".csv"),
           row.names=F)
	###

	st   = SortingAndType(score.sort,tets)
	if(length(tets)==1) st = list(tets,1)
	ord  = st[[1]]
	rcol = rainbow(max( st[[2]] ))[ st[[2]] ]
	names(rcol) = ord

	df_ord = do.call(cbind, st)
	df_enriched_tetramers$cluster_id = as.numeric(df_ord[match(df_enriched_tetramers$tetramer,df_ord[,1]),2])

	write.csv(df_enriched_tetramers,
		file=paste0(pp,format(Sys.time(),"%Y%m%d_%H_%M"),"-tets-",pr,"-emp-",cEmp,"-fisher-",cFisher,"-nBoot-",bootstrapN,".csv"),
		row.names=F)

	#--------------------------------
	# plot RNA MAPS
	#--------------------------------
	cat("Plotting RNA maps...\n")

	datar=cbind()
	for( i in ord )	datar = rbind(	datar,  cbind(  "s" = ES[,i], "e" = e[,i], "l" = rown, "g" = rep(i,nrow(ES))))

	datar = as.data.frame(datar,stringsAsFactors=F);rownames(datar)=NULL

	datar[,c("s","e","l")] = apply(datar[,c("s","e","l")],2,as.numeric)

	ms      = max(datar$s)
	datar$z = datar$e/datar$s
	datar$g = factor(datar$g,rev(ord))
	datar$s = datar$s/ms
	datar$z = datar$z/ms

	###

	write.table(df_ord[,1], file = paste0(pp,"enriched_tetramers.txt"), sep="\n", quote=F, row.names=F, col.names=F)
	write.table(df_ord, file = paste0(pp,"cluster_ids.txt"), sep = "\n", quote = F, row.names = F, col.names = F)

	prn = rnaScoreMap(datar,rcol,ylabels=c("",as.character(floor(ms))), exon = inExon, intron = inIntron, regions = region_spacing)
	prn = update(prn,lwd=0.8)

	suppressWarnings({pdf(file=paste(pp,format(Sys.time(),"%Y%m%d_%H_%M"),"-tets-",pr,"-emp-",cEmp,"-fisher-",cFisher,"-nBoot-",bootstrapN,".pdf",sep="",collapse=""),
		height=10, width = 10)
	print(prn,panel.height=list(0.4,"cm"),panel.width=list(10,"cm"))
	dev.off()
        })
	save(prn, file=paste(pp,format(Sys.time(),"%Y%m%d_%H_%M"),"-tets-",pr,"-emp-",cEmp,"-fisher-",cFisher,"-nBoot-",bootstrapN,".Rdata",sep="",collapse=""))

}else{
	print("NO significant tetramers")
}
