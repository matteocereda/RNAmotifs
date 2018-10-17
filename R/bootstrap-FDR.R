args = commandArgs(TRUE)

if(length(args)<3) stop("Error: you must specify your working directory, your protein folder, and number of boostraps")


wd = args[1]
pr = args[2]
pp = paste0(pr,"/")

## Number of permutations
#-------------------------
bootstrapN  =  as.numeric(args[3])

ncores = 1
if(length(args)==4){
	ncores=   as.numeric(args[4])
}

# ------------------
## CONFIGURATION
# ------------------
source("config.R")
require(bootstrap)
require(parallel)

# Folders with results of tetramer analysis
wr = c("r/","nr/")
ff = c("fisher-redundant-RSYW-percent-","fisher-not-redundant-percent-")

# Set your working directory here
setwd(wd)

#---------------
# RUN BOOTSTRAPS
#---------------


print.logo()

files     = read.table(paste(pp,wr[1],"filelist_count.tsv",sep="",coll=""),stringsAsFactors=F,header=F)[,1]
tet_names = sapply(strsplit(files,"_"),function(x) x[1])
fNList    = paste0( pp,wr[1], files)

# not 
files     = read.table(paste(pp,wr[2],"filelist_count.tsv",sep="",coll=""),stringsAsFactors=F,header=F)[,1]
tet_names = c(tet_names,sapply(strsplit(files,"_"),function(x) x[1]))
fNList    = c(fNList,paste0( pp,wr[2], files))

if(length(grep("filelist",fNList))>0) fNList = fNList[-grep("filelist",fNList)]

outFName  =  paste0(pp,"bootstrap_",bootstrapN,".tsv")
outFRName =  paste0(pp,"bootstrap_",bootstrapN,".Rdata")


# BOOTSTRAP 

cat("Reading tetramer counts....\n")
fNList = as.list(fNList)
myHitsData=lapply_pb( fNList, function(x){ dd<-read.delim(x); return(dd[,c("type","hits_region1","hits_region2","hits_region3")])})
names(myHitsData) = tet_names

numTet <- length(myHitsData)

if (numTet > 0) {
	bootInd <- 1:dim(myHitsData[[1]])[[1]]

	cat("Calculating enrichment....\n")
	pFisher <- fisher.region(bootInd, myHitsData)

	# reshape orig results and test if correctly reshaped

	pFisherV <- as.vector(pFisher)
	stopifnot(all(matrix(pFisherV, nrow=numTet) == pFisher))

	# bootstrap

	cat("Bootstrapping results....\n")

	# cl = makeCluster(no_cores)

	# tmp = clusterEvalQ(cl, library(bootstrap))
 
 #    bootRes = parLapply(cl, 
 #    	)

	ptm <- proc.time()

	bootRes <- bootstrap(bootInd, bootstrapN, fisher.region.boot, myHitsData)
	
	print(proc.time() - ptm)

	# compare and reshape results; rows correspond to tetramers

	rsM <- matrix(rowSums(bootRes$thetastar<=pFisherV), nrow=numTet)
	pEmpirical<-(1+rsM)/(1+bootstrapN)

	# prepare output

	rownames(pEmpirical) <- rownames(pFisher)
	colnames(pEmpirical) <- paste(colnames(pFisher), "pEmp", sep="_")
	colnames(pFisher) <- paste(colnames(pFisher), "pFis", sep="_")
	cat("Saving results....\n")
	outRes<-data.frame('tetramer'=rownames(pFisher))
	for (cIdx in 1:6) {
		outRes<-cbind(outRes, pFisher[,cIdx], pEmpirical[,cIdx])
 		colnames(outRes)[(ncol(outRes)-1):ncol(outRes)] <- c( colnames(pFisher)[[cIdx]], colnames(pEmpirical)[[cIdx]] )
	}
  	save(outRes, file=outFRName)
	write.table(format(outRes, digits=3), outFName, sep="\t", quote=FALSE, row.names=FALSE)

}else{
  stop("No Tetramer counts")
}
