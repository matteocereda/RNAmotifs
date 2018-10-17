#######################################################################################
## Selection of tetramers and RNA maps by Matteo Cereda <matteo.cereda@hugef-torino.org>
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

## CONFIGURATION
# ------------------
source("config.R")
# Set your working directory here
setwd(wd)

# Folders with results of tetramer analysis
wr = c("r/","nr/")
ff = c("fisher-redundant-RSYW-percent-","fisher-not-redundant-percent-")

rown = c(950:1200,1800:2050,2950:3200,3800:4050)

#---------------
# RUN SELECTION
#---------------


print.logo()

cat("\nSelecting significantly enriched motifs...\n")

load(paste0(pp,"bootstrap_",bootstrapN,".Rdata"))

res     = outRes
# res[,2:13]= apply(res[,2:13],2,as.numeric)
# res[,1] = sapply(strsplit(sapply(strsplit(res[,1],"/"),function(x) x[length(x)]),"_region"),function(x) x[1])

sig = subset(res,

	(r1enh_pFis<=cFisher & r1enh_pEmp<=cEmp) |
	(r2enh_pFis<=cFisher & r2enh_pEmp<=cEmp) |
	(r3enh_pFis<=cFisher & r3enh_pEmp<=cEmp) |

	(r1sil_pFis<=cFisher & r1sil_pEmp<=cEmp) |
	(r2sil_pFis<=cFisher & r2sil_pEmp<=cEmp) |
	(r3sil_pFis<=cFisher & r3sil_pEmp<=cEmp)

	)


# x = grep("S",aa[,1]); if(length(x)>0) aa = aa[-x,]
# x = grep("W",aa[,1]); if(length(x)>0) aa = aa[-x,]

# sig  = subset(aa, r1enh_pEmp<=cEmp | r2enh_pEmp<=cEmp | r3enh_pEmp<=cEmp | r1sil_pEmp<=cEmp | r2sil_pEmp<=cEmp | r3sil_pEmp<=cEmp )

if(nrow(sig)>0){
	cat("\nSignificantly enriched motifs:\t",nrow(sig),"\n")

	sig[,1]   = as.character(sig[,1])

	s = matrix(0,nr=2*nrow(sig),nc=11,dimnames=list(NULL,c("tetramer","full_motifs","exonType","is.sign","where.sign","r1_pf","r1_pe","r2_pf","r2_pe","r3_pf","r3_pe")))
	s[,1] = rep(sig[,1],2)
	s[,3] = c(rep("enh",nrow(sig)),rep("sil",nrow(sig)))
	s[1:nrow(sig),c("r1_pf","r1_pe","r2_pf","r2_pe","r3_pf","r3_pe")] = as.matrix(sig[,c("r1enh_pFis","r1enh_pEmp","r2enh_pFis","r2enh_pEmp","r3enh_pFis","r3enh_pEmp")])
	s[(nrow(sig)+1):nrow(s),c("r1_pf","r1_pe","r2_pf","r2_pe","r3_pf","r3_pe")] = as.matrix(sig[,c("r1sil_pFis","r1sil_pEmp","r2sil_pFis","r2sil_pEmp","r3sil_pFis","r3sil_pEmp")])
	s = as.data.frame(s,stringsAsFactors=F)
	s[,4:ncol(s)]=apply(s[,4:ncol(s)],2,as.numeric)
	x1=	apply(s[,c("r1_pe","r2_pe","r3_pe")],1,function(x) x<=cEmp)
	x2=	apply(s[,c("r1_pf","r2_pf","r3_pf")],1,function(x) x<=cFisher)


	s$where.sign = apply(x1+x2,2,function(x) paste(which(x==2),collapse=","))
	s$is.sign=s$where.sign!=""
	tmp = strsplit(s[,1],"")
	s$full_motifs=sapply( tmp, function(x,y) paste(y[x], collapse="") , y=IUPAC_CODE)
	# s$where.sign[which(s$where.sign=='')] = NA
	# sig= subset(s,where!="")

	df_enriched_tetramers=s

	tets = unique(sig[,1])

	fileCounts = paste(pp,wr[1],"filelist_count.tsv",sep="",coll="")

	fn = read.delim( paste( pp,wr[1], read.table(fileCounts,stringsAsFactors=F,header=F)[1,1], sep="",coll=""))
	et = table(fn$type)

	CEone    = et["1"]
	CEminone = et["-1"]
	CEzero   = et["0"]

	# calulate fisher at single positions
	cat("\nCalulating fisher's test at single positions...\n")

	p.ena = p.sil = cbind()
	for( f in 1:2){
		ll = getTables(pp,wr[f],tets)
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

	# e =
	# s = (-2)*log(p.sil)


	#----------------------------------
	# sort and cluster tetramers on ES
	#-----------------------------------
	cat("\nSorting and clustering tetramers on ES...\n")

	auc  = apply(ES,2,sum)
	tets = names(sort(auc,dec=T))

	# ks = matrix(0, nc=length(tets), nr=length(tets), dimnames=list(tets,tets))
	# for(i in tets) for(j in tets) ks[i,j] = ks.test(ES[,i],ES[,j])$p.value
	# ks[which(row(ks)+col(ks)>length(tets))] = NA
	# heatmap(ks)

	score.sort =  (-2)*log(p.ena) - (-2)*log(p.sil)

	# x = melt(score.sort, id.var=colnames(score.sort))
	# 	ggplot(x, aes(x=Var1,y=value, fill=as.numeric(value>0)+1))+geom_area()+scale_fill_manual(values=c("blue",'red'))+facet_grid(~Var2)+theme_bw()

	st   = SortingAndType(score.sort,tets)
	if(length(tets)==1) st = list(tets,1)
	ord  = st[[1]]
	rcol = rainbow(max( st[[2]] ))[ st[[2]] ]
	names(rcol) = ord

	df_ord = do.call(cbind, st)
	s$cluster_id = as.numeric(df_ord[match(s$tetramer,df_ord[,1]),2])

	write.csv(s,
		file=paste0("RNAmotifs/",format(Sys.time(),"%Y%m%d_%H_%M"),"-tets-emp-",cEmp,"-fisher-",cFisher,"-nBoot-10000.csv"),
		row.names=F)

	#--------------------------------
	# plot RNA MAPS
	#--------------------------------
	cat("Plotting RNA maps...\n")

	x = melt(s, id.vars = c("tetramer", "full_motifs", 'cluster_id', 'exonType'),measure.vars = c('r1_pf','r2_pf','r3_pf'))
  x$name=paste0(x$cluster_id," ", x$tetramer," [", x$full_motifs,"]")
  # x$name = factor(x$name, levels=x$name[x$cluster_id])
  x$variable = toupper(x$variable); x$variable = sapply(strsplit(x$variable,"\\_"),"[[",1)
  x$value[which(x$value>0.05)] = NA
  ggplot(x, aes(y=name,x=variable, fill=log10(value)))+geom_tile(color="black")+coord_equal()+facet_grid(~exonType)+scale_fill_gradientn(colours = myPalette(100))+
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.border = element_blank(),  panel.background = element_blank(), axis.text.x = element_text(angle=90, vjust=0.5,hjust=1)) + xlab("") +ylab("")

  x$value[which(x$value>cFisher)] = NA
  ggplot(x, aes(y=name,x=variable, fill=log10(value)))+geom_tile(color="black")+coord_equal()+facet_grid(~exonType)+scale_fill_gradientn(colours = myPalette(100))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.border = element_blank(),  panel.background = element_blank(), axis.text.x = element_text(angle=90, vjust=0.5,hjust=1)) + xlab("") +ylab("")


	df = s
	e = (-2)*log(p.ena)
	s = (-2)*log(p.sil)
	datar=cbind()

	for( i in ord )	datar = rbind(	datar,  cbind(  "s" = ES[,i], "e" = e[,i], "l" = rown, "g" = rep(i,nrow(ES))))

	datar = as.data.frame(datar,stringsAsFactors=F);rownames(datar)=NULL

	datar[,c("s","e","l")] = apply(datar[,c("s","e","l")],2,as.numeric)

	ms      = max(datar$s)
	datar$z = datar$e/datar$s
	datar$g = factor(datar$g,rev(ord))
	datar$s = datar$s/ms
	datar$z = datar$z/ms

	prn = rnaScoreMap(datar,rcol,ylabels=c("",as.character(floor(ms))))
	prn = update(prn,lwd=0.8)

	pdf(file=paste(pp,format(Sys.time(),"%Y%m%d_%H_%M"),"-tets-",pr,"-emp-",cEmp,"-fisher-",cFisher,"-nBoot-",bootstrapN,".pdf",sep="",collapse=""),
		height=10)
	print(prn,panel.height=list(0.4,"cm"),panel.width=list(10,"cm"))
	dev.off()
	save(prn,df_enriched_tetramers, file=paste(pp,format(Sys.time(),"%Y%m%d_%H_%M"),"-tets-",pr,"-emp-",cEmp,"-fisher-",cFisher,"-nBoot-",bootstrapN,".Rdata",sep="",collapse=""))

}else{
	print("NO significant tetrames")
}






















