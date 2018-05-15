
BOOTSTRAP_FILE="~/UV2_remote/matteo/RNAmotifs/results/NOVA2/bootstrap_10000.Rdata"
SPLICING_FILE ="AS_events/NOVA2.SE.rnamotifs.txt"
RESULT_DIR="~/UV2_remote//matteo/RNAmotifs/results/NOVA2/"
XLSX_FILE="NOVA2.RNAmotifs.xlsx"
P.EMP=0.05
P.FISH=0.01

source("~/UV2_remote/matteo/RNAmotifs/R/config.R")

load(BOOTSTRAP_FILE)

mats=read.table(SPLICING_FILE, sep=";")

nr = lapply(list.files(paste0(RESULT_DIR,"nr"), "_region_count.tsv", full.names = T), read.delim )
names(nr) = gsub( "_region_count.tsv", "", list.files(paste0(RESULT_DIR,"nr"), "_region_count.tsv"))

r = lapply(list.files(paste0(RESULT_DIR,"r"), "_region_count.tsv", full.names = T), read.delim)
names(r) = gsub( "_region_count.tsv", "", list.files(paste0(RESULT_DIR,"r"), "_region_count.tsv"))

xx = c(nr, r)
snr = lapply(xx, function(x) ddply(x, .(type), summarise,
                                   R1=sum(hits_region1),
                                   R2=sum(hits_region2),
                                   R3=sum(hits_region3),
                                   tot=sum(hits_region1)+sum(hits_region2)+sum(hits_region3)))

s = mapply(function(x,y){ x$tet=y; return(x)}, x=snr, names(snr), SIMPLIFY = F )
s = do.call(rbind, s)
m=s
t = table(mats$V9)

m$n_exons = t[as.character(m$type)]

mt=melt(m, id.vars = c('tet','type','tot','n_exons') )

ct = subset(mt, type==0); ct$key=paste0(ct$tet,".",ct$variable)
es = subset(mt, type!=0); es$key=paste0(es$tet,".",es$variable)

es$n_exons.ct = ct$n_exons[match(es$key,ct$key)]
es$value.ct = ct$value[match(es$key,ct$key)]

mf = ddply(es, .(tet,type,variable), mutate,
           pv = fisher.test(matrix(c(value, n_exons,
                                     value.ct, n_exons.ct),nr=2, byrow=T), alternative = 'g')$p.value
)

boot = melt(outRes, id.vars = 'tetramer')

nn = c('r1enh_pFis' ,'r1enh_pEmp' ,'r1sil_pFis', 'r1sil_pEmp', 'r2enh_pFis' ,'r2enh_pEmp', 'r2sil_pFis', 'r2sil_pEmp', 'r3enh_pFis', 'r3enh_pEmp', 'r3sil_pFis', 'r3sil_pEmp')
rr = c(rep('R1',4),rep('R2',4),rep('R3',4))
tt = rep(c(1,1,-1,-1),3)
te = rep(c('Fisher',"Emp",'Fisher',"Emp"),3)

boot$region=rr[boot$variable]
boot$type=tt[boot$variable]
boot$test=te[boot$variable]

names(rr)=nn
names(tt)=nn
names(te)=nn


mf$key=paste0(mf$tet,".",mf$variable,'.',mf$type)
boot$key=paste0(boot$tetramer,".",boot$region,'.',boot$type)

mf$Emp=with(subset(boot, test=='Emp'), value[match(mf$key, key)])

write.xlsx(mf, file=XLSX_FILE, row.names = F)


e = subset(s, type==1)
c = subset(s, type==(0))
s = subset(s, type==(-1))

smf = subset(mf, pv<=P.FISH & Emp<=P.EMP)

# RNAMAPS
pp = RESULT_DIR #"~/UV2_remote/matteo/RNAmotifs/results/NOVA2/"

fileCounts = list.files(pp, "_region_count.tsv", recursive = T, full.names = T)
bed        = list.files(pp, "bed", recursive = T, full.names = T)

nf = gsub("_region_count.tsv","",sapply(strsplit(fileCounts, "\\/"), function(x) x[length(x)]))
nb = gsub(".bed","",sapply(strsplit(bed, "\\/"), function(x) x[length(x)]))

tet =  lapply(smf$tet, function(x,y) grep(x,y), y=nf) ; ids  = unique(unlist(tet))
tetb = lapply(smf$tet, function(x,y) grep(x,y), y=nb) ; idsb = unique(unlist(tetb))

counts = lapply(fileCounts[ids], read.delim); names(counts)=nf[ids]
et = table(counts[[1]]$type)
CEone    = et["1"]
CEminone = et["-1"]
CEzero   = et["0"]

tet =  lapply(smf$tet, function(x,y) grep(x,y), y=nf[ids])
names(tet) = smf$tet

fbed = lapply(bed[idsb], read.delim,fill=TRUE,skip=1,header=F); names(fbed)=nb[idsb]
fbed = lapply(fbed, function(x,y){colnames(x)=y; return(x)}, y= c("pos","cat","no"))
fbed = mapply(function(x,y) {x$tet=y; return(x)}, fbed, names(fbed), SIMPLIFY = F)

# counts = mapply(function(x,y) {x$tet=y; return(x)}, counts, names(counts), SIMPLIFY = F)

tet.fisher = function(a,b,c,d){
  unlist(mapply(
    function(x,y,w,z)  fisher.test(matrix( data=c(x, w-x, y,z-y ), nrow=2, byrow=T))$p.value,
    as.list(a),as.list(b),
    c,d, SIMPLIFY = F))
}

p.ena = p.sil = cbind()
for(x in 1:length(tet)){
  print(names(tet)[x])
  dd = as.data.frame(fbed[tet[[x]]]);
  colnames(dd)= c('pos','cat', 'no', 'tet')
  merco = ddply(dd, .(pos,cat), summarise,
                no = sum(no, na.rm=T)
  )
  l = dlply(merco, ~cat)
  l = lapply(l, function(x){rownames(x)=x[,1]; return(x)})

  tab = matrix( 0, nrow=1004, ncol=1 )

  ll = list( 'enh'  = l[['1']][ match(rown, rownames(l[['1']])),"no"],
             'sil'  = l[['-1']][ match(rown, rownames(l[['-1']])),'no'],
             'cont' = l[['0']][ match(rown, rownames(l[['0']])),'no']
  )

  for(i in 1:3) ll[[i]][which(is.na(ll[[i]]))]=0
  for(i in 1:3) names(ll[[i]]) = rown

  # calulate fisher at single positions
  p.ena = cbind(p.ena,tet.fisher(ll[["enh"]],ll[["cont"]],CEone,CEzero))
  p.sil = cbind(p.sil,tet.fisher(ll[["sil"]],ll[["cont"]],CEminone,CEzero))
}

tets = names(tet)
ES = matrix(0,nc=length(tets),nr=nrow(p.ena), dimnames=list(rownames(p.ena),tets))
for( i in 1:length(tet) )	 ES[,i] = (-2)*( log(p.ena[,i])+log( p.sil[,i] ) )


auc  = apply(ES,2,sum)
tets = names(sort(auc,dec=T))

score.sort =  (-2)*log(p.ena) - (-2)*log(p.sil)
colnames(score.sort) = names(tet)
st   = SortingAndType(score.sort,tets)

if(length(tets)==1) st = list(tets,1)
ord  = st[[1]]
rcol = rainbow(max( st[[2]] ))[ st[[2]] ]
names(rcol) = ord

df_ord = do.call(cbind, st)
smf$cluster_id = as.numeric(df_ord[match(smf$tet,df_ord[,1]),2])

write.xlsx(smf, file="NOVA2.RNAmotifs.xlsx", "Sig Tets", append=T, row.names = F)

# plot RNA MAPS
cat("Plotting RNA maps...\n")

e = (-2)*log(p.ena)
s = (-2)*log(p.sil)
colnames(e) = names(tet)
colnames(s) = names(tet)

datar=cbind()
for( i in smf$tet )
  datar = rbind(	datar,  cbind(  "s" = ES[,i], "e" = e[,i], "l" = rown, "g" = rep(i,nrow(ES))))

datar = as.data.frame(datar,stringsAsFactors=F);rownames(datar)=NULL

datar[,c("s","e","l")] = apply(datar[,c("s","e","l")],2,as.numeric)


ms      = max(datar$s)
datar$z = datar$e/datar$s
datar$g = factor(datar$g,rev(ord))
datar$s = datar$s/ms
datar$z = datar$z/ms

prn = rnaScoreMap(datar,rcol,ylabels=c("",as.character(floor(ms))))
prn = update(prn,lwd=0.8)

pdf(file=paste0("Figure/RNAmaps-",format(Sys.time(),"%Y%m%d_%H_%M"),"-tets-NOVA2-emp-",P.EMP,"-fisher-",P.FISH,"-nBoot-10000.pdf"), height=10)
    print(prn,panel.height=list(0.4,"cm"),panel.width=list(10,"cm"))
dev.off()

# save(prn,df_enriched_tetramers, file=paste(pp,format(Sys.time(),"%Y%m%d_%H_%M"),"-tets-",pr,"-emp-",cEmp,"-fisher-",cFisher,"-nBoot-",bootstrapN,".Rdata",sep="",collapse=""))








cEmp=10^(-4) # minimum with 10,000 simulation
cFisher=10^(-8)

# Summary pvalue =
s = apply(s[,2:ncol(outRes)],2, quantile, probs=seq(0,1,.05) )
ds = melt(s, id.vars=rownames(s))
ds$region = substr(ds$Var2,1,2)
ds$as = substr(ds$Var2,3,5)
ds$method=sapply(strsplit(as.character(ds$Var2),"\\_"),"[[",2)
ds$n = as.numeric(gsub("%","",as.character(ds$Var1)))

ggplot(ds, aes(x=n, y=value, group=Var2, col=region)) + geom_line() + geom_point(aes(shape=method))+
  geom_hline(yintercept = cEmp, col="grey", lty='dashed')+
  geom_hline(yintercept = cFisher, col="black", lty='dashed')+
  theme_bw()+
  scale_y_log10()+facet_wrap(as~region)

