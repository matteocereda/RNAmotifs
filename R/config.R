#######################################################################################
## Configuration file by Matteo Cereda <matteo.cereda@kcl.ac.uk>
#######################################################################################


pkgsR  =  c("lattice","latticeExtra","bootstrap")
for (pkgR in pkgsR) 
  if (!pkgR %in% rownames(installed.packages())) { 
    install.packages(pkgR)
    library(pkgR, character.only = TRUE, quietly = TRUE)
  } else { 
    library(pkgR, character.only = TRUE, quietly = TRUE)
  }

colorFun = colorRampPalette(c("blue2","yellow","red2"))
ggrey     = "grey95"
grey.line = rgb(94,90,90,max=255)

options(warn=-1, stringsAsFactors=F)

IUPAC_CODE=c('A','C','T','G','(A/G)','(C/T)','(G/C)','(A/T)','(G/T)','(A/C)','(C/G/T)','(A/G/T)','(A/C/T)','(A/C/G)')
names(IUPAC_CODE) = c('A','C','T','G','R','Y','S','W','K','M','B','D','H','V')

rown = c(950:1200,1800:2050,2950:3200,3800:4050)

print.logo = function(){
  cat(  "                                                                              \n")
  cat( "██████╗  ███╗   ██╗  █████╗  ███╗   ███╗  ██████╗  ████████╗ ██╗ ███████╗ ███████╗    \n")
  cat(  "██╔══██╗ ████╗  ██║ ██╔══██╗ ████╗ ████║ ██╔═══██╗ ╚══██╔══╝ ██║ ██╔════╝ ██╔════╝    \n")
  cat(  "██████╔╝ ██╔██╗ ██║ ███████║ ██╔████╔██║ ██║   ██║    ██║    ██║ █████╗   ███████╗    \n")
  cat(  "██╔══██╗ ██║╚██╗██║ ██╔══██║ ██║╚██╔╝██║ ██║   ██║    ██║    ██║ ██╔══╝   ╚════██║    \n")
  cat(  "██║  ██║ ██║ ╚████║ ██║  ██║ ██║ ╚═╝ ██║ ╚██████╔╝    ██║    ██║ ██║      ███████║    \n")
  cat(  "╚═╝  ╚═╝ ╚═╝  ╚═══╝ ╚═╝  ╚═╝ ╚═╝     ╚═╝  ╚═════╝     ╚═╝    ╚═╝ ╚═╝      ╚══════╝    \n")
  cat(  "                                                                              \n")
}

lapply_pb = function(X, FUN, ...){
  env <- environment()
  pb_Total <- length(X)
  counter <- 0
  pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)   
  
  # wrapper around FUN
  wrapper <- function(...){
    curVal <- get("counter", envir = env)
    assign("counter", curVal +1 ,envir=env)
    setTxtProgressBar(get("pb", envir=env), curVal +1)
    FUN(...)
  }
  res <- lapply(X, wrapper, ...)
  close(pb)
  res
}


get_fisher <- function (dd, exon) {
  outP <- c()
  selE <- dd[,]$type == 1  # 1: Enh, 0: Cont, -1: Sil
  selC <- dd[,]$type == 0  # 1: Enh, 0: Cont, -1: Sil
  selS <- dd[,]$type == -1  # 1: Enh, 0: Cont, -1: Sil
  regList <- list("hits_region1","hits_region2","hits_region3")
  
  for (reg in regList) {
    selP <- dd[exon, reg] > 0    # present, if hits_region1 > 0
    selA <- dd[exon, reg] == 0   # absent, if hits_region1 == 0
    
    cEnhPres <- sum(selE & selP)
    cEnhAbs  <- sum(selE & selA)
    cConPres <- sum(selC & selP)
    cConAbs  <- sum(selC & selA)
    cSilPres <- sum(selS & selP)
    cSilAbs  <- sum(selS & selA)
    
    contEnh <- matrix( data=c(cEnhPres, cConPres, cEnhAbs, cConAbs), nrow=2)
    contSil <- matrix( data=c(cSilPres, cConPres, cSilAbs, cConAbs), nrow=2)
    FEnh <- fisher.test(contEnh, alternative = "greater")
    FSil <- fisher.test(contSil, alternative = "greater")
    
    outP <- c(outP, FEnh$p.value, FSil$p.value)
  }
  outP
}


fisher.region <- function(exon, myHitsData) {
  ## input: named list of matrices with 3 columns named "type","hits_region1","hits_region2","hits_region3" and rows corresponding to exons
  ## output: matrix with 6 cols (r1enh, r1sil, r2enh, r2sil, r3enh, r3sil) and rows corresponding to exons
  numTet <- length(myHitsData)
  namesTet <- names(myHitsData)
  pMatrix <- matrix(ncol=6, nrow=numTet)
  
  outP = lapply_pb(myHitsData,get_fisher, exon=exon)

  for (tIdx in 1:numTet) {#     dd<-myHitsData[[tIdx]]
    pMatrix[tIdx,] <- outP[[tIdx]]
  }
  pAdjMatrix <- matrix(ncol=6, nrow=numTet)
  for (colIdx in 1:6) pAdjMatrix[,colIdx] <- p.adjust(pMatrix[,colIdx], method="BH")
  rownames(pAdjMatrix) <- namesTet
  colnames(pAdjMatrix) <- c("r1enh", "r1sil", "r2enh", "r2sil", "r3enh", "r3sil")
  pAdjMatrix
}

fisher.region.boot <- function(exon, myHitsData) {
  ## input: named list of matrices with 3 columns named "type","hits_region1","hits_region2","hits_region3" and rows corresponding to exons
  ## output: matrix with 6 cols (r1enh, r1sil, r2enh, r2sil, r3enh, r3sil) and rows corresponding to exons
  numTet <- length(myHitsData)
  namesTet <- names(myHitsData)
  pMatrix <- matrix(ncol=6, nrow=numTet)
  cat("=")
  outP = lapply(myHitsData,get_fisher, exon=exon)
  
  for (tIdx in 1:numTet) {#     dd<-myHitsData[[tIdx]]
    pMatrix[tIdx,] <- outP[[tIdx]]
  }
  pAdjMatrix <- matrix(ncol=6, nrow=numTet)
  for (colIdx in 1:6) pAdjMatrix[,colIdx] <- p.adjust(pMatrix[,colIdx], method="BH")
  rownames(pAdjMatrix) <- namesTet
  colnames(pAdjMatrix) <- c("r1enh", "r1sil", "r2enh", "r2sil", "r3enh", "r3sil")
  pAdjMatrix
}


readTab = function(filename){
  read.table(file=filename,header=T,sep="\t",stringsAsFactors=F)
}

writeTab = function(what,filename){
  write.table(what,file=filename, row.names=F,col.names=T,quote=F,sep="\t")	
}

getNApos = function( pos ){
  myget   = function(x){ if(!is.na(x[1])) return(x[2]) else return(NA) }
  mysplit = function(x){ if(!is.na(x[1])) return(unlist(strsplit(x, split="-"))) else return(NA)}
  myna    = function(x){ if(!is.na(x[1])){ if(!is.na(sum(match(x,"NA")))) return(NA) else return(x)} else return(NA)}
  strpin = strsplit(pos,split=":")
  a = lapply(strpin,myget)
  b = lapply( a, mysplit)
  x = lapply(b,myna)
  return( which(is.na(x)) )
}


# ------------------
## FUNCTIONs
# ------------------

getTables = function(protein,place,tets){
  # print(place)
  res = list()
  filelist = read.delim(paste(protein,place,"filelist.txt",sep="",coll=""),header=F,stringsAsFactors=F)[,1] 
  overlap  = pmatch(tets,filelist) 
  if(!(length(overlap)==1 & is.na(overlap))){
    filelist = filelist[na.omit(overlap)]
    sel.tets = tets[which(!is.na(overlap))]
    rown = as.character(c(950:1200,1800:2050,2950:3200,3800:4050))
    tab = matrix( 0, nrow=1004, ncol=length(sel.tets) )
    rownames(tab) = rown 
    colnames(tab) = sel.tets
    res = list("enh"=tab,"sil"=tab,"cont"=tab)
    for (i in 1:length(filelist)){
      if(!is.na(i)){
        df = read.delim(paste(protein,place,filelist[i], sep=""),fill=TRUE,skip=1,header=F)
        if(!is.null(df) && nrow(df)>1){
          colnames(df) = c("pos","cat","no")
          tet = sel.tets[i];
          chosen = subset(df, cat==1); rownames(chosen) = as.character(chosen$pos);	res[["enh"]][,tet] = chosen[rown,"no"]
          chosen = subset(df, cat==(-1)); rownames(chosen) = as.character(chosen$pos);res[["sil"]][,tet] = chosen[rown,"no"]
          chosen = subset(df, cat==0); rownames(chosen) = as.character(chosen$pos);	res[["cont"]][,tet] = chosen[rown,"no"]
        }
        for(i in 1:3) res[[i]][which(is.na(res[[i]]))]=0
      }
      
    }
  }
  res
}


lmb.cluster.fisher.test = function(d, tot.d, contr, tot.contr ){
  pvd = matrix(0,nrow=nrow(d),ncol=ncol(d))
  colnames(pvd) = colnames(d)
  for(j in 1:ncol(pvd))
    for(i in 1:nrow(pvd)){
      test = matrix( data=c( d[i,j], contr[i,j], tot.d - d[i,j] , tot.contr - contr[i,j] ), nrow=2 )
      pvd[i,j] = fisher.test(test, alternative = "greater")$p.value
    }
  pvd
}

idAlignMot = function(m,mall){
  mot = getTrims(m)
  lall = vector("list",length(mall))
  for( i in 1:length(mall)) lall[[i]] = getTrims(mall[i])
  l    = sapply(lall,length)
  hits = sapply( lall, function(x,m) length(na.omit(match(x,m))), m=mot )
  res = hits/l
  which(res>=0.5)
}

getTrims = function(mot){
  sp = unlist(strsplit(mot,split=""))
  nr = sum(match(sp,c("Y","R","W","S")),na.rm=T)
  if(nr==0){
    trms= c(paste(sp[1:3],collapse=""),paste(sp[2:4],collapse=""))
  }else{
    trimers = vector("list",nr*2)
    switch(sp[1],
           Y = { trimers[[1]] = paste("C",sp[2],sp[3],sep="",collapse=""); trimers[[2]] = paste("T",sp[2],sp[3],sep="",collapse="");},
           R = { trimers[[1]] = paste("A",sp[2],sp[3],sep="",collapse=""); trimers[[2]] = paste("G",sp[2],sp[3],sep="",collapse="");},
           W = { trimers[[1]] = paste("T",sp[2],sp[3],sep="",collapse=""); trimers[[2]] = paste("A",sp[2],sp[3],sep="",collapse="");},
           S = { trimers[[1]] = paste("C",sp[2],sp[3],sep="",collapse=""); trimers[[2]] = paste("G",sp[2],sp[3],sep="",collapse="");},
           {paste(sp[1],sp[2],sp[3],sep="",collapse="");}
           
    )
    switch(sp[4],
           Y = { trimers[[3]] = paste(sp[2],sp[3],"C",sep="",collapse=""); trimers[[4]] = paste(sp[2],sp[3],"T",sep="",collapse="");},
           R = { trimers[[3]] = paste(sp[2],sp[3],"A",sep="",collapse=""); trimers[[4]] = paste(sp[2],sp[3],"G",sep="",collapse="");},
           W = { trimers[[3]] = paste(sp[2],sp[3],"T",sep="",collapse=""); trimers[[4]] = paste(sp[2],sp[3],"A",sep="",collapse="");},
           S = { trimers[[3]] = paste(sp[2],sp[3],"C",sep="",collapse=""); trimers[[4]] = paste(sp[2],sp[3],"G",sep="",collapse="");},
           {paste(sp[2],sp[3],sp[4],sep="",collapse="");})
    trms = unlist(trimers)
  }
  trms
}

SortingAndType = function(score,nn){
  ord=type=c()
  ii = 1 
  while(length(nn)!=0){
    if(length(nn)>1){
      m  = nn[1]
      nn = nn[2:length(nn)]
      
      id  = idAlignMot(m,nn)
      lid = length(id)
      add =c()
      
      if(lid == 0){
        add = m
      }
      if(lid == 1){
        add = c(m,nn[id])
        nn  = nn[-id]
      }	 
      if(lid > 1) {
        cc = cp  = c()
        for(i in id){
          tmp = cor.test( score[,m], score[,nn[i]], method="pearson")
          cc  = c( cc, tmp$estimate )
          cp  = c( cp, tmp$p.value )
        }
        
        xx  = cbind(cc,cp)
        xs  = order(xx[,1],xx[,2],decreasing=T)
        add = c( m, nn[id[xs]] )
        nn  = nn[-id]
      }
      
      ord = c(ord,add)
      type=c(type, rep(ii,length(add) ))
      
    }else{
      ord  = c(ord,nn)
      type = c(type,ii)
      nn=NULL
    }
    ii=ii+1
  }
  res=list(ord,type)
}

# ------------
# RNA PLOT
# ------------

panel.rnaplot <-function(x, y, ...,ry=29,rc="red", intensity=1){
  colorFun <- colorRampPalette(c("blue2","yellow","red2"))  # cols<-colorFun(101)[ cut(intensity,breaks=101,label=FALSE) ]
  cols<-colorFun(101)[ findInterval(intensity,seq(0,max(intensity,na.rm=T),by=max(intensity,na.rm=T)/100)) ]
  for(i in 1:length(y)) panel.polygon(c(x[i],x[i]),c(0,y[i]),border=cols[i],cex=0.5)
  panel.rect(-30,0,-10,ry,col=rc)
}

rnaScoreMap=function(data,rcol, ylabels=c("0.5","1"),main="",exon=50,intron=200,gap=10,regions=1000,ylim.redox=F){
  offset = data.frame(	reg   = c("1","2","3","4"),
                       start = c( 0, exon + intron + gap, 2*exon + 2*intron + 2*gap, 3*exon + 3*intron + 3*gap ),
                       plus  = c( exon, exon + 2*intron + gap, 3*exon + 2*intron + 2*gap, 4*intron + 3*exon + 3*gap ),
                       end   = c( exon + intron, 2*exon + 2*intron + gap, 3*exon + 3*intron + 2*gap, 4*exon + 4*intron + 3*gap))
  pos_region = c( (1*regions-exon):(1*regions+intron), (2*regions-intron):(2*regions+exon), (3*regions-exon):(3*regions+intron), (4*regions-intron):(4*regions+exon))
  
  aa = as.character(c( (1*regions-exon),(1*regions+intron/2),(1*regions+intron), (2*regions-intron),(2*regions-intron/2),(2*regions+exon), (3*regions-exon),(3*regions+intron/2),(3*regions+intron), (4*regions-intron),(4*regions-intron/2),(4*regions+exon)))
  
  rownames(data)=NULL
  
  pf = cbind(data,reg=NA,val=NA);
  rownames(pf)=NULL
  pf=as.data.frame(pf)
  pf$reg = floor( ( pf$l + intron ) / regions )
  pf$val = pf$l - ( pf$reg * regions)
  pf = cbind(pf,plot=NA); 
  for (reg in c("1","2","3","4")){ pf$plot[ pf$reg==reg ] = pf$val[ pf$reg==reg ] + offset$plus[ offset$reg==reg ] }
  pf_real=pf[ pf$l %in% pos_region,]
  
  if(ylim.redox == F){
    ylim = c(0,max(pf_real$s))
  }else{
    ylim=c(0,round(max(pf_real$s)/3,digits=1))
    ylabels[2] = as.character(round(as.numeric(ylabels[2])/3,digits=1))
  }
  
  pf_real$rcol = rcol[as.character(pf_real$g)]
  
  p1=xyplot(s~plot|g,data=pf_real,
            intensity  = pf_real$z,
            rcol       = pf_real$rcol,
            hrect      = ylim[2],
            aspect="fill",
            type="l",
            layout=c(1,length(levels(data$g))), 
            ylim=ylim,
            xlim = c(-30, offset$end[4]+10 ),
            subscripts=TRUE,
            strip = F,
            strip.left = strip.custom(horizontal = TRUE,var.name=unique(pf_real$g),bg=ggrey),
            strip.left.lines = 0.5,
            panel = function(x,y,hrect,rcol,intensity,groups,subscripts,...){
              panel.xblocks(x=offset$start[1]:offset$end[4],
                            c(rep("0",50),rep(NA,201),rep(NA,9),rep(NA,201),rep("3",50),
                              rep(NA,9),rep("4",50),rep(NA,201),rep(NA,9),rep(NA,201),rep("7",50)), 
                            col = "gray96",border = FALSE,lwd=0.5)
              panel.xyplot(x,y,...)
              panel.xblocks(x=offset$start[1]:offset$end[4],
                            c(rep(NA,251),rep("0",9),rep(NA,251),rep("0",9),rep(NA,251),rep("0",9),rep(NA,251)), 
                            col = "white",border = FALSE)
              panel.abline(v=c(offset$plus[1],offset$plus[2]+1,offset$plus[3],offset$plus[4]+1,
                               offset$start[1],offset$start[2],offset$start[3],offset$start[4],
                               offset$end[1]+1,offset$end[2]+1,offset$end[3]+1,offset$end[4]+1),lty="dotted",col="black")
              panel.rnaplot(x,y,ry=hrect,rc=rcol[subscripts],intensity=intensity[subscripts],...)
              
            },
            scales = list( x = list( rot = 90, at=reg, tick.number = length(reg), labels=NULL, tck = c(0,1)),
                           y = list( tck = c(0,1),at=c(ylim[2]/2,ylim[2]),tick.number=2,
                                     labels=ylabels,
                                     alternating=2),
                           cex=0.5),
            ylab="",
            xlab="",
            main=main,
            par.strip.text = list(cex=0.8),
            par.settings = list(axis.line     = list(lwd = 0.5,col=grey.line),
                                layout.widths = list(key.ylab.padding=1, strip.left=4),
                                strip.border  = list(col=grey.line,lwd=0.5),
                                plot.line     = list(col=grey.line),
                                add.line      = list(lwd=0.5)
            ),
            legend = list(left = list(fun = draw.colorkey,
                                      args = list(key = list( col = colorFun(100),
                                                              at=0:100,
                                                              tick.number=3, 
                                                              labels=list(labels=c("100% S","50% E,S","100% E"),cex=0.5),
                                                              raster = T,
                                                              width=1,height=0.5,
                                                              space="left"), 
                                                  draw = FALSE)))
            
  )
  
  p1	
}

# ------------
# GATHER RESULTS
# ------------
gather_RNAmotifs_results = function(boot_file, splice_infile, res_path, P.EMP, P.FISH,file_table_out){
  print('Loading results ...')
  load(boot_file)
  mats=read.table(splice_infile, sep=";")
  nr = lapply(list.files(paste0(res_path,"/nr"), "_region_count.tsv", full.names = T),
              read.delim)
  names(nr) = gsub( "_region_count.tsv", "", list.files(paste0(res_path,"/nr"), "_region_count.tsv"))
  
  r = lapply(list.files(paste0(res_path,"/r"), "_region_count.tsv", full.names = T),
             read.delim)
  names(r) = gsub( "_region_count.tsv", "", list.files(paste0(res_path,"/r"), "_region_count.tsv"))
  
  print('Merging results ...')
  
  xx = c(nr, r)
  snr = lapply(xx, function(x) ddply(x, .(type), summarise,
                                     R1=sum(hits_region1),
                                     R2=sum(hits_region2),
                                     R3=sum(hits_region3),
                                     tot_hits_in_regions=sum(hits_region1)+sum(hits_region2)+sum(hits_region3)))
  
  s = mapply(function(x,y){ x$tet=y; return(x)}, x=snr, names(snr), SIMPLIFY = F )
  s = do.call(rbind, s)
  m=s
  t = table(mats$V9)
  
  
  m$n_exons = t[as.character(m$type)]
  
  mt=melt(m, id.vars = c('tet','type','tot_hits_in_regions','n_exons') )
  
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
  print('Printing table ...')
  
  mf = mf[,c('key','tet','type', 'tot_hits_in_regions','variable', 'value', 'n_exons', 'value.ct', 'n_exons.ct', 'pv', 'Emp' )]
  
  write.xlsx(mf, file=file_table_out, row.names = F)
  return(list(mf,s))
}  

draw_RNAmaps = function(ires, P.EMP,P.FISH, config_file, res_path, RNAmaps_out,hh=10){
  e = subset(ires[[2]], type==1)
  c = subset(ires[[2]], type==(0))
  s = subset(ires[[2]], type==(-1))
  
  smf = subset(ires[[1]], pv<=P.FISH & Emp<=P.EMP)
  
  source(config_file)
  pp = res_path
  rown = c(950:1200,1800:2050,2950:3200,3800:4050)
  
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
  
  print("Calculating positional enrichment ...")
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
  # print("Calculating positional enrichment ...")
  # 
  # write.xlsx(smf, file="RNAmotifs/PTEN.human.RNAmotifs.xlsx", "Sig Tets", append=T, row.names = F)
  
  # plot RNA MAPS
  print("Plotting RNA maps...\n")
  
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
  
  pdf(file=RNAmaps_out, height=hh)
  print(prn,panel.height=list(0.4,"cm"),panel.width=list(10,"cm"))
  dev.off()
  prn
}


from_MATS_to_RNAmotifs = function(exons, fname){
  #row_id;second_id;chrom;strand;upstream_exon_end_position;exon_start_position;exon_end_position;dwstream_exon_start_position;diRank
  exons$DI = 0
  exons$DI[exons$type=="Silenced"]=(-1)
  exons$DI[exons$type=="Enhanced"]=1
  cn = c("ID",'chr','strand','upstreamEE','exonStart_0base','exonEnd','downstreamES',"DI")
  res = cbind(1:nrow(exons), exons[,cn])
  write.table(res, file=fname, col.names = F, row.names = F, quote = F, sep=";")
}
