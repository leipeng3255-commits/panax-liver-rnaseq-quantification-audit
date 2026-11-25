o <- normalizePath(commandArgs(trailingOnly=TRUE)[1], mustWork=TRUE)
r <- file.path(o,'root')
f <- 'closure/inputs'
p <- 'models/flower'
library(edgeR);library(limma)
stopifnot(as.character(getRversion())=='4.6.1',as.character(packageVersion('edgeR'))=='4.10.7',as.character(packageVersion('limma'))=='3.68.5')
sessionInfo()
tab <- function(x) read.delim(x,check.names=FALSE)
mat <- function(file){a<-tab(file);stopifnot(!anyDuplicated(a$gene_id));z<-as.matrix(a[,-1]);rownames(z)<-a$gene_id;z}
tp <- mat(file.path(f,'import_abundance.tsv.gz'))
raw <- mat(file.path(f,'import_counts.tsv.gz'))
stopifnot(identical(dimnames(raw),dimnames(tp)))
s <- readRDS(file.path(p,'inputs/filtered_model_inputs.rds'))
j <- read.csv(file.path(f,'flower_samples.csv.gz'))
groups <- split(j$run,j$condition)
groups <- lapply(groups,sort)
a <- tab(file.path(p,'context_inputs/root/gene_annotation.tsv'));rownames(a)<-a$gene_id
g <- read.csv(file.path(f,'flower_descriptive.csv.gz'),check.names=FALSE)
z <- read.csv(file.path(f,'root_descriptive.csv.gz'));rownames(z)<-z$gene_id
# Historical descriptive TPM comes from quant.genes.sf (round37); use the same stored TPM values for both contrasts.
glog <- log2(1+as.matrix(g[,j$run]));rownames(glog)<-g$gene_id
delta <- function(x,hi,lo) rowMeans(x[,hi,drop=FALSE])-rowMeans(x[,lo,drop=FALSE])
rho <- function(a,b) cor(a,b,method='spearman')
rows<-list()
record <- function(comp,stat,x,y,expression,domain,scale){
 stopifnot(length(x)==length(y),all(is.finite(c(x,y,expression))))
 low<-expression<=median(expression)
 for (nm in c('all','lower','higher')){
  k<-switch(nm,all=rep(TRUE,length(x)),lower=low,higher=!low)
  rows[[length(rows)+1]]<<-data.frame(contrast=comp,statistic=stat,stratum=nm,n_genes=sum(k),rho=rho(x[k],y[k]),domain=domain,count_scale=scale)
 }
}
ex <- rowMeans(glog)
dt <- delta(glog,groups$PNF,groups$DM)
stopifnot(max(abs(dt-g$flower_delta))<1e-12)
record('treatment','log2(1+TPM) mean difference',g$root_delta,dt,ex,'round37 common descriptive','TPM')
record('disease','log2(1+TPM) mean difference',z[g$gene_id,'delta_model_control'],delta(glog,groups$DM,groups$NC),ex,'round37 common descriptive','TPM')
rt<-tab(file.path(r,'gene_results_all.tsv'));ft<-tab(file.path(o,'flower/gene_results_all.tsv'))
id<-intersect(rt$gene_id,ft$gene_id);rt<-rt[match(id,rt$gene_id),];ft<-ft[match(id,ft$gene_id),]
expr<-rowMeans(log2(1+tp[id,,drop=FALSE]))
for(stat in c('edgeR_logFC','voom_t'))record('treatment',stat,rt[[stat]],ft[[stat]],expr,'intersection of native coding model domains','historical lengthScaledTPM')
# Additional secondary disease model C: exact scaledTPM formula from cached gene-level abundance + estimated library totals.
C <- sweep(tp,2,colSums(raw)/colSums(tp),'*')
hi <- groups$DM;lo<-groups$NC
ct<-C[a$gene_id,c(lo,hi),drop=FALSE];ct<-ct[rowSums(ct>=10)>=4,,drop=FALSE]
voomfit<-function(m,lo,hi,onlyV=FALSE){
 m<-m[,c(lo,hi),drop=FALSE];grp<-factor(rep(c(0,1),c(length(lo),length(hi))));X<-model.matrix(~grp)
 y<-calcNormFactors(DGEList(counts=m,group=grp))
 V<-voom(y,X,plot=FALSE);mod<-eBayes(lmFit(V,X),robust=TRUE)
 list(t=mod$t[,2],v=V,design=X,y=y)
}
cf<-voomfit(ct,lo,hi)
u<-which(a[rownames(ct),'unique_symbol_in_domain'])
gmt<-strsplit(readLines(file.path(p,'resources/mh.all.v2025.1.Mm.symbols.gmt')),'\t')
idx<-setNames(lapply(gmt,function(q)which(a[rownames(ct)[u],'symbol']%in%q[-(1:2)])),vapply(gmt,'[','',1))
stopifnot(length(idx)==50);eligible<-lengths(idx)>=15 & lengths(idx)<=500
cam<-do.call(rbind,lapply(c('estimated','zero'),function(tag){
 obj<-camera(cf$v[u,],idx[eligible],cf$design,contrast=2,inter.gene.cor=if(tag=='estimated') NA_real_ else 0,allow.neg.cor=FALSE,sort=FALSE)
 obj<-obj[match(names(idx),rownames(obj)),,drop=FALSE]
 data.frame(study='flower',scale='scaledTPM_secondary',correlation=tag,gene_set=names(idx),NGenes=lengths(idx),eligible=eligible,direction=obj$Direction,p=obj$PValue,BH=obj$FDR,significant=obj$FDR<.05,overlap_other_study=NA)
}))
cam$status<-ifelse(cam$eligible,'computed_secondary','below_15_gene_coverage')
write.csv(cam,file.path(o,'disease_positive_control_camera.csv'),row.names=FALSE)
cross<-do.call(rbind,rows)
write.csv(cross,file.path(o,'cross_context_snr_calibration.csv'),row.names=FALSE)
# Each labelled treated pair occurs once if it contains the first ordered run. All 6 choices for comparator group.
rel<-list()
for (cmp in c('treatment','disease')){
 hi<-if(cmp=='treatment')groups$PNF else groups$DM
 lo<-if(cmp=='treatment')groups$DM else groups$NC
 counts<-if(cmp=='treatment')s$counts else ct
 pairs<-0
 for (h in combn(hi,2,simplify=FALSE)) if(hi[1]%in%h) for(l in combn(lo,2,simplify=FALSE)){
  hh<-setdiff(hi,h);ll<-setdiff(lo,l);pairs<-pairs+1
  corrTP<-rho(delta(glog,h,l),delta(glog,hh,ll))
  corrV<-rho(voomfit(counts,l,h)$t,voomfit(counts,ll,hh)$t)
  rel[[length(rel)+1]]<-data.frame(contrast=cmp,statistic=c('log2(1+TPM) mean difference','voom_t'),scale=c('TPM',if(cmp=='treatment')'lengthScaledTPM' else 'scaledTPM_secondary'),n_genes=c(nrow(glog),nrow(counts)),pair=pairs,half1_hi=paste(h,collapse='|'),half1_lo=paste(l,collapse='|'),half2_hi=paste(hh,collapse='|'),half2_lo=paste(ll,collapse='|'),rho=c(corrTP,corrV))
 }
 stopifnot(pairs==18)
}
rl<-do.call(rbind,rel);rl$row_type<-'pair'
summ<-do.call(rbind,lapply(split(rl,interaction(rl$contrast,rl$statistic)),function(x){
 stopifnot(nrow(x)==18,all(is.finite(x$rho)));y<-x[1,];y[,c('pair','half1_hi','half1_lo','half2_hi','half2_lo','rho')]<-NA
 y$row_type<-'summary';y$min<-min(x$rho);y$q25<-unname(quantile(x$rho,.25));y$median<-median(x$rho);y$mean<-mean(x$rho);y$q75<-unname(quantile(x$rho,.75));y$IQR<-IQR(x$rho);y$max<-max(x$rho);y
}))
for(nm in setdiff(names(summ),names(rl)))rl[[nm]]<-NA
write.csv(rbind(rl,summ),file.path(o,'flower_split_half_reliability.csv'),row.names=FALSE)
print(cross);print(summ[,c('contrast','statistic','n_genes','min','median','mean','IQR','max')]);print(tapply(cam$significant,cam$correlation,sum,na.rm=TRUE))
