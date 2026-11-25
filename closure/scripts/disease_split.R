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

len <- mat(file.path(f,'import_length.tsv.gz'))
stopifnot(identical(dimnames(tp),dimnames(len)))
W <- tp*rowMeans(len); A <- sweep(W,2,colSums(raw)/colSums(W),'*')
err <- max(abs(A[rownames(s$counts),colnames(s$counts)]-s$counts)); print(err); stopifnot(err<1e-7)
ct <- A[a$gene_id,c(groups$NC,groups$DM)]
ct <- ct[rowSums(ct>=10)>=4,,drop=FALSE]
voomfit<-function(m,lo,hi,onlyV=FALSE){
 m<-m[,c(lo,hi),drop=FALSE];grp<-factor(rep(c(0,1),c(length(lo),length(hi))));X<-model.matrix(~grp)
 y<-calcNormFactors(DGEList(counts=m,group=grp))
 V<-voom(y,X,plot=FALSE);mod<-eBayes(lmFit(V,X),robust=TRUE)
 list(t=mod$t[,2],v=V,design=X,y=y)
}
rel<-list()
for (cmp in 'disease'){
 hi<-if(cmp=='treatment')groups$PNF else groups$DM
 lo<-if(cmp=='treatment')groups$DM else groups$NC
 counts<-if(cmp=='treatment')s$counts else ct
 pairs<-0
 for (h in combn(hi,2,simplify=FALSE)) if(hi[1]%in%h) for(l in combn(lo,2,simplify=FALSE)){
  hh<-setdiff(hi,h);ll<-setdiff(lo,l);pairs<-pairs+1
  corrTP<-rho(delta(glog,h,l),delta(glog,hh,ll))
  corrV<-rho(voomfit(counts,l,h)$t,voomfit(counts,ll,hh)$t)
  rel[[length(rel)+1]]<-data.frame(contrast=cmp,statistic=c('log2(1+TPM) mean difference','voom_t'),scale=c('TPM',if(cmp=='treatment')'lengthScaledTPM' else 'lengthScaledTPM'),n_genes=c(nrow(glog),nrow(counts)),pair=pairs,half1_hi=paste(h,collapse='|'),half1_lo=paste(l,collapse='|'),half2_hi=paste(hh,collapse='|'),half2_lo=paste(ll,collapse='|'),rho=c(corrTP,corrV))
 }
 stopifnot(pairs==18)
}
rl<-do.call(rbind,rel);rl$row_type<-'pair'
summ<-do.call(rbind,lapply(split(rl,interaction(rl$contrast,rl$statistic)),function(x){
 stopifnot(nrow(x)==18,all(is.finite(x$rho)));y<-x[1,];y[,c('pair','half1_hi','half1_lo','half2_hi','half2_lo','rho')]<-NA
 y$row_type<-'summary';y$min<-min(x$rho);y$q25<-unname(quantile(x$rho,.25));y$median<-median(x$rho);y$mean<-mean(x$rho);y$q75<-unname(quantile(x$rho,.75));y$IQR<-IQR(x$rho);y$max<-max(x$rho);y
}))
for(nm in setdiff(names(summ),names(rl)))rl[[nm]]<-NA
write.csv(rbind(rl,summ),file.path(o,'round71b/root_qa/flower_disease_A_split_half.csv'),row.names=FALSE)
